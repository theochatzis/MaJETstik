"""Jet clustering and substructure, via a just-in-time compiled FastJet kernel.

Why a C++ kernel?
-----------------
FastJet is a C++ library.  The LCG view we build on ships the library, its
headers and the ``fastjet-contrib`` add-ons (N-subjettiness, soft drop) but not
the optional Python bindings.  Rather than add a compile step to the project,
we hand a small C++ source string to ROOT's built-in C++ interpreter (Cling),
which compiles it in memory the first time it is needed.  From then on the
inner loop -- the part that is called once per event -- runs at full C++ speed,
while everything around it stays ordinary Python.

The result is one importable class::

    from majetstik.fastjet import JetClusterer
    clusterer = JetClusterer(algorithm="antikt", R=0.4, pt_min=20.0)
    jets = clusterer.cluster(px, py, pz, E)      # numpy arrays, one event
    print(jets["pt"], jets["tau21"], jets["softdrop_mass"])

Physics recap
-------------
*Sequential recombination.*  All the algorithms here repeatedly merge the pair
of objects with the smallest distance

    d_ij = min(kT_i^{2p}, kT_j^{2p}) * DeltaR_ij^2 / R^2 ,    d_iB = kT_i^{2p}

with ``p = -1`` for anti-kT, ``p = 1`` for kT and ``p = 0`` for Cambridge/Aachen.
When ``d_iB`` is the smallest, object *i* is called a jet and removed.  Anti-kT
(arXiv:0802.1189) is the LHC default because the negative power makes hard
particles cluster first, giving jets with a regular, cone-like boundary that is
insensitive to soft radiation -- which is what makes them calibratable.

*N-subjettiness* tau_N (arXiv:1011.2268) measures how well a jet's radiation is
described by N subjets.  A jet from a single parton has small tau_1; a jet from
a two-pronged decay (W, Z, H) has tau_1 large but tau_2 small.  The ratio
tau_21 = tau_2/tau_1 is therefore the standard two-prong tagger, and
tau_32 = tau_3/tau_2 the three-prong (top) tagger.

*Soft drop* (arXiv:1402.2657) walks the jet's clustering history backwards and
drops the softer branch whenever

    min(pT1, pT2) / (pT1 + pT2)  <  z_cut * (DeltaR_12 / R)^beta .

This removes soft, wide-angle radiation (underlying event, pileup), leaving a
"groomed" mass that is much closer to the mass of the originating particle and
far less sensitive to the detector environment.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Optional

import numpy as np

_KERNEL_LOADED = False
_FASTJET_PREFIX: Optional[str] = None


# ---------------------------------------------------------------------------
# The C++ side.  Kept in one string so the whole kernel is visible in one place.
# ---------------------------------------------------------------------------
_KERNEL_SOURCE = r"""
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include <string>
#include <vector>
#include <cmath>

namespace majetstik {

// Flat, parallel arrays: easy to turn into numpy on the Python side, and it
// matches the column-oriented layout of the NanoAOD output.
struct ClusterResult {
  std::vector<float> pt, eta, phi, mass;
  std::vector<float> softdrop_mass, softdrop_pt;
  std::vector<float> tau1, tau2, tau3;
  std::vector<int>   nconstituents;
  // Jet -> constituent map, stored the way CMS PFNano does it: two parallel
  // arrays, one entry per (jet, constituent) pair.  See docs/output_format.md.
  std::vector<int>   constituent_jet_index;
  std::vector<int>   constituent_input_index;
};

class Clusterer {
 public:
  Clusterer(const std::string& algorithm, double R, double pt_min,
            double eta_max, double nsub_beta, double sd_beta, double sd_zcut)
      : R_(R), pt_min_(pt_min), eta_max_(eta_max),
        nsub_beta_(nsub_beta), sd_beta_(sd_beta), sd_zcut_(sd_zcut),
        jet_def_(resolve(algorithm), R),
        soft_drop_(sd_beta, sd_zcut, R) {}

  // One event in, one event out.  The caller passes four parallel arrays of
  // length n holding the four-momenta of the detector-level particle-flow
  // candidates (or truth particles, for generator-level jets).
  ClusterResult cluster(const double* px, const double* py,
                        const double* pz, const double* E, int n) const {
    std::vector<fastjet::PseudoJet> inputs;
    inputs.reserve(n);
    for (int i = 0; i < n; ++i) {
      fastjet::PseudoJet p(px[i], py[i], pz[i], E[i]);
      // user_index lets us trace every constituent back to the input array
      // after clustering, which is how Jet_constituents is filled.
      p.set_user_index(i);
      inputs.push_back(p);
    }

    ClusterSequencePtr cs(new fastjet::ClusterSequence(inputs, jet_def_));
    std::vector<fastjet::PseudoJet> jets =
        fastjet::sorted_by_pt(cs->inclusive_jets(pt_min_));

    // Axes are found with one pass of kT axes refined by a single
    // minimisation step: the standard, well-tested choice.
    fastjet::contrib::OnePass_KT_Axes axes;
    fastjet::contrib::NormalizedMeasure measure(nsub_beta_, R_);
    fastjet::contrib::Nsubjettiness n1(1, axes, measure);
    fastjet::contrib::Nsubjettiness n2(2, axes, measure);
    fastjet::contrib::Nsubjettiness n3(3, axes, measure);

    ClusterResult out;
    int kept = 0;
    for (const fastjet::PseudoJet& jet : jets) {
      if (eta_max_ > 0 && std::fabs(jet.eta()) > eta_max_) continue;

      out.pt.push_back(jet.pt());
      out.eta.push_back(jet.eta());
      out.phi.push_back(jet.phi_std());          // phi in (-pi, pi]
      out.mass.push_back(jet.m());

      const std::vector<fastjet::PseudoJet> cons = jet.constituents();
      out.nconstituents.push_back(static_cast<int>(cons.size()));
      for (const fastjet::PseudoJet& c : cons) {
        out.constituent_jet_index.push_back(kept);
        out.constituent_input_index.push_back(c.user_index());
      }

      // N-subjettiness is undefined for a jet with fewer constituents than
      // requested axes; FastJet returns 0 there, which we keep as a sentinel.
      out.tau1.push_back(cons.size() >= 1 ? n1(jet) : -1.f);
      out.tau2.push_back(cons.size() >= 2 ? n2(jet) : -1.f);
      out.tau3.push_back(cons.size() >= 3 ? n3(jet) : -1.f);

      const fastjet::PseudoJet groomed = soft_drop_(jet);
      if (groomed != 0) {
        out.softdrop_mass.push_back(groomed.m());
        out.softdrop_pt.push_back(groomed.pt());
      } else {                                   // grooming removed the jet
        out.softdrop_mass.push_back(-1.f);
        out.softdrop_pt.push_back(-1.f);
      }
      ++kept;
    }
    return out;
  }

 private:
  typedef fastjet::SharedPtr<fastjet::ClusterSequence> ClusterSequencePtr;

  static fastjet::JetAlgorithm resolve(const std::string& name) {
    if (name == "antikt")    return fastjet::antikt_algorithm;
    if (name == "kt")        return fastjet::kt_algorithm;
    if (name == "cambridge") return fastjet::cambridge_algorithm;
    return fastjet::antikt_algorithm;
  }

  double R_, pt_min_, eta_max_, nsub_beta_, sd_beta_, sd_zcut_;
  fastjet::JetDefinition jet_def_;
  fastjet::contrib::SoftDrop soft_drop_;
};

inline std::string version() { return fastjet::fastjet_version_string(); }

}  // namespace majetstik
"""


def _fastjet_prefix() -> str:
    """Locate the FastJet installation via ``fastjet-config``."""
    global _FASTJET_PREFIX
    if _FASTJET_PREFIX is not None:
        return _FASTJET_PREFIX
    if not shutil.which("fastjet-config"):
        raise RuntimeError(
            "fastjet-config not found on PATH.\n"
            "Did you run  `source setup_env.sh`  first?"
        )
    _FASTJET_PREFIX = subprocess.check_output(
        ["fastjet-config", "--prefix"], text=True).strip()
    return _FASTJET_PREFIX


def load_kernel(verbose: bool = False):
    """Compile the FastJet kernel into the running process (idempotent).

    The first call takes a second or two; later calls are free.
    """
    global _KERNEL_LOADED
    import ROOT

    if _KERNEL_LOADED:
        return ROOT

    prefix = _fastjet_prefix()
    ROOT.gInterpreter.AddIncludePath(f"-I{prefix}/include")

    # libfastjetcontribfragile carries Nsubjettiness and SoftDrop.
    for lib in ("libfastjet", "libfastjettools", "libfastjetplugins",
                "libfastjetcontribfragile"):
        if ROOT.gSystem.Load(lib) < 0:
            raise RuntimeError(
                f"could not load {lib}.so from the FastJet installation at "
                f"{prefix}. Check that `source setup_env.sh` succeeded."
            )

    if not ROOT.gInterpreter.Declare(_KERNEL_SOURCE):
        raise RuntimeError(
            "failed to JIT-compile the MaJETstik FastJet kernel; "
            "the compiler diagnostics above say why."
        )
    _KERNEL_LOADED = True
    if verbose:
        print(f"[majetstik.fastjet] kernel compiled against FastJet "
              f"{fastjet_version()} at {prefix}")
    return ROOT


def fastjet_version() -> str:
    """Return the FastJet version number, e.g. ``'3.4.3'``.

    FastJet itself reports ``'FastJet version 3.4.3'``; we keep just the number
    so it reads well in logs and in the provenance record.
    """
    ROOT = load_kernel()
    raw = str(ROOT.majetstik.version())
    return raw.rsplit(" ", 1)[-1] if " " in raw else raw


class JetClusterer:
    """Cluster one event at a time and return jets with substructure.

    Parameters
    ----------
    algorithm : {'antikt', 'kt', 'cambridge'}
        Sequential recombination algorithm.  'antikt' is the LHC standard.
    R : float
        Jet radius.  0.4 is the standard "small-radius" choice used for quark
        and gluon jets; 0.8 or 1.0 is used for boosted W/Z/H/top jets.
    pt_min : float
        Minimum jet pT in GeV.
    eta_max : float
        Maximum |eta|; pass 0 to disable the cut.
    nsub_beta, softdrop_beta, softdrop_zcut : float
        Substructure parameters, documented in the module docstring.

    Example
    -------
    >>> import numpy as np                                  # doctest: +SKIP
    >>> c = JetClusterer(R=0.4, pt_min=20.0)                # doctest: +SKIP
    >>> out = c.cluster(px, py, pz, E)                      # doctest: +SKIP
    >>> out["pt"][:3]                                       # doctest: +SKIP
    array([142.3,  98.1,  24.6], dtype=float32)
    """

    #: field name -> numpy dtype of everything :meth:`cluster` returns
    FIELDS = {
        "pt": np.float32, "eta": np.float32, "phi": np.float32,
        "mass": np.float32, "softdrop_mass": np.float32,
        "softdrop_pt": np.float32, "tau1": np.float32, "tau2": np.float32,
        "tau3": np.float32, "nconstituents": np.int32,
        "constituent_jet_index": np.int32, "constituent_input_index": np.int32,
    }

    def __init__(self, algorithm: str = "antikt", R: float = 0.4,
                 pt_min: float = 20.0, eta_max: float = 4.7,
                 nsub_beta: float = 1.0, softdrop_beta: float = 0.0,
                 softdrop_zcut: float = 0.1):
        ROOT = load_kernel()
        self.algorithm = algorithm
        self.R = R
        self.pt_min = pt_min
        self.eta_max = eta_max
        self._impl = ROOT.majetstik.Clusterer(
            algorithm, R, pt_min, eta_max, nsub_beta, softdrop_beta, softdrop_zcut)

    def cluster(self, px, py, pz, E) -> Dict[str, np.ndarray]:
        """Cluster the four-vectors of a single event.

        Parameters
        ----------
        px, py, pz, E : array_like
            Parallel arrays of equal length, in GeV.

        Returns
        -------
        dict of numpy arrays
            Per-jet arrays (``pt``, ``eta``, ``phi``, ``mass``, ``tau1..3``,
            ``softdrop_mass``, ``softdrop_pt``, ``nconstituents``) plus the two
            flat constituent-map arrays ``constituent_jet_index`` and
            ``constituent_input_index``.  Jets are sorted by decreasing pT.
        """
        px = np.ascontiguousarray(px, dtype=np.float64)
        py = np.ascontiguousarray(py, dtype=np.float64)
        pz = np.ascontiguousarray(pz, dtype=np.float64)
        E = np.ascontiguousarray(E, dtype=np.float64)
        n = px.size
        if not (py.size == pz.size == E.size == n):
            raise ValueError("px, py, pz and E must have the same length")

        if n == 0:
            return {k: np.empty(0, dtype=d) for k, d in self.FIELDS.items()}

        res = self._impl.cluster(px, py, pz, E, n)
        return {name: np.asarray(getattr(res, name), dtype=dtype).copy()
                for name, dtype in self.FIELDS.items()}

    def __repr__(self) -> str:
        return (f"JetClusterer(algorithm={self.algorithm!r}, R={self.R}, "
                f"pt_min={self.pt_min}, eta_max={self.eta_max})")


def ptetaphim_to_pxpypze(pt, eta, phi, mass):
    """Convert (pt, eta, phi, m) to (px, py, pz, E).

    Delphes stores objects in the first form; FastJet wants the second.
    Vectorised over numpy arrays.
    """
    pt = np.asarray(pt, dtype=np.float64)
    eta = np.asarray(eta, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    mass = np.asarray(mass, dtype=np.float64)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    E = np.sqrt(px * px + py * py + pz * pz + mass * mass)
    return px, py, pz, E
