#!/usr/bin/env python3
"""Stage 3 -- anti-kT jet clustering with FastJet, plus substructure.

What this stage does
--------------------
Quarks and gluons are never observed directly: confinement turns each one into
a collimated spray of hadrons.  A *jet algorithm* is the operational definition
that turns that spray back into a single object we can compare with theory.

We take the particle-flow candidates Delphes produced -- the detector's best
guess at the list of particles it saw -- and cluster them with the anti-kT
algorithm at R = 0.4, the LHC standard.  For each resulting jet we also compute

* **N-subjettiness** tau1, tau2, tau3: how well the jet is described by 1, 2 or
  3 subjets.  The ratios tau21 = tau2/tau1 and tau32 = tau3/tau2 are the
  standard two-prong (W/Z/H) and three-prong (top) taggers.
* **Soft-drop mass**: the jet mass after recursively removing soft, wide-angle
  radiation.  Much closer to the mass of whatever produced the jet than the
  raw mass, and far less sensitive to the underlying event.

The algorithms themselves are explained in the docstring of
:mod:`majetstik.fastjet`, which also holds the C++ kernel.

Why cluster ourselves when Delphes already makes jets?
------------------------------------------------------
Three reasons, all pedagogical: you see exactly which inputs go in, you control
the algorithm and radius without editing a detector card, and you get
substructure observables that the stock card does not write out.  Stage 4
compares our jets with the Delphes ones by matching them in (eta, phi).

Input
-----
``data/<run_name>_delphes.root`` from stage 2.  We read three collections and
merge them into one list of particle-flow candidates per event:

===================== ============================ ===================
collection            what it is                   four-vector from
===================== ============================ ===================
EFlowTrack            charged particles, from      PT, Eta, Phi, Mass
                      tracks (best at low pT)
EFlowPhoton           photons, from ECAL towers    ET, Eta, Phi (m=0)
EFlowNeutralHadron    neutral hadrons, HCAL towers ET, Eta, Phi (m=0)
===================== ============================ ===================

Output
------
``data/<run_name>_jets.root`` with a ``Jets`` tree holding the clustered jets,
the particle-flow candidates, and the jet->constituent index map.  Stage 4
merges this with the leptons and MET from Delphes into the final NanoAOD file.

Run standalone
--------------
    source setup_env.sh
    python3 analysis/cluster_jets.py --run-name smoke --jet-r 0.8

How to modify
-------------
* **Large-radius jets for boosted objects**: ``--jet-r 0.8``.  Substructure is
  much more informative there, since an R=0.4 jet rarely contains two prongs.
* **A different algorithm**: ``--algorithm kt`` or ``cambridge``.
* **Cluster truth particles instead**, to study the effect of the detector:
  point ``read_pf_candidates`` at the ``Particle`` collection with status 1.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Dict

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from majetstik.config import PipelineConfig
from majetstik.fastjet import JetClusterer, ptetaphim_to_pxpypze, fastjet_version
from majetstik.logging_setup import banner, get_logger
from majetstik.provenance import build_provenance, merge_provenance

#: PDG-like identifiers we assign to particle-flow candidates.  Delphes does
#: not tell us the true species of a neutral tower, so these are *categories*:
#: 22 = photon, 130 = K0_L standing in for "neutral hadron", and for charged
#: candidates we keep the PID Delphes assigned to the track.
PF_PHOTON_ID = 22
PF_NEUTRAL_HADRON_ID = 130


def read_pf_candidates(tree, entry_stop=None) -> Dict[str, "awkward.Array"]:
    """Read and merge the three particle-flow collections from a Delphes tree.

    Returns a dict of jagged (per-event) awkward arrays: ``pt``, ``eta``,
    ``phi``, ``mass``, ``charge``, ``pdgId``.  The concatenation order --
    tracks, then photons, then neutral hadrons -- defines the candidate index
    used by the jet constituent map, so it must stay fixed.
    """
    import awkward as ak

    branches = [
        "EFlowTrack.PT", "EFlowTrack.Eta", "EFlowTrack.Phi",
        "EFlowTrack.Mass", "EFlowTrack.Charge", "EFlowTrack.PID",
        "EFlowPhoton.ET", "EFlowPhoton.Eta", "EFlowPhoton.Phi",
        "EFlowNeutralHadron.ET", "EFlowNeutralHadron.Eta",
        "EFlowNeutralHadron.Phi",
    ]
    arrays = tree.arrays(branches, entry_stop=entry_stop)

    trk_pt = arrays["EFlowTrack.PT"]
    pho_pt = arrays["EFlowPhoton.ET"]
    nh_pt = arrays["EFlowNeutralHadron.ET"]

    # Calorimeter towers carry no mass information; treating them as massless
    # is the standard convention and is what Delphes itself does internally.
    zeros_like = lambda a: ak.zeros_like(a)

    return {
        "pt": ak.concatenate([trk_pt, pho_pt, nh_pt], axis=1),
        "eta": ak.concatenate([arrays["EFlowTrack.Eta"],
                               arrays["EFlowPhoton.Eta"],
                               arrays["EFlowNeutralHadron.Eta"]], axis=1),
        "phi": ak.concatenate([arrays["EFlowTrack.Phi"],
                               arrays["EFlowPhoton.Phi"],
                               arrays["EFlowNeutralHadron.Phi"]], axis=1),
        "mass": ak.concatenate([arrays["EFlowTrack.Mass"],
                                zeros_like(pho_pt), zeros_like(nh_pt)], axis=1),
        "charge": ak.concatenate([arrays["EFlowTrack.Charge"],
                                  zeros_like(pho_pt), zeros_like(nh_pt)], axis=1),
        "pdgId": ak.concatenate([
            arrays["EFlowTrack.PID"],
            ak.full_like(pho_pt, PF_PHOTON_ID),
            ak.full_like(nh_pt, PF_NEUTRAL_HADRON_ID)], axis=1),
    }


#: Neutrinos are invisible to any detector, so they are excluded from
#: generator-level jets -- otherwise the truth jet would contain energy the
#: reconstructed jet could never have seen.
NEUTRINO_IDS = (12, 14, 16)


def read_truth_particles(tree, entry_stop=None):
    """Read the stable generator particles that a detector could see.

    Returns jagged arrays of ``pt``, ``eta``, ``phi``, ``mass`` for particles
    with HepMC status 1 (final state), excluding neutrinos.

    These are clustered with *the same* jet definition as the detector-level
    candidates, which is what makes the generator-vs-reconstructed comparisons
    in stage 5 meaningful: matching R = 0.8 jets against R = 0.4 truth jets
    would bias both the energy response and the efficiency.
    """
    import awkward as ak

    arrays = tree.arrays(["Particle.PT", "Particle.Eta", "Particle.Phi",
                          "Particle.Mass", "Particle.PID", "Particle.Status"],
                         entry_stop=entry_stop)
    visible = (arrays["Particle.Status"] == 1)
    for nu in NEUTRINO_IDS:
        visible = visible & (abs(arrays["Particle.PID"]) != nu)
    return {
        "pt": arrays["Particle.PT"][visible],
        "eta": arrays["Particle.Eta"][visible],
        "phi": arrays["Particle.Phi"][visible],
        "mass": arrays["Particle.Mass"][visible],
    }


def cluster(cfg: PipelineConfig, log=None) -> dict:
    """Cluster every event in the Delphes file and write ``<run>_jets.root``."""
    import awkward as ak
    import uproot

    log = log or get_logger("cluster_jets", cfg.log_dir)
    banner(log, f"STAGE 3/5  Jet clustering  ({cfg.jet_algorithm}, R={cfg.jet_r})")

    if not cfg.delphes_path.is_file():
        raise FileNotFoundError(
            f"{cfg.delphes_path} not found -- run stage 2 first:\n"
            f"    python3 simulation/run_delphes.py --run-name {cfg.run_name}")

    clusterer = JetClusterer(
        algorithm=cfg.jet_algorithm, R=cfg.jet_r, pt_min=cfg.jet_pt_min,
        eta_max=cfg.jet_eta_max, nsub_beta=cfg.nsub_beta,
        softdrop_beta=cfg.softdrop_beta, softdrop_zcut=cfg.softdrop_zcut)

    log.info("fastjet      : %s", fastjet_version())
    log.info("clusterer    : %r", clusterer)
    log.info("substructure : nsub_beta=%.2f  softdrop(beta=%.2f, zcut=%.2f)",
             cfg.nsub_beta, cfg.softdrop_beta, cfg.softdrop_zcut)
    log.info("input        : %s", cfg.delphes_path)

    t0 = time.time()
    with uproot.open(cfg.delphes_path) as f:
        tree = f["Delphes"]
        n_events = tree.num_entries
        log.info("events       : %d", n_events)
        pf = read_pf_candidates(tree)
        truth = read_truth_particles(tree)

    def flat_four_vectors(collection):
        """Flatten a jagged collection into contiguous arrays + event offsets."""
        px, py, pz, energy = ptetaphim_to_pxpypze(
            ak.flatten(collection["pt"]), ak.flatten(collection["eta"]),
            ak.flatten(collection["phi"]), ak.flatten(collection["mass"]))
        counts = np.asarray(ak.num(collection["pt"]), dtype=np.int64)
        return (px, py, pz, energy), np.concatenate([[0], np.cumsum(counts)])

    pf_vectors, pf_offsets = flat_four_vectors(pf)
    truth_vectors, truth_offsets = flat_four_vectors(truth)

    # Accumulate per-event results, then build jagged arrays in one go.
    per_jet: Dict[str, list] = {k: [] for k in (
        "pt", "eta", "phi", "mass", "softdrop_mass", "softdrop_pt",
        "tau1", "tau2", "tau3", "nconstituents")}
    cons_jet, cons_pf = [], []
    gen_jet: Dict[str, list] = {k: [] for k in ("pt", "eta", "phi", "mass")}

    report_every = max(1, n_events // 10)
    for i in range(n_events):
        lo, hi = pf_offsets[i], pf_offsets[i + 1]
        out = clusterer.cluster(*(v[lo:hi] for v in pf_vectors))
        for key in per_jet:
            per_jet[key].append(out[key])
        cons_jet.append(out["constituent_jet_index"])
        cons_pf.append(out["constituent_input_index"])

        # The same jet definition applied to the truth particles.
        glo, ghi = truth_offsets[i], truth_offsets[i + 1]
        gout = clusterer.cluster(*(v[glo:ghi] for v in truth_vectors))
        for key in gen_jet:
            gen_jet[key].append(gout[key])

        if (i + 1) % report_every == 0:
            log.info("  clustered %6d / %d events", i + 1, n_events)

    elapsed = time.time() - t0
    jets = {f"Jet_{k}": ak.Array(v) for k, v in per_jet.items()}
    n_jets = int(ak.sum(ak.num(jets["Jet_pt"])))

    gen_jets = {f"GenJet_{k}": ak.Array(v) for k, v in gen_jet.items()}
    n_gen_jets = int(ak.sum(ak.num(gen_jets["GenJet_pt"])))

    data = {
        **jets,
        **gen_jets,
        "JetConstituent_jetIdx": ak.Array(cons_jet),
        "JetConstituent_pfIdx": ak.Array(cons_pf),
        "PFCand_pt": pf["pt"], "PFCand_eta": pf["eta"], "PFCand_phi": pf["phi"],
        "PFCand_mass": pf["mass"], "PFCand_charge": pf["charge"],
        "PFCand_pdgId": pf["pdgId"],
    }
    if not cfg.store_constituents:
        data.pop("JetConstituent_jetIdx")
        data.pop("JetConstituent_pfIdx")

    with uproot.recreate(cfg.jets_path) as fout:
        fout["Jets"] = data

    mean_jets = n_jets / max(n_events, 1)
    mean_pf = float(ak.mean(ak.num(pf["pt"]))) if n_events else 0.0
    log.info("clustered %d jets in %d events (%.2f jets/event) in %.1f s",
             n_jets, n_events, mean_jets, elapsed)
    log.info("mean PF candidates/event : %.1f", mean_pf)
    log.info("generator jets (same jet definition, truth particles): %d (%.2f/event)",
             n_gen_jets, n_gen_jets / max(n_events, 1))
    if n_jets:
        pts = ak.flatten(jets["Jet_pt"])
        log.info("jet pT: min %.1f  median %.1f  max %.1f GeV",
                 ak.min(pts), np.median(np.asarray(pts)), ak.max(pts))
    else:
        log.warning("no jets passed pt > %.1f GeV -- lower --jet-pt-min or "
                    "use a harder process (generators/cards/dijet.cmnd)",
                    cfg.jet_pt_min)
    log.info("output       : %s (%.1f MB)",
             cfg.jets_path, cfg.jets_path.stat().st_size / 1e6)

    result = {
        "n_events": int(n_events),
        "n_jets": n_jets,
        "n_gen_jets": n_gen_jets,
        "jets_per_event": round(mean_jets, 3),
        "mean_pf_candidates": round(mean_pf, 2),
        "jets_file": str(cfg.jets_path),
        "fastjet_version": fastjet_version(),
        "elapsed_s": round(elapsed, 2),
    }
    merge_provenance(cfg.provenance_path,
                     build_provenance(cfg, "cluster_jets", **result))
    return result


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Stage 3: cluster Delphes particle-flow candidates into jets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--config", help="JSON config file")
    ap.add_argument("--run-name", help="label for input/output files")
    ap.add_argument("--algorithm", choices=["antikt", "kt", "cambridge"],
                    help="jet algorithm")
    ap.add_argument("--jet-r", type=float, help="jet radius R")
    ap.add_argument("--jet-pt-min", type=float, help="minimum jet pT in GeV")
    args = ap.parse_args(argv)

    cfg = PipelineConfig.from_file(args.config) if args.config else PipelineConfig()
    if args.run_name:   cfg.run_name = args.run_name
    if args.algorithm:  cfg.jet_algorithm = args.algorithm
    if args.jet_r:      cfg.jet_r = args.jet_r
    if args.jet_pt_min is not None: cfg.jet_pt_min = args.jet_pt_min
    cfg.validate()
    cfg.ensure_dirs()

    cluster(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
