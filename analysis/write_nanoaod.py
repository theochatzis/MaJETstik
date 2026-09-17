#!/usr/bin/env python3
"""Stage 4 -- assemble a flat, NanoAOD-like ROOT file.

What NanoAOD is, and why we imitate it
--------------------------------------
CMS's analysis formats form a ladder.  RAW and AOD hold full detector
information and C++ objects with methods; reading them requires the experiment's
software framework.  **NanoAOD** is the bottom rung: one flat ``Events`` tree in
which every branch is a plain number or a variable-length array of plain
numbers.  No custom classes, no dictionaries, no framework -- ``uproot``,
``ROOT``, ``RDataFrame`` or even a C++ program with nothing but ROOT can read
it, and a typical event costs only ~1-2 kB.

The naming convention is ``<Collection>_<variable>``: all branches beginning
``Jet_`` have the same length in a given event (the number of jets), all
branches beginning ``Electron_`` have the number of electrons, and so on.  That
is the whole format.  This file follows the same convention so that anything
you learn here transfers directly to real CMS open data.

What this stage does
--------------------
It merges the two halves of the pipeline into that single flat tree:

* from ``<run>_jets.root`` (stage 3): our anti-kT jets, their substructure, the
  particle-flow candidates and the jet -> constituent map;
* from ``<run>_delphes.root`` (stage 2): electrons, muons, missing transverse
  energy, generator-level jets and the truth particle record.

It also does two pieces of *matching*, each in (eta, phi) using

    DeltaR = sqrt(Delta_eta^2 + Delta_phi^2)

1. **Our jets -> Delphes jets**, to borrow the b-tag decision and the parton
   flavour, which Delphes computes but FastJet knows nothing about.
2. **Our jets -> generator jets**, so a reconstruction efficiency and a jet
   energy response can be measured in stage 5.  Both collections were
   clustered in stage 3 with the same algorithm and radius, which is what
   makes that comparison meaningful at any jet_r.

A third, simpler association flags jets that are really an isolated electron or
muon: particle flow puts leptons into the candidate list, so a Z -> ee event
genuinely produces two "jets" that are electrons.  ``Jet_isLepton`` lets an
analysis remove them, which is what a real experiment's overlap removal does.

Input / Output
--------------
in : ``data/<run_name>_jets.root``, ``data/<run_name>_delphes.root``
out: ``data/<run_name>_nano.root`` -- tree ``Events``, plus ``ReadMe`` and
     ``Metadata`` TObjStrings describing every branch and the full provenance.

Run standalone
--------------
    source setup_env.sh
    python3 analysis/write_nanoaod.py --run-name smoke

How to modify
-------------
* **Add a branch**: compute the array, add it to ``columns`` in
  :func:`build_events`, and add a line to :data:`BRANCH_DOC` so the embedded
  ReadMe stays honest.
* **Shrink the file**: set ``store_truth_particles: false`` or raise
  ``truth_pt_min`` in the config; the truth record dominates the size.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from majetstik.config import PipelineConfig
from majetstik.logging_setup import banner, get_logger
from majetstik.provenance import build_provenance, merge_provenance

# ---------------------------------------------------------------------------
# Branch documentation.  This single dict is the source of truth for the
# ReadMe embedded in the ROOT file AND for docs/output_format.md, so the
# documentation cannot drift away from the data.
# ---------------------------------------------------------------------------
BRANCH_DOC: Dict[str, str] = {
    # ---- event level ----
    "run": "Run number. Always 1; present so the file looks like real NanoAOD.",
    "event": "Event number, 0-based, in generation order.",
    "nJet": "Number of reconstructed jets in this event.",
    "nElectron": "Number of reconstructed electrons.",
    "nMuon": "Number of reconstructed muons.",
    "nGenJet": "Number of generator-level (truth) jets.",
    "nParticle": "Number of stored truth particles.",
    "nPFCand": "Number of particle-flow candidates.",

    # ---- jets ----
    "Jet_pt": "Jet transverse momentum [GeV]. Anti-kT, R as recorded in Metadata.",
    "Jet_eta": "Jet pseudorapidity eta = -ln(tan(theta/2)).",
    "Jet_phi": "Jet azimuthal angle [rad], in (-pi, pi].",
    "Jet_mass": "Jet invariant mass [GeV], from the sum of constituent four-vectors.",
    "Jet_nConstituents": "Number of particle-flow candidates clustered into the jet.",
    "Jet_btag": ("1 if the matched Delphes jet is b-tagged, else 0; -1 if no "
                 "Delphes jet was found within DeltaR < 0.2. Delphes applies a "
                 "parameterised efficiency (~70% for true b jets, ~1% for light)."),
    "Jet_flavor": ("PDG id of the parton the matched Delphes jet came from "
                   "(5 = b, 4 = c, 1-3 = light, 21 = gluon, 0 = undefined); "
                   "-1 if unmatched."),
    "Jet_tau1": ("1-subjettiness: how well the jet is described by a single "
                 "subjet. Small for an ordinary quark or gluon jet. "
                 "-1 if undefined (too few constituents)."),
    "Jet_tau2": "2-subjettiness. tau2/tau1 is the standard two-prong (W/Z/H) tagger.",
    "Jet_tau3": "3-subjettiness. tau3/tau2 is the standard three-prong (top) tagger.",
    "Jet_softdrop_mass": ("Jet mass [GeV] after soft-drop grooming, which "
                          "recursively removes soft wide-angle radiation. "
                          "-1 if grooming rejected the jet."),
    "Jet_softdrop_pt": "Transverse momentum [GeV] of the groomed jet.",
    "Jet_isLepton": ("1 if an isolated reconstructed electron or muon lies "
                     "within DeltaR < 0.2 of the jet axis, i.e. this 'jet' is "
                     "really a lepton that particle flow also handed to the "
                     "jet algorithm. Require Jet_isLepton == 0 for hadronic jets."),
    "Jet_genJetIdx": ("Index into the GenJet_* arrays of the closest generator "
                      "jet within DeltaR < 0.2, or -1 if there is none. Use it "
                      "to measure jet energy response and reconstruction efficiency."),

    # ---- jet constituents ----
    "JetConstituent_jetIdx": ("Flat jet->constituent map, part 1: for each "
                              "(jet, constituent) pair, the index into Jet_*. "
                              "Two flat arrays are used instead of a nested "
                              "array because a NanoAOD tree stores only "
                              "one level of variable length."),
    "JetConstituent_pfIdx": ("Flat jet->constituent map, part 2: the index "
                             "into the PFCand_* arrays of that constituent. "
                             "Example: the constituents of jet j are "
                             "PFCand_pt[JetConstituent_pfIdx[JetConstituent_jetIdx == j]]."),

    # ---- particle-flow candidates ----
    "PFCand_pt": "Particle-flow candidate transverse momentum [GeV].",
    "PFCand_eta": "Particle-flow candidate pseudorapidity.",
    "PFCand_phi": "Particle-flow candidate azimuthal angle [rad].",
    "PFCand_mass": "Particle-flow candidate mass [GeV]; 0 for calorimeter towers.",
    "PFCand_charge": "Electric charge in units of e; 0 for neutral candidates.",
    "PFCand_pdgId": ("Species of the candidate: the track's PDG id for charged "
                     "candidates, 22 for photon towers, 130 for neutral-hadron "
                     "towers. Neutral towers carry no true species information."),

    # ---- leptons ----
    "Electron_pt": "Electron transverse momentum [GeV], after detector smearing.",
    "Electron_eta": "Electron pseudorapidity.",
    "Electron_phi": "Electron azimuthal angle [rad].",
    "Electron_charge": "Electron charge in units of e (-1 for e-, +1 for e+).",
    "Electron_iso": ("Relative isolation: scalar sum of the pT of other "
                     "particles in a cone around the electron, divided by the "
                     "electron pT. Prompt electrons from a Z have small values "
                     "(< 0.1); electrons inside jets have large ones."),
    "Muon_pt": "Muon transverse momentum [GeV], after detector smearing.",
    "Muon_eta": "Muon pseudorapidity.",
    "Muon_phi": "Muon azimuthal angle [rad].",
    "Muon_charge": "Muon charge in units of e.",
    "Muon_iso": "Relative isolation, defined as for Electron_iso.",

    # ---- MET ----
    "MET_pt": ("Missing transverse momentum [GeV]: the magnitude of the "
               "negative vector sum of all reconstructed objects. Large when "
               "neutrinos (or new invisible particles) escape."),
    "MET_phi": "Azimuthal direction [rad] of the missing transverse momentum.",
    "GenMET_pt": "Truth-level missing transverse momentum [GeV], from neutrinos only.",
    "GenMET_phi": "Truth-level missing transverse momentum direction [rad].",

    # ---- generator level ----
    "GenJet_pt": ("Generator-level jet pT [GeV]. Clustered from the stable "
                  "truth particles (status 1, neutrinos removed) with exactly "
                  "the same jet definition as Jet_*, so the two may be "
                  "compared directly whatever the jet radius is set to."),
    "GenJet_eta": "Generator-level jet pseudorapidity.",
    "GenJet_phi": "Generator-level jet azimuthal angle [rad].",
    "GenJet_mass": "Generator-level jet mass [GeV].",
    "Particle_pt": "Truth particle transverse momentum [GeV].",
    "Particle_eta": "Truth particle pseudorapidity.",
    "Particle_phi": "Truth particle azimuthal angle [rad].",
    "Particle_mass": "Truth particle mass [GeV].",
    "Particle_pdgId": ("PDG identification code (11 = e-, 13 = mu-, 22 = gamma, "
                       "211 = pi+, 23 = Z, 2212 = proton, 21 = gluon, ...)."),
    "Particle_status": ("HepMC status code: 1 = stable final-state particle, "
                        "2 = decayed, 4 = beam particle, 21-29 = hard process."),
}


def delta_phi(a, b):
    """Signed azimuthal difference wrapped into (-pi, pi]."""
    return (np.asarray(a) - np.asarray(b) + np.pi) % (2 * np.pi) - np.pi


def delta_r(eta1, phi1, eta2, phi2):
    """DeltaR = sqrt(Delta_eta^2 + Delta_phi^2), the standard angular distance.

    It is used everywhere in collider physics because it is approximately
    invariant under boosts along the beam axis, which is the one direction the
    collision's rest frame is unknown in.
    """
    return np.hypot(np.asarray(eta1) - np.asarray(eta2), delta_phi(phi1, phi2))


def match_by_dr(eta1, phi1, eta2, phi2, max_dr: float) -> np.ndarray:
    """For every object in list 1, the index of the closest object in list 2.

    Returns -1 where nothing lies within ``max_dr``.  Written as a small dense
    loop because the collections here hold at most a few tens of objects.
    """
    eta1 = np.asarray(eta1, dtype=np.float64)
    out = np.full(len(eta1), -1, dtype=np.int32)
    if len(eta1) == 0 or len(np.asarray(eta2)) == 0:
        return out
    for i in range(len(eta1)):
        dr = delta_r(eta1[i], phi1[i], eta2, phi2)
        j = int(np.argmin(dr))
        if dr[j] < max_dr:
            out[i] = j
    return out


def _take(values, idx, default: int = -1) -> np.ndarray:
    """Gather ``values[idx]``, substituting ``default`` wherever ``idx`` is -1.

    Guards the case of an empty source collection: an event can contain jets of
    ours but no Delphes jet at all, and plain fancy indexing would then fail.
    """
    values = np.asarray(values)
    out = np.full(len(idx), default, dtype=np.int32)
    good = idx >= 0
    if good.any() and len(values):
        out[good] = values[idx[good]].astype(np.int32)
    return out


def _truth_mask(status, pt, cfg: PipelineConfig):
    """Which truth particles to keep, per ``cfg.truth_status_filter``."""
    status = np.asarray(status)
    keep_pt = np.asarray(pt) > cfg.truth_pt_min
    if cfg.truth_status_filter == "all":
        return np.ones(len(status), dtype=bool)
    if cfg.truth_status_filter == "final":
        return (status == 1) & keep_pt
    # 'final_or_hard': stable particles plus the hard-process record (21-29),
    # which is what you need to ask "which parton did this jet come from?".
    return ((status == 1) & keep_pt) | ((status >= 21) & (status <= 29))


def build_events(cfg: PipelineConfig, log) -> Dict[str, object]:
    """Read both inputs and return the dict of columns for the Events tree."""
    import awkward as ak
    import uproot

    with uproot.open(cfg.jets_path) as fj:
        jets = fj["Jets"].arrays()
    with uproot.open(cfg.delphes_path) as fd:
        tree = fd["Delphes"]
        delphes = tree.arrays([
            "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.BTag", "Jet.Flavor",
            "Electron.PT", "Electron.Eta", "Electron.Phi", "Electron.Charge",
            "Electron.IsolationVar",
            "Muon.PT", "Muon.Eta", "Muon.Phi", "Muon.Charge", "Muon.IsolationVar",
            "MissingET.MET", "MissingET.Phi",
            "GenMissingET.MET", "GenMissingET.Phi",
            "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass",
            "Particle.PID", "Particle.Status",
        ])

    n_events = len(jets["Jet_pt"])
    if len(delphes["MissingET.MET"]) != n_events:
        raise RuntimeError(
            f"event count mismatch: {cfg.jets_path} has {n_events} events but "
            f"{cfg.delphes_path} has {len(delphes['MissingET.MET'])}. "
            "Re-run stages 2 and 3 with the same --run-name.")

    # ---- per-event matching -------------------------------------------
    btag, flavor, is_lepton, genjet_idx = [], [], [], []
    n_matched_delphes = n_matched_gen = n_lepton_jets = 0

    for i in range(n_events):
        jeta = np.asarray(jets["Jet_eta"][i], dtype=np.float64)
        jphi = np.asarray(jets["Jet_phi"][i], dtype=np.float64)

        # 1. our jets -> Delphes jets, to inherit b-tag and parton flavour.
        #    Half the jet radius is the usual matching cone: close enough that
        #    it is the same jet, loose enough to absorb calibration shifts.
        idx = match_by_dr(jeta, jphi, delphes["Jet.Eta"][i], delphes["Jet.Phi"][i],
                          max_dr=0.5 * cfg.jet_r)
        btag.append(_take(delphes["Jet.BTag"][i], idx))
        flavor.append(_take(delphes["Jet.Flavor"][i], idx))
        n_matched_delphes += int((idx >= 0).sum())

        # 2. our jets -> generator jets, for efficiency and energy response.
        #    Both were clustered in stage 3 with the *same* jet definition, so
        #    this comparison stays valid whatever jet_r is set to.
        gidx = match_by_dr(jeta, jphi, jets["GenJet_eta"][i],
                           jets["GenJet_phi"][i], max_dr=0.5 * cfg.jet_r)
        genjet_idx.append(gidx)
        n_matched_gen += int((gidx >= 0).sum())

        # 3. flag jets that are really an isolated lepton.
        lep_eta = np.concatenate([np.asarray(delphes["Electron.Eta"][i]),
                                  np.asarray(delphes["Muon.Eta"][i])])
        lep_phi = np.concatenate([np.asarray(delphes["Electron.Phi"][i]),
                                  np.asarray(delphes["Muon.Phi"][i])])
        lidx = match_by_dr(jeta, jphi, lep_eta, lep_phi, max_dr=0.5 * cfg.jet_r)
        flag = (lidx >= 0).astype(np.int32)
        is_lepton.append(flag)
        n_lepton_jets += int(flag.sum())

    log.info("matched %d/%d jets to a Delphes jet (b-tag, flavour)",
             n_matched_delphes, int(ak.sum(ak.num(jets["Jet_pt"]))))
    log.info("matched %d jets to a generator jet", n_matched_gen)
    log.info("flagged %d jets as isolated leptons (Jet_isLepton == 1)", n_lepton_jets)

    # ---- assemble the flat tree ---------------------------------------
    columns: Dict[str, object] = {
        "run": np.ones(n_events, dtype=np.int32),
        "event": np.arange(n_events, dtype=np.int64),

        "Jet_pt": jets["Jet_pt"],
        "Jet_eta": jets["Jet_eta"],
        "Jet_phi": jets["Jet_phi"],
        "Jet_mass": jets["Jet_mass"],
        "Jet_nConstituents": jets["Jet_nconstituents"],
        "Jet_tau1": jets["Jet_tau1"],
        "Jet_tau2": jets["Jet_tau2"],
        "Jet_tau3": jets["Jet_tau3"],
        "Jet_softdrop_mass": jets["Jet_softdrop_mass"],
        "Jet_softdrop_pt": jets["Jet_softdrop_pt"],
        "Jet_btag": ak.Array(btag),
        "Jet_flavor": ak.Array(flavor),
        "Jet_isLepton": ak.Array(is_lepton),
        "Jet_genJetIdx": ak.Array(genjet_idx),

        "Electron_pt": delphes["Electron.PT"],
        "Electron_eta": delphes["Electron.Eta"],
        "Electron_phi": delphes["Electron.Phi"],
        "Electron_charge": delphes["Electron.Charge"],
        "Electron_iso": delphes["Electron.IsolationVar"],

        "Muon_pt": delphes["Muon.PT"],
        "Muon_eta": delphes["Muon.Eta"],
        "Muon_phi": delphes["Muon.Phi"],
        "Muon_charge": delphes["Muon.Charge"],
        "Muon_iso": delphes["Muon.IsolationVar"],

        # MissingET is a one-element collection in Delphes; flatten it to a
        # scalar per event, which is how NanoAOD stores MET.
        "MET_pt": ak.to_numpy(ak.firsts(delphes["MissingET.MET"], axis=1)),
        "MET_phi": ak.to_numpy(ak.firsts(delphes["MissingET.Phi"], axis=1)),
        "GenMET_pt": ak.to_numpy(ak.firsts(delphes["GenMissingET.MET"], axis=1)),
        "GenMET_phi": ak.to_numpy(ak.firsts(delphes["GenMissingET.Phi"], axis=1)),

        "GenJet_pt": jets["GenJet_pt"],
        "GenJet_eta": jets["GenJet_eta"],
        "GenJet_phi": jets["GenJet_phi"],
        "GenJet_mass": jets["GenJet_mass"],

        "PFCand_pt": jets["PFCand_pt"],
        "PFCand_eta": jets["PFCand_eta"],
        "PFCand_phi": jets["PFCand_phi"],
        "PFCand_mass": jets["PFCand_mass"],
        "PFCand_charge": jets["PFCand_charge"],
        "PFCand_pdgId": jets["PFCand_pdgId"],
    }

    if cfg.store_constituents and "JetConstituent_jetIdx" in jets.fields:
        columns["JetConstituent_jetIdx"] = jets["JetConstituent_jetIdx"]
        columns["JetConstituent_pfIdx"] = jets["JetConstituent_pfIdx"]

    if cfg.store_truth_particles:
        mask = [_truth_mask(delphes["Particle.Status"][i],
                            delphes["Particle.PT"][i], cfg)
                for i in range(n_events)]
        mask = ak.Array(mask)
        columns.update({
            "Particle_pt": delphes["Particle.PT"][mask],
            "Particle_eta": delphes["Particle.Eta"][mask],
            "Particle_phi": delphes["Particle.Phi"][mask],
            "Particle_mass": delphes["Particle.Mass"][mask],
            "Particle_pdgId": delphes["Particle.PID"][mask],
            "Particle_status": delphes["Particle.Status"][mask],
        })

    # NanoAOD stores an explicit multiplicity next to every collection so that
    # a reader can size its buffers without decompressing the arrays.
    for prefix, ref in (("nJet", "Jet_pt"), ("nElectron", "Electron_pt"),
                        ("nMuon", "Muon_pt"), ("nGenJet", "GenJet_pt"),
                        ("nPFCand", "PFCand_pt"), ("nParticle", "Particle_pt")):
        if ref in columns:
            columns[prefix] = ak.to_numpy(ak.num(columns[ref])).astype(np.int32)

    return columns


def readme_text(cfg: PipelineConfig, columns) -> str:
    """The human-readable branch guide embedded in the ROOT file."""
    lines = [
        "MaJETstik NanoAOD-like output",
        "=" * 70,
        "",
        "One flat tree, 'Events', with one entry per collision event.",
        "Branches are named <Collection>_<variable>; all branches sharing a",
        "prefix have the same length within an event (n<Collection>).",
        "",
        f"Jet algorithm : {cfg.jet_algorithm}, R = {cfg.jet_r}, "
        f"pT > {cfg.jet_pt_min} GeV, |eta| < {cfg.jet_eta_max}",
        f"Detector      : {Path(cfg.delphes_card).name}",
        f"Generator     : Pythia 8, card {Path(cfg.pythia_card).name}",
        f"Random seed   : {cfg.seed}",
        "",
        "Read it with, for example:",
        "    import uproot",
        f"    t = uproot.open('{cfg.nanoaod_path.name}')['Events']",
        "    pt = t['Jet_pt'].array()",
        "",
        "Full provenance is in the 'Metadata' object of this file and in",
        f"{cfg.provenance_path.name}.",
        "",
        "BRANCHES",
        "-" * 70,
    ]
    for name in columns:
        doc = BRANCH_DOC.get(name, "(undocumented)")
        lines.append(f"{name}")
        for chunk in _wrap(doc, 66):
            lines.append(f"    {chunk}")
    return "\n".join(lines)


def _wrap(text: str, width: int):
    words, line = text.split(), ""
    for w in words:
        if len(line) + len(w) + 1 > width:
            yield line
            line = w
        else:
            line = f"{line} {w}".strip()
    if line:
        yield line


def write(cfg: PipelineConfig, log=None) -> dict:
    """Build and write ``data/<run_name>_nano.root``."""
    import awkward as ak
    import uproot

    log = log or get_logger("write_nanoaod", cfg.log_dir)
    banner(log, "STAGE 4/5  Writing NanoAOD-like ROOT output")

    for path, stage in ((cfg.jets_path, "3 (analysis/cluster_jets.py)"),
                        (cfg.delphes_path, "2 (simulation/run_delphes.py)")):
        if not path.is_file():
            raise FileNotFoundError(f"{path} not found -- run stage {stage} first")

    t0 = time.time()
    columns = build_events(cfg, log)
    n_events = len(columns["event"])

    provenance = build_provenance(cfg, "write_nanoaod", n_events=n_events)

    with uproot.recreate(cfg.nanoaod_path) as fout:
        fout["Events"] = columns
        # Metadata travels inside the file so it can never be separated from
        # the data: generator, seed, cards, software versions, jet definition.
        fout["ReadMe"] = readme_text(cfg, columns)
        fout["Metadata"] = json.dumps(provenance, indent=2)

    elapsed = time.time() - t0
    size_mb = cfg.nanoaod_path.stat().st_size / 1e6

    log.info("events       : %d", n_events)
    log.info("branches     : %d", len(columns))
    log.info("output       : %s (%.2f MB, %.2f kB/event)",
             cfg.nanoaod_path, size_mb, 1000 * size_mb / max(n_events, 1))
    undocumented = [c for c in columns if c not in BRANCH_DOC]
    if undocumented:
        log.warning("branches missing from BRANCH_DOC: %s", undocumented)
    log.info("wrote ReadMe and Metadata objects into the file")
    log.info("elapsed      : %.1f s", elapsed)

    result = {
        "n_events": int(n_events),
        "n_branches": len(columns),
        "branches": sorted(columns),
        "nanoaod_file": str(cfg.nanoaod_path),
        "nanoaod_size_mb": round(size_mb, 3),
        "elapsed_s": round(elapsed, 2),
    }
    merge_provenance(cfg.provenance_path,
                     build_provenance(cfg, "write_nanoaod", **result))
    return result


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Stage 4: merge jets + Delphes objects into a flat NanoAOD-like tree.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--config", help="JSON config file")
    ap.add_argument("--run-name", help="label for input/output files")
    ap.add_argument("--no-truth", action="store_true",
                    help="omit the Particle_* truth branches (smaller file)")
    args = ap.parse_args(argv)

    cfg = PipelineConfig.from_file(args.config) if args.config else PipelineConfig()
    if args.run_name: cfg.run_name = args.run_name
    if args.no_truth: cfg.store_truth_particles = False
    cfg.validate()
    cfg.ensure_dirs()

    write(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
