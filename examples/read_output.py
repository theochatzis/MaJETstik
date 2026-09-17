#!/usr/bin/env python3
"""A minimal, complete example: open the output file and look at it.

This is the smallest useful thing you can do with a MaJETstik NanoAOD file.
Run it, read it, then start changing it.

    source setup_env.sh
    python3 examples/read_output.py                      # default run
    python3 examples/read_output.py --run-name dijet     # another run
    python3 examples/read_output.py --readme             # branch documentation

Everything here is plain uproot, numpy and awkward -- no MaJETstik imports are
required to read the file, which is the whole point of the NanoAOD format.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from majetstik.config import PipelineConfig


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run-name", default="zee_jets")
    ap.add_argument("--file", help="path to a NanoAOD file (overrides --run-name)")
    ap.add_argument("--readme", action="store_true",
                    help="print the branch documentation stored in the file")
    ap.add_argument("--metadata", action="store_true",
                    help="print the provenance record stored in the file")
    args = ap.parse_args(argv)

    import awkward as ak
    import uproot

    path = Path(args.file) if args.file else (
        PipelineConfig(run_name=args.run_name).nanoaod_path)
    if not path.is_file():
        print(f"{path} does not exist.\nRun the pipeline first:\n"
              f"    python3 main.py --run-name {args.run_name}")
        return 1

    with uproot.open(path) as f:
        if args.readme:
            print(str(f["ReadMe"]))
            return 0
        if args.metadata:
            print(str(f["Metadata"]))
            return 0

        events = f["Events"]
        print(f"file    : {path}")
        print(f"events  : {events.num_entries}")
        print()

        # ------------------------------------------------------------------
        # 1. The jets of the first few events.
        # ------------------------------------------------------------------
        data = events.arrays([
            "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_nConstituents",
            "Jet_tau1", "Jet_tau2", "Jet_softdrop_mass", "Jet_isLepton",
            "Jet_btag", "Jet_flavor", "MET_pt", "Electron_pt",
        ])

        print("First 3 events, jet by jet")
        print("-" * 78)
        header = (f"{'evt':>3} {'jet':>3} {'pT':>8} {'eta':>7} {'phi':>7} "
                  f"{'mass':>7} {'nCons':>6} {'tau21':>6} {'mSD':>7} "
                  f"{'btag':>5} {'flav':>5} {'lep?':>5}")
        print(header)
        for i in range(min(3, events.num_entries)):
            for j in range(len(data["Jet_pt"][i])):
                t1 = data["Jet_tau1"][i][j]
                t2 = data["Jet_tau2"][i][j]
                tau21 = t2 / t1 if t1 > 0 else float("nan")
                print(f"{i:>3} {j:>3} {data['Jet_pt'][i][j]:>8.1f} "
                      f"{data['Jet_eta'][i][j]:>7.2f} {data['Jet_phi'][i][j]:>7.2f} "
                      f"{data['Jet_mass'][i][j]:>7.2f} "
                      f"{data['Jet_nConstituents'][i][j]:>6d} "
                      f"{tau21:>6.2f} {data['Jet_softdrop_mass'][i][j]:>7.2f} "
                      f"{data['Jet_btag'][i][j]:>5d} {data['Jet_flavor'][i][j]:>5d} "
                      f"{data['Jet_isLepton'][i][j]:>5d}")
        print()

        # ------------------------------------------------------------------
        # 2. Simple whole-sample numbers.  Note the Jet_isLepton cut: particle
        #    flow hands electrons to the jet algorithm too, so a Z -> ee event
        #    really does contain two "jets" that are electrons.
        # ------------------------------------------------------------------
        hadronic = data["Jet_isLepton"] == 0
        jet_pt = ak.flatten(data["Jet_pt"][hadronic])
        print("Whole sample")
        print("-" * 78)
        print(f"  hadronic jets        : {len(jet_pt)}")
        if len(jet_pt):
            pts = np.asarray(jet_pt)
            print(f"  jet pT mean / median : {pts.mean():.1f} / "
                  f"{np.median(pts):.1f} GeV")
            print(f"  hardest jet          : {pts.max():.1f} GeV")
        print(f"  mean MET             : {float(ak.mean(data['MET_pt'])):.1f} GeV")
        print(f"  events with >= 2 electrons : "
              f"{int(ak.sum(ak.num(data['Electron_pt']) >= 2))}")
        print()

        # ------------------------------------------------------------------
        # 3. Follow the jet -> constituent map.  NanoAOD trees hold only one
        #    level of variable length, so the map is stored as two flat,
        #    parallel arrays rather than as a nested list.
        # ------------------------------------------------------------------
        if "JetConstituent_jetIdx" in events.keys():
            cons = events.arrays(["JetConstituent_jetIdx", "JetConstituent_pfIdx",
                                  "PFCand_pt", "PFCand_pdgId"], entry_stop=1)
            jet_idx = np.asarray(cons["JetConstituent_jetIdx"][0])
            pf_idx = np.asarray(cons["JetConstituent_pfIdx"][0])
            pf_pt = np.asarray(cons["PFCand_pt"][0])
            pf_id = np.asarray(cons["PFCand_pdgId"][0])
            if len(jet_idx):
                which = pf_idx[jet_idx == 0]                # constituents of jet 0
                print("Constituents of jet 0 in event 0")
                print("-" * 78)
                print(f"  {len(which)} candidates, "
                      f"scalar sum pT = {pf_pt[which].sum():.1f} GeV")
                order = which[np.argsort(-pf_pt[which])][:5]
                for k in order:
                    print(f"    pdgId {int(pf_id[k]):>6}   pT {pf_pt[k]:>7.2f} GeV")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
