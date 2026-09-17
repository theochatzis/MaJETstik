#!/usr/bin/env python3
"""MaJETstik -- run the whole simulation pipeline with one command.

    source setup_env.sh
    python3 main.py

That generates 1000 pp -> Z/gamma* -> e+e- events, pushes them through a
CMS-like detector simulation, clusters anti-kT R = 0.4 jets with substructure,
writes a flat NanoAOD-like ROOT file and produces summary plots -- about three
minutes end to end.

The five stages
---------------
    1  generators/gen_pythia.py     hard scatter, shower, hadronisation  -> .hepmc
    2  simulation/run_delphes.py    detector response                    -> _delphes.root
    3  analysis/cluster_jets.py     anti-kT clustering + substructure     -> _jets.root
    4  analysis/write_nanoaod.py    flatten and merge                     -> _nano.root
    5  analysis/analyze_jets.py     summary table and figures             -> plots/

Each stage is an independent script that reads the previous stage's file, so
you can run them one at a time while learning, and re-run only the stage you
changed.  ``--from`` and ``--to`` do exactly that:

    python3 main.py --from 3            # re-cluster and redo everything after
    python3 main.py --to 2              # stop after the detector simulation
    python3 main.py --only 5            # just remake the plots

Everything a run depends on -- seed, cards, software versions -- is written to
``data/<run_name>_provenance.json`` and embedded in the output ROOT file, so a
run can always be reproduced or audited.

Examples
--------
    # a quick 200-event test
    python3 main.py --events 200 --run-name quicktest

    # QCD dijets instead of Drell-Yan, with large-radius jets
    python3 main.py --pythia-card generators/cards/dijet.cmnd \
                    --run-name dijet --events 2000 --jet-r 0.8

    # an ATLAS-like detector
    python3 main.py --delphes-card $DELPHES_DIR/cards/delphes_card_ATLAS.tcl
"""

# NOTE: this file deliberately avoids `from __future__ import annotations` and
# any syntax newer than Python 3.6, so that it still *parses* under an old
# system interpreter and can print the message below instead of dying with a
# SyntaxError.  The rest of the project targets Python 3.9+.
import sys

if sys.version_info < (3, 9):
    sys.stderr.write(
        "\n"
        "MaJETstik needs Python 3.9 or newer; this is Python {}.{}.{}.\n"
        "\n"
        "You are most likely using the system interpreter because the project\n"
        "environment has not been set up.  Run:\n"
        "\n"
        "    source setup_env.sh\n"
        "\n"
        "and try again.  See README.md -> Installation Guide.\n"
        "\n".format(*sys.version_info[:3]))
    raise SystemExit(2)

import argparse
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from majetstik.config import PipelineConfig
from majetstik.logging_setup import get_logger
from majetstik.provenance import software_versions

#: stage number -> (label, import path, function name)
STAGES = {
    1: ("Pythia 8 event generation", "generators.gen_pythia", "generate"),
    2: ("Delphes detector simulation", "simulation.run_delphes", "simulate"),
    3: ("FastJet jet clustering", "analysis.cluster_jets", "cluster"),
    4: ("NanoAOD-like output", "analysis.write_nanoaod", "write"),
    5: ("Analysis and plots", "analysis.analyze_jets", "analyze"),
}


def _import_stage(module_path: str, function: str):
    """Import a stage function.

    The stage directories are plain folders rather than packages, so we import
    by file location instead of relying on ``__init__.py`` files -- that keeps
    each stage script runnable on its own with ``python3 path/to/stage.py``.
    """
    import importlib.util

    folder, module_name = module_path.split(".")
    path = Path(__file__).resolve().parent / folder / f"{module_name}.py"
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return getattr(module, function)


def check_environment(log) -> None:
    """Fail early and helpfully if the LCG environment was not set up."""
    import shutil

    missing = [tool for tool in ("DelphesHepMC3", "fastjet-config")
               if shutil.which(tool) is None]
    try:
        import pythia8  # noqa: F401
    except ImportError:
        missing.append("pythia8 (Python module)")
    try:
        import uproot  # noqa: F401
    except ImportError:
        missing.append("uproot")

    if missing:
        log.error("missing from the environment: %s", ", ".join(missing))
        log.error("Run  `source setup_env.sh`  first "
                  "(see README.md -> Installation Guide).")
        raise SystemExit(2)


def run(cfg: PipelineConfig, first: int = 1, last: int = 5) -> dict:
    """Run stages ``first`` .. ``last`` in order."""
    log = get_logger("pipeline", cfg.log_dir)

    log.info("#" * 68)
    log.info("#  MaJETstik pipeline")
    log.info("#" * 68)
    log.info("run name     : %s", cfg.run_name)
    log.info("events       : %d", cfg.n_events)
    log.info("seed         : %d", cfg.seed)
    log.info("pythia card  : %s", cfg.pythia_card)
    log.info("delphes card : %s", cfg.delphes_card)
    log.info("jets         : %s R=%.1f, pT > %.1f GeV",
             cfg.jet_algorithm, cfg.jet_r, cfg.jet_pt_min)
    log.info("output dir   : %s", cfg.resolve(cfg.data_dir))
    versions = software_versions()
    log.info("software     : pythia %s | delphes %s | fastjet %s | root %s",
             versions["pythia8"], versions["delphes"],
             versions["fastjet"], versions["root"])
    log.info("stages       : %d .. %d", first, last)

    results, t_start = {}, time.time()
    for number in range(first, last + 1):
        label, module_path, function = STAGES[number]
        t0 = time.time()
        try:
            stage_fn = _import_stage(module_path, function)
            results[number] = stage_fn(cfg)
        except Exception as exc:
            log.error("stage %d (%s) FAILED: %s", number, label, exc)
            raise
        log.info("stage %d (%s) done in %.1f s", number, label, time.time() - t0)

    total = time.time() - t_start
    log.info("#" * 68)
    log.info("#  pipeline finished in %.1f s", total)
    log.info("#" * 68)
    if last >= 4:
        log.info("NanoAOD file : %s", cfg.nanoaod_path)
    if last >= 5 and cfg.make_plots:
        log.info("plots        : %s", cfg.plot_dir)
    log.info("provenance   : %s", cfg.provenance_path)
    log.info("")
    log.info("Next: inspect the output with")
    log.info("    python3 examples/read_output.py --run-name %s", cfg.run_name)
    return results


def build_config(args) -> PipelineConfig:
    cfg = PipelineConfig.from_file(args.config) if args.config else PipelineConfig()
    overrides = {
        "n_events": args.events, "seed": args.seed, "run_name": args.run_name,
        "pythia_card": args.pythia_card, "delphes_card": args.delphes_card,
        "jet_r": args.jet_r, "jet_pt_min": args.jet_pt_min,
        "jet_algorithm": args.algorithm,
    }
    for key, value in overrides.items():
        if value is not None:
            setattr(cfg, key, value)
    if args.no_plots:
        cfg.make_plots = False
    if args.no_truth:
        cfg.store_truth_particles = False
    if args.drop_hepmc:
        cfg.keep_hepmc = False
    return cfg


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        epilog="See README.md for the full guide.",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    g = ap.add_argument_group("run")
    g.add_argument("--config", help="JSON config file with all settings")
    g.add_argument("--run-name", help="label for all output files")
    g.add_argument("--events", type=int, help="number of events to generate")
    g.add_argument("--seed", type=int, help="master random seed")

    g = ap.add_argument_group("physics")
    g.add_argument("--pythia-card", help="Pythia 8 command file (.cmnd)")
    g.add_argument("--delphes-card", help="Delphes detector card (.tcl)")
    g.add_argument("--algorithm", choices=["antikt", "kt", "cambridge"],
                   help="jet algorithm")
    g.add_argument("--jet-r", type=float, help="jet radius R")
    g.add_argument("--jet-pt-min", type=float, help="minimum jet pT in GeV")

    g = ap.add_argument_group("stage selection")
    g.add_argument("--from", dest="first", type=int, default=1, choices=range(1, 6),
                   help="first stage to run (default 1)")
    g.add_argument("--to", dest="last", type=int, default=5, choices=range(1, 6),
                   help="last stage to run (default 5)")
    g.add_argument("--only", type=int, choices=range(1, 6),
                   help="run exactly one stage")

    g = ap.add_argument_group("output")
    g.add_argument("--no-plots", action="store_true", help="skip the figures")
    g.add_argument("--no-truth", action="store_true",
                   help="omit truth particles from the NanoAOD file")
    g.add_argument("--drop-hepmc", action="store_true",
                   help="delete the large intermediate HepMC file after stage 2")
    g.add_argument("--save-config", metavar="PATH",
                   help="write the resolved configuration to PATH and exit")

    args = ap.parse_args(argv)
    first, last = (args.only, args.only) if args.only else (args.first, args.last)
    if first > last:
        ap.error(f"--from {first} is after --to {last}")

    cfg = build_config(args)
    cfg.validate()
    cfg.ensure_dirs()

    if args.save_config:
        cfg.save(args.save_config)
        print(f"wrote {args.save_config}")
        return 0

    log = get_logger("pipeline", cfg.log_dir)
    check_environment(log)
    run(cfg, first, last)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
