#!/usr/bin/env python3
"""Stage 1 -- event generation with Pythia 8.

What this stage does
--------------------
Pythia 8 simulates a proton-proton collision from first principles, in four
conceptual steps:

1. **Hard scatter.**  Two partons (quarks or gluons), one from each proton, are
   drawn from the parton distribution functions and scattered using the exact
   matrix element for the process selected in the command card -- here
   ``f fbar -> gamma*/Z -> e+ e-``.  This is the only part computed at fixed
   order in perturbation theory.

2. **Parton shower.**  The incoming and outgoing coloured partons radiate
   gluons, which radiate further.  Pythia generates this as a Markov chain
   ordered in decreasing transverse momentum, resumming the large logarithms
   that fixed-order calculations miss.  Initial-state radiation is what turns
   "Z -> ee" into "Z + jets".

3. **Multiparton interactions.**  The other partons in the two protons scatter
   too, producing the soft "underlying event" on top of the hard process.

4. **Hadronisation.**  QCD confinement means free quarks cannot survive.  Pythia
   uses the Lund string model: colour-connected partons are joined by strings
   whose potential energy grows with separation until the string fragments into
   hadrons.  Unstable hadrons are then decayed.  The output is a list of
   long-lived particles -- pions, kaons, photons, electrons, muons, neutrons,
   protons, neutrinos -- which is what a real detector would see.

Input
-----
A Pythia command file (``generators/cards/*.cmnd``) plus the run parameters in
:class:`majetstik.config.PipelineConfig` (number of events, random seed).

Output
------
``data/<run_name>.hepmc`` -- a HepMC3 ASCII file, the standard interchange
format, containing the full event record (all four steps above, not just the
final-state particles).  Stage 2 feeds this to Delphes.

Also appends a provenance entry (seed, versions, cross section) to
``data/<run_name>_provenance.json``.

Run standalone
--------------
    source setup_env.sh
    python3 generators/gen_pythia.py --events 200 --seed 42

How to modify
-------------
* **Different physics**: write a new ``.cmnd`` card and pass ``--pythia-card``.
  ``generators/cards/dijet.cmnd`` is a worked second example.
* **Turn physics off to see its effect**: set ``PartonLevel:MPI = off`` in the
  card and compare the jet multiplicity; or ``HadronLevel:Hadronize = off`` to
  get parton-level "jets".
* **More statistics**: raise ``--events``.  Generation is the slowest stage and
  scales linearly.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

# Allow `python3 generators/gen_pythia.py` without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from majetstik.config import PipelineConfig
from majetstik.hepmc import HepMC3Writer, count_events
from majetstik.logging_setup import banner, get_logger
from majetstik.provenance import build_provenance, merge_provenance


def generate(cfg: PipelineConfig, log=None) -> dict:
    """Run Pythia and write a HepMC3 file.

    Returns a dict with the number of events written and the cross section.
    """
    import pythia8

    log = log or get_logger("gen_pythia", cfg.log_dir)
    banner(log, f"STAGE 1/5  Pythia 8 event generation  ({cfg.n_events} events)")

    card = cfg.resolve(cfg.pythia_card)
    log.info("command card : %s", card)
    log.info("random seed  : %d", cfg.seed)
    log.info("output       : %s", cfg.hepmc_path)

    pythia = pythia8.Pythia()
    pythia.readFile(str(card))

    # The config always wins over the card, so that a run is defined by the
    # config alone.  Pythia applies settings in the order they are read.
    pythia.readString("Random:setSeed = on")
    pythia.readString(f"Random:seed = {cfg.seed}")
    pythia.readString(f"Beams:eCM = {cfg.beam_energy_cm}")

    if not pythia.init():
        raise RuntimeError(f"Pythia failed to initialise with card {card}")

    t0 = time.time()
    n_failed = 0
    report_every = max(1, cfg.n_events // 10)

    with HepMC3Writer(cfg.hepmc_path) as writer:
        for i in range(cfg.n_events):
            if not pythia.next():
                # A failed event is normal at a low rate (phase-space edge
                # cases); we skip it rather than abort the run.
                n_failed += 1
                continue
            writer.write(pythia.event, i)
            if (i + 1) % report_every == 0:
                rate = (i + 1) / (time.time() - t0)
                log.info("  generated %6d / %d events  (%.1f evt/s)",
                         i + 1, cfg.n_events, rate)
        n_written = writer.n_written

    elapsed = time.time() - t0
    pythia.stat()                                   # cross-section summary table

    # sigmaGen() is in mb; convert to pb (1 mb = 1e9 pb), the usual LHC unit.
    sigma_pb = pythia.infoPython().sigmaGen() * 1e9
    sigma_err_pb = pythia.infoPython().sigmaErr() * 1e9

    size_mb = cfg.hepmc_path.stat().st_size / 1e6
    log.info("wrote %d events (%d failed) in %.1f s  [%.1f evt/s]",
             n_written, n_failed, elapsed, n_written / max(elapsed, 1e-9))
    log.info("HepMC3 file  : %.1f MB", size_mb)
    log.info("cross section: %.4g +- %.2g pb", sigma_pb, sigma_err_pb)

    # Cheap insurance: a HepMC file that cannot be parsed makes Delphes produce
    # an empty ROOT file with no error message at all.  Catch it here.
    n_readable = count_events(cfg.hepmc_path)
    if n_readable != n_written:
        raise RuntimeError(
            f"HepMC3 self-check failed: wrote {n_written} events but only "
            f"{n_readable} can be read back from {cfg.hepmc_path}. "
            "Delphes would silently produce an empty file."
        )
    log.info("HepMC3 self-check passed (%d events readable)", n_readable)

    result = {
        "n_events_requested": cfg.n_events,
        "n_events_written": n_written,
        "n_events_failed": n_failed,
        "cross_section_pb": sigma_pb,
        "cross_section_err_pb": sigma_err_pb,
        "hepmc_file": str(cfg.hepmc_path),
        "hepmc_size_mb": round(size_mb, 2),
        "elapsed_s": round(elapsed, 2),
    }
    merge_provenance(cfg.provenance_path,
                     build_provenance(cfg, "gen_pythia", **result))
    return result


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Stage 1: generate events with Pythia 8 -> HepMC3.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--config", help="JSON config file (see majetstik/config.py)")
    ap.add_argument("--events", type=int, help="number of events to generate")
    ap.add_argument("--seed", type=int, help="random seed")
    ap.add_argument("--pythia-card", help="Pythia 8 command file")
    ap.add_argument("--run-name", help="label for output files")
    args = ap.parse_args(argv)

    cfg = PipelineConfig.from_file(args.config) if args.config else PipelineConfig()
    if args.events:      cfg.n_events = args.events
    if args.seed:        cfg.seed = args.seed
    if args.pythia_card: cfg.pythia_card = args.pythia_card
    if args.run_name:    cfg.run_name = args.run_name
    cfg.validate()
    cfg.ensure_dirs()

    generate(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
