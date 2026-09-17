"""MaJETstik -- an end-to-end particle-physics simulation pipeline.

The package bundles the pieces that every stage of the pipeline shares:

    config      -- the PipelineConfig dataclass and its JSON/YAML-ish loader
    logging_setup -- consistent console + file logging
    provenance  -- records seeds, software versions and card checksums
    hepmc       -- converts a Pythia 8 event record into a valid HepMC3 file
    fastjet     -- a just-in-time compiled FastJet kernel usable from Python

The physics stages themselves live in top-level directories so they read as
standalone, runnable examples:

    generators/gen_pythia.py    stage 1  hard scatter + shower + hadronisation
    simulation/run_delphes.py   stage 2  fast detector simulation
    analysis/cluster_jets.py    stage 3  anti-kt jet clustering + substructure
    analysis/write_nanoaod.py   stage 4  flat NanoAOD-like ROOT output
    analysis/analyze_jets.py    stage 5  plots and summary tables
"""

__version__ = "0.1.0"
__all__ = ["config", "logging_setup", "provenance", "hepmc", "fastjet"]
