"""Central configuration for the MaJETstik pipeline.

Every stage of the pipeline reads the *same* :class:`PipelineConfig` object, so
there is exactly one place that decides how many events to generate, which
random seed to use, which Delphes card describes the detector and how jets are
clustered.  That is what makes a run reproducible: save the config, re-run,
get identical numbers.

The config can be loaded from a JSON file::

    from majetstik.config import PipelineConfig
    cfg = PipelineConfig.from_file("my_run.json")

or built in Python and overridden per field::

    cfg = PipelineConfig(n_events=5000, seed=42)

Only fields you want to change need to appear in the JSON file; everything else
falls back to the defaults below.

Example JSON (see examples/config_dijet.json for a real one)::

    {
      "n_events": 2000,
      "seed": 12345,
      "pythia_card": "generators/cards/dijet.cmnd",
      "jet_r": 0.4,
      "jet_pt_min": 20.0
    }
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field, asdict, fields
from pathlib import Path
from typing import Any, Dict

# Project root = the directory containing this package.
PROJECT_ROOT = Path(__file__).resolve().parent.parent


@dataclass
class PipelineConfig:
    """All knobs for one end-to-end pipeline run.

    Paths may be given relative to the project root; :meth:`resolve` turns them
    into absolute paths.
    """

    # ---- run identity -----------------------------------------------------
    run_name: str = "zee_jets"
    """Label used for output file names, e.g. data/zee_jets_nano.root."""

    n_events: int = 1000
    """Number of events to generate.  1000 takes roughly a minute."""

    seed: int = 12345
    """Master random seed.

    Pythia is seeded with exactly this value (``Random:seed``).  Because
    Delphes' smearing uses ROOT's gRandom, we also seed that from the same
    number, so the whole chain is reproducible from this one integer.
    Pythia only accepts seeds in [0, 900000000]; 0 means "use a time-based
    seed", which would destroy reproducibility, so we forbid it.
    """

    # ---- stage 1: event generation ---------------------------------------
    pythia_card: str = "generators/cards/zee_jets.cmnd"
    """Pythia 8 command file describing the physics process."""

    beam_energy_cm: float = 13600.0
    """Centre-of-mass energy in GeV.  13600 = LHC Run 3 (13.6 TeV)."""

    # ---- stage 2: detector simulation ------------------------------------
    delphes_card: str = "simulation/cards/majetstik_cms.tcl"
    """Delphes detector description (a Tcl 'card')."""

    # ---- stage 3: jet clustering -----------------------------------------
    jet_algorithm: str = "antikt"
    """Sequential recombination algorithm: 'antikt', 'kt' or 'cambridge'."""

    jet_r: float = 0.4
    """Jet radius parameter R in (eta, phi) space."""

    jet_pt_min: float = 20.0
    """Only keep clustered jets above this pT (GeV).  Below ~15 GeV jets are
    dominated by the underlying event and are not calibratable."""

    jet_eta_max: float = 4.7
    """Drop jets outside the detector acceptance."""

    # N-subjettiness: tau_N measures how well a jet is described by N subjets.
    nsub_beta: float = 1.0
    """Angular exponent of the N-subjettiness measure.  beta=1 is the standard
    choice for quark/gluon and boosted-object tagging."""

    # Soft drop grooming: recursively drop soft, wide-angle radiation.
    softdrop_beta: float = 0.0
    """beta=0 makes soft drop identical to the modified Mass Drop Tagger."""

    softdrop_zcut: float = 0.1
    """Momentum-fraction threshold for the soft-drop condition."""

    # ---- stage 4: NanoAOD output -----------------------------------------
    store_truth_particles: bool = True
    """Write the generator-level Particle_* branches."""

    truth_status_filter: str = "final_or_hard"
    """Which truth particles to keep: 'final' (status 1 only), 'final_or_hard'
    (status 1 plus the hard-process partons) or 'all'.  'all' makes the output
    file several times larger."""

    truth_pt_min: float = 0.5
    """Drop truth particles below this pT to keep the file small."""

    store_constituents: bool = True
    """Write the jet->constituent index map (see docs/output_format.md)."""

    # ---- stage 5: analysis ------------------------------------------------
    make_plots: bool = True

    # ---- bookkeeping ------------------------------------------------------
    data_dir: str = "data"
    keep_hepmc: bool = True
    """HepMC files are large (~10 MB per 100 events).  Set False to delete the
    intermediate file once Delphes has consumed it."""

    extra: Dict[str, Any] = field(default_factory=dict)
    """Free-form extras, recorded in the provenance block but not interpreted."""

    # ------------------------------------------------------------------
    # derived paths
    # ------------------------------------------------------------------
    @property
    def hepmc_path(self) -> Path:
        return self.resolve(self.data_dir) / f"{self.run_name}.hepmc"

    @property
    def delphes_path(self) -> Path:
        return self.resolve(self.data_dir) / f"{self.run_name}_delphes.root"

    @property
    def jets_path(self) -> Path:
        return self.resolve(self.data_dir) / f"{self.run_name}_jets.root"

    @property
    def nanoaod_path(self) -> Path:
        return self.resolve(self.data_dir) / f"{self.run_name}_nano.root"

    @property
    def provenance_path(self) -> Path:
        return self.resolve(self.data_dir) / f"{self.run_name}_provenance.json"

    @property
    def plot_dir(self) -> Path:
        return self.resolve(self.data_dir) / "plots"

    @property
    def log_dir(self) -> Path:
        return self.resolve(self.data_dir) / "logs"

    # ------------------------------------------------------------------
    @staticmethod
    def resolve(path: str | os.PathLike) -> Path:
        """Turn a project-relative path into an absolute one."""
        p = Path(path)
        return p if p.is_absolute() else (PROJECT_ROOT / p)

    def validate(self) -> None:
        """Fail fast on values that would silently produce nonsense."""
        if self.n_events <= 0:
            raise ValueError(f"n_events must be positive, got {self.n_events}")
        if not 0 < self.seed <= 900_000_000:
            raise ValueError(
                f"seed must be in (0, 900000000]; got {self.seed}. "
                "Pythia reserves 0 for a time-based (non-reproducible) seed."
            )
        if not 0 < self.jet_r <= 2.0:
            raise ValueError(f"jet_r must be in (0, 2], got {self.jet_r}")
        if self.jet_algorithm not in ("antikt", "kt", "cambridge"):
            raise ValueError(
                f"unknown jet_algorithm {self.jet_algorithm!r}; "
                "expected 'antikt', 'kt' or 'cambridge'"
            )
        if not 0.0 <= self.softdrop_zcut < 1.0:
            raise ValueError(f"softdrop_zcut must be in [0, 1), got {self.softdrop_zcut}")
        if self.truth_status_filter not in ("final", "final_or_hard", "all"):
            raise ValueError(
                f"unknown truth_status_filter {self.truth_status_filter!r}"
            )
        for label, p in (("pythia_card", self.pythia_card),
                         ("delphes_card", self.delphes_card)):
            if not self.resolve(p).is_file():
                raise FileNotFoundError(f"{label} not found: {self.resolve(p)}")

    def ensure_dirs(self) -> None:
        for d in (self.resolve(self.data_dir), self.plot_dir, self.log_dir):
            d.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def save(self, path: str | os.PathLike) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2) + "\n")

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "PipelineConfig":
        known = {f.name for f in fields(cls)}
        unknown = set(data) - known
        if unknown:
            raise ValueError(
                f"unknown config keys: {sorted(unknown)}. "
                f"Valid keys: {sorted(known)}"
            )
        return cls(**data)

    @classmethod
    def from_file(cls, path: str | os.PathLike) -> "PipelineConfig":
        return cls.from_dict(json.loads(Path(path).read_text()))
