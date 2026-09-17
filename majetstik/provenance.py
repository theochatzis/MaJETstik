"""Record *how* a data file was made, so it can be remade identically.

A ROOT file full of jets is only scientifically useful if you can answer: which
generator version, which seed, which detector card, which jet radius?  This
module collects that information into a plain dictionary which is

  * written next to the data as ``<run>_provenance.json``, and
  * embedded inside the NanoAOD ROOT file itself (see
    :mod:`analysis.write_nanoaod`) so the two can never drift apart.

Nothing here affects the physics; it is pure bookkeeping.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional


def file_checksum(path: str | Path, algorithm: str = "sha256") -> Optional[str]:
    """Return a hex digest of a file, or None if it does not exist.

    Used for the Pythia and Delphes cards: two runs with the same checksums and
    the same seed must produce the same events.
    """
    p = Path(path)
    if not p.is_file():
        return None
    h = hashlib.new(algorithm)
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return f"{algorithm}:{h.hexdigest()}"


def _safe(fn, default="unavailable"):
    try:
        return fn()
    except Exception:
        return default


def _pythia_version() -> str:
    def go():
        import pythia8
        p = pythia8.Pythia("", False)          # False = don't print the banner
        return f"{p.settings.parm('Pythia:versionNumber'):.3f}"
    return _safe(go)


def _fastjet_version() -> str:
    def go():
        from majetstik.fastjet import fastjet_version
        return fastjet_version()
    return _safe(go)


def _delphes_version() -> str:
    """Delphes version, e.g. '3.5.1pre09'.

    The LCG ``version.txt`` records the full dependency set as
    ``delphes=3.5.1pre09:d92c8/Boost=.../ROOT=...``; we keep only the Delphes
    version itself so the log stays readable.  The full string is recoverable
    from the installation path, which is also recorded.
    """
    def go():
        import os
        d = os.environ.get("DELPHES_DIR")
        version_file = Path(d) / "version.txt" if d else None
        if version_file and version_file.is_file():
            raw = version_file.read_text().strip()
            head = raw.split("/")[0]                   # 'delphes=3.5.1pre09:d92c8'
            if "=" in head:
                return head.split("=", 1)[1].split(":")[0]
            return head
        exe = Path(subprocess.check_output(["which", "DelphesHepMC3"],
                                           text=True).strip()).resolve()
        return exe.parts[-3]                           # '<version>-<hash>' dir
    return _safe(go)


def _root_version() -> str:
    return _safe(lambda: subprocess.check_output(["root-config", "--version"],
                                                 text=True).strip())


def _git_commit() -> str:
    def go():
        root = Path(__file__).resolve().parent.parent
        return subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "--short", "HEAD"],
            text=True, stderr=subprocess.DEVNULL).strip()
    return _safe(go, default="not-a-git-checkout")


def software_versions() -> Dict[str, str]:
    """Version strings for every piece of software the pipeline touches."""
    import importlib

    versions = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "pythia8": _pythia_version(),
        "delphes": _delphes_version(),
        "fastjet": _fastjet_version(),
        "root": _root_version(),
        "majetstik_git": _git_commit(),
    }
    for mod in ("numpy", "uproot", "awkward", "pandas", "matplotlib"):
        versions[mod] = _safe(
            lambda m=mod: importlib.import_module(m).__version__)
    return versions


def build_provenance(cfg, stage: str, **extra: Any) -> Dict[str, Any]:
    """Assemble the provenance record for one stage.

    Parameters
    ----------
    cfg : PipelineConfig
    stage : str
        Which pipeline stage produced the file ('gen_pythia', 'run_delphes', ...).
    **extra
        Stage-specific facts worth recording, e.g. the Pythia cross section.
    """
    return {
        "stage": stage,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "config": cfg.to_dict(),
        "seed": cfg.seed,
        "cards": {
            "pythia_card": str(cfg.resolve(cfg.pythia_card)),
            "pythia_card_checksum": file_checksum(cfg.resolve(cfg.pythia_card)),
            "delphes_card": str(cfg.resolve(cfg.delphes_card)),
            "delphes_card_checksum": file_checksum(cfg.resolve(cfg.delphes_card)),
        },
        "software": software_versions(),
        **extra,
    }


def merge_provenance(path: str | Path, record: Dict[str, Any]) -> Dict[str, Any]:
    """Append ``record`` to the run's provenance JSON file and return the whole.

    The file accumulates one entry per stage under ``stages``, so a finished run
    carries the full history from generation to analysis.
    """
    p = Path(path)
    doc: Dict[str, Any] = {"run": record["config"]["run_name"], "stages": {}}
    if p.is_file():
        try:
            doc = json.loads(p.read_text())
            doc.setdefault("stages", {})
        except json.JSONDecodeError:
            pass                                   # corrupt file: start fresh
    doc["stages"][record["stage"]] = record
    doc["updated_utc"] = record["created_utc"]
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(doc, indent=2) + "\n")
    return doc


def summarize(record: Dict[str, Any]) -> str:
    """One-line-per-item human summary, for the log."""
    sw = record.get("software", {})
    lines = [
        f"stage        : {record.get('stage')}",
        f"created      : {record.get('created_utc')}",
        f"seed         : {record.get('seed')}",
        f"pythia       : {sw.get('pythia8')}",
        f"delphes      : {sw.get('delphes')}",
        f"fastjet      : {sw.get('fastjet')}",
        f"root         : {sw.get('root')}",
    ]
    return "\n".join(lines)
