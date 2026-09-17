"""Consistent logging for every pipeline stage.

Each stage writes to the console *and* to ``data/logs/<stage>.log`` so that a
run can be audited afterwards.  Reproducibility is one of the project's goals
(see README "Reproducibility"), and a log that records the seed, the card files
and the software versions is half of that story -- the other half is
:mod:`majetstik.provenance`.

Usage::

    from majetstik.logging_setup import get_logger
    log = get_logger("gen_pythia", log_dir=cfg.log_dir)
    log.info("generating %d events", cfg.n_events)
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

_FORMAT = "%(asctime)s [%(levelname)-7s] %(name)-14s | %(message)s"
_DATEFMT = "%H:%M:%S"


def get_logger(name: str,
               log_dir: Optional[Path] = None,
               level: int = logging.INFO) -> logging.Logger:
    """Return a logger that prints to stderr and optionally tees to a file.

    Calling this twice with the same ``name`` returns the same logger without
    adding duplicate handlers, so stages are safe to import from ``main.py``.
    """
    log = logging.getLogger(name)
    log.setLevel(level)
    log.propagate = False

    if log.handlers:                       # already configured
        return log

    formatter = logging.Formatter(_FORMAT, datefmt=_DATEFMT)

    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(formatter)
    log.addHandler(console)

    if log_dir is not None:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_dir / f"{name}.log", mode="w")
        fh.setFormatter(formatter)
        log.addHandler(fh)

    return log


def banner(log: logging.Logger, title: str) -> None:
    """Print a visually obvious stage separator."""
    log.info("=" * 68)
    log.info(title)
    log.info("=" * 68)
