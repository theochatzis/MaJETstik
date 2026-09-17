#!/usr/bin/env python3
"""Installation script for MaJETstik.

MaJETstik is normally used *without* installing: the pipeline scripts are
designed to be run in place, and ``source setup_env.sh`` already puts the
project on ``PYTHONPATH``.  Installing is useful when you want to
``import majetstik`` from your own analysis code elsewhere on the machine.

    # editable install, so edits to the source take effect immediately
    pip install --user -e .

    # optional extras
    pip install --user -e ".[notebook]"

Note that the heavy C++ dependencies -- Pythia 8, Delphes, FastJet and HepMC3 --
cannot be installed by pip.  They come from the CVMFS LCG view that
``setup_env.sh`` sources, or from your own build; see README.md.
"""

from pathlib import Path

from setuptools import find_packages, setup

ROOT = Path(__file__).parent
LONG_DESCRIPTION = (ROOT / "README.md").read_text(encoding="utf-8")

# Only the pip-installable requirements; the C++ stack is documented in
# requirements.txt and README.md.
INSTALL_REQUIRES = [
    "numpy>=1.24,<2.0",
    "uproot>=5.0",
    "awkward>=2.0",
    "pandas>=2.0",
    "matplotlib>=3.6",
]

setup(
    name="majetstik",
    version="0.1.0",
    description=("End-to-end particle-physics simulation pipeline: "
                 "Pythia 8 -> Delphes -> FastJet -> NanoAOD-like ROOT output"),
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    author="Theo Chatzistavrou",
    url="https://github.com/your-org/MaJETstik",
    license="MIT",
    packages=find_packages(include=["majetstik", "majetstik.*"]),
    python_requires=">=3.9",
    install_requires=INSTALL_REQUIRES,
    extras_require={
        "notebook": ["jupyterlab>=4.0", "particle>=0.23", "vector>=1.0"],
    },
    entry_points={
        "console_scripts": [
            "majetstik = main:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords="physics hep pythia delphes fastjet jets nanoaod lhc",
)
