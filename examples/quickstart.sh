#!/bin/bash
# ---------------------------------------------------------------------------
# MaJETstik quick start -- the five-minute path, as a script.
#
#     bash examples/quickstart.sh
#
# Generates 200 pp -> Z/gamma* -> e+e- events, runs them through the full
# pipeline, and shows you the result.  Uses a small event count so it finishes
# in well under a minute; drop the --events flag for the default 1000.
# ---------------------------------------------------------------------------
# Note: no `set -u` here.  The CVMFS LCG setup script that setup_env.sh sources
# references unset variables internally, which would abort the script.
set -eo pipefail

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
cd "$HERE"

echo "=== 1/3  setting up the environment ==============================="
# shellcheck disable=SC1091
source setup_env.sh

echo
echo "=== 2/3  running the pipeline ====================================="
python3 main.py --run-name quickstart --events 200

echo
echo "=== 3/3  looking at the output ===================================="
python3 examples/read_output.py --run-name quickstart

cat <<MSG

------------------------------------------------------------------------
Done.  What you have now:

  data/quickstart_nano.root       the analysis file (flat NanoAOD-like tree)
  data/plots/                     six figures
  data/quickstart_provenance.json seed, cards and software versions
  data/logs/                      one log per stage

Next:
  python3 examples/read_output.py --run-name quickstart --readme
      ... what every branch means

  python3 main.py --config examples/config_dijet.json
      ... QCD dijets with large-radius jets, where substructure does real work

  docs/theory.md        the physics behind each stage
  docs/extending.md     recipes for changing the process, detector or observables
------------------------------------------------------------------------
MSG
