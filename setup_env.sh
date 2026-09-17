#!/bin/bash
# ---------------------------------------------------------------------------
# MaJETstik environment setup
#
# Everything MaJETstik needs (Pythia 8, Delphes, FastJet, ROOT, uproot, numpy,
# ...) already exists as a prebuilt "LCG view" on CVMFS.  A view is just a
# directory tree of symlinks; sourcing its setup.sh puts all of it on PATH,
# LD_LIBRARY_PATH and PYTHONPATH at once.  There is nothing to compile.
#
#     source setup_env.sh
#
# If you are NOT on a CVMFS machine (lxplus, CERN/WLCG grid nodes, or any host
# with cvmfs mounted), see README.md -> "Installation without CVMFS".
# ---------------------------------------------------------------------------

# Pick the newest LCG view that matches this machine's OS.  LCG_107 is the last
# release built for EL8; newer ones (LCG_110) are EL9-only.
MAJETSTIK_LCG_VIEW="${MAJETSTIK_LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el8-gcc11-opt}"

if [ ! -d "$MAJETSTIK_LCG_VIEW" ]; then
    echo "ERROR: LCG view not found: $MAJETSTIK_LCG_VIEW" >&2
    echo "       Set MAJETSTIK_LCG_VIEW to a view that exists on this machine," >&2
    echo "       e.g. ls /cvmfs/sft.cern.ch/lcg/views/" >&2
    return 1 2>/dev/null || exit 1
fi

source "$MAJETSTIK_LCG_VIEW/setup.sh"

# Delphes ships its executables (DelphesHepMC3, ...) and its stock detector
# cards in the same release area.  We export the card directory so scripts can
# offer the stock cards (delphes_card_CMS.tcl, delphes_card_ATLAS.tcl, ...)
# alongside the MaJETstik one in simulation/cards/.
_delphes_bin="$(command -v DelphesHepMC3 2>/dev/null)"
if [ -n "$_delphes_bin" ]; then
    export DELPHES_DIR="$(dirname "$(dirname "$(readlink -f "$_delphes_bin")")")"
fi
unset _delphes_bin

# Make `import majetstik` work from anywhere.
MAJETSTIK_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export MAJETSTIK_ROOT
export PYTHONPATH="$MAJETSTIK_ROOT${PYTHONPATH:+:$PYTHONPATH}"

echo "MaJETstik environment ready"
echo "  LCG view   : $MAJETSTIK_LCG_VIEW"
echo "  project    : $MAJETSTIK_ROOT"
echo "  python     : $(python3 --version 2>&1)"
echo "  Delphes    : ${DELPHES_DIR:-<not found>}"
