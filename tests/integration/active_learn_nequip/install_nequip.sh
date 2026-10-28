#!/bin/sh

# This script installs pytorch, nequip, and allegro
# into python3 environment
#
# NOTE for developers:
# nequip is pinned to 0.6.2 (and nequip-allegro to 0.3.0) because the abICS
# nequip solver/trainer target the pre-0.7 API. That series uses numpy.in1d,
# which was removed in NumPy 2, so this test only works with numpy<2 and is
# therefore restricted to Python 3.9 in the CI matrix (Python 3.13 has no
# numpy 1.x wheels). Lifting this requires migrating abICS to nequip >=0.7.

set -ue

echo "python3 points to the following:"
which python3

echo

# Pin numpy<2 up front: nequip 0.6.2 declares "numpy" with no upper bound but
# uses numpy.in1d (removed in NumPy 2). Setting it before the nequip stack keeps
# the subsequent installs resolving against a compatible numpy.
python3 -m pip install "numpy<2"

python3 -m pip install torch==2.5.1
python3 -m pip install nequip==0.6.2
python3 -m pip install nequip-allegro==0.3.0

# Re-assert numpy<2 in case any of the above pulled in numpy>=2.
python3 -m pip install "numpy<2"
