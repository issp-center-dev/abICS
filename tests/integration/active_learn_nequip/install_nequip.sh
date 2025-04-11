#!/bin/sh

# This script installs pytorch, nequip, and allegro
# into python3 environment

set -ue

echo "python3 points to the following:"
which python3

echo

python3 -m pip install torch==2.5.1
python3 -m pip install nequip
python3 -m pip install git+https://github.com/mir-group/allegro.git
