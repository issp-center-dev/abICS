#!/bin/sh

# This script installs MACE into the python3 environment.
# The pretrained MACE-MP foundation model is downloaded on first use.

set -ue

echo "python3 points to the following:"
which python3

echo

python3 -m pip install mace-torch
