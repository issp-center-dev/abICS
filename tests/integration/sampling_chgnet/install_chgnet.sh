#!/bin/sh

# This script installs CHGNet into the python3 environment.
# CHGNet ships a pretrained model, so the smoke test needs no network at runtime.

set -ue

echo "python3 points to the following:"
which python3

echo

python3 -m pip install chgnet
