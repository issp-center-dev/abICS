#!/bin/sh

# This script installs SevenNet into the python3 environment.
# The pretrained 7net-0 potential ships with the sevenn package.

set -ue

echo "python3 points to the following:"
which python3

echo

python3 -m pip install sevenn
