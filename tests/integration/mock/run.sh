#!/bin/sh

export OMP_NUM_THREADS=1

set -e

rm -f kTs.npy
mpiexec -np 2 --oversubscribe abics_sampling input.toml
python3 ./check.py
