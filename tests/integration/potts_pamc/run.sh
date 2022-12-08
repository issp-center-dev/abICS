#!/bin/sh

export OMP_NUM_THREADS=1

set -e

rm -f result.dat
mpiexec -np 4 --oversubscribe abics_sampling input.toml
python3 ./check.py
