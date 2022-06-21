#!/bin/sh

export OMP_NUM_THREADS=1

set -e

mpiexec -np 2 --oversubscribe abics_sampling input.toml

if [ -e kTs.npy ] ; then
  echo OK
  exit 0
else
  echo FAILED
  exit 1
fi
