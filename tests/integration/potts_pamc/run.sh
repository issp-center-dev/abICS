#!/bin/sh

export OMP_NUM_THREADS=1

set -e

rm -f result.dat
mpiexec -np 4 --oversubscribe abics_sampling input.toml

if [ -e result.dat ] ; then
  echo OK
  exit 0
else
  echo FAILED
  exit 1
fi
