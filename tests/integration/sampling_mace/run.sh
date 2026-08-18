#!/bin/sh

export OMP_NUM_THREADS=1

set -e

rm -f kTs.npy
rm -rf 0 1 baseinput

# Pre-cache the pretrained MACE-MP model in a single process before launching
# the MPI run. mace_mp downloads the model on first use; if that download is
# left to the parallel ranks, a failure in one rank crashes only that process
# and the surviving rank deadlocks on the replica-exchange MPI Recv (the run
# then hangs until the CI job's hard timeout). Downloading once up front makes
# the subsequent ranks load the cached model with no network access.
python3 -c "from mace.calculators import mace_mp; mace_mp(device='cpu', default_dtype='float64')"

mpiexec -np 2 --oversubscribe abics_sampling input.toml

python3 ./check.py
