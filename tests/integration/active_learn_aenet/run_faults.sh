#!/bin/sh
# Fault injection into the integrated abics_mlref and the abics_al driver
# (requires libEnsemble): failing, hanging and killed solver runs, the driver
# killed in the middle of a stage, the wall-time stop and inconsistent
# directories.  See run_faults.py for the scenarios.

export OMP_NUM_THREADS=1

set -e

python3 ./run_faults.py
echo OK
