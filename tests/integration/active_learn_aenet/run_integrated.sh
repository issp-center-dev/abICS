#!/bin/sh
# Same loop as run.sh, with abics_mlref in integrated mode (requires libEnsemble)

export OMP_NUM_THREADS=1

set -e

rm -f ./ALloop.progress
rm -f ./rundirs.txt
rm -f ./runtask.log

rm -rf ./AL0
rm -rf ./AL1
rm -rf ./MC0
rm -rf ./MC1
rm -rf ./aenetXSF
rm -rf ./train0
rm -rf ./generate0
rm -rf ./baseinput

sh ./AL_integrated.sh
sh ./MC.sh
sh ./AL_integrated.sh
sh ./MC.sh

# Verify the outputs of AL0 -> train -> MC0 -> AL1 -> train -> MC1
python3 ../check_active_learn.py input_integrated.toml 2
echo OK
