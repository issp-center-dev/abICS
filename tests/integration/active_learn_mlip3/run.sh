#!/bin/sh

export OMP_NUM_THREADS=1

set -e

rm -f ./ALloop.progress
rm -f ./rundirs.txt
rm -f ./runtask.log

rm -rf ./AL0
rm -rf ./AL1
rm -rf ./MC0
rm -rf ./MC1
rm -rf ./mlip-3_XSF
rm -rf ./train0
rm -rf ./generate0
rm -rf ./baseinput

sh ./AL.sh
sh ./MC.sh
sh ./AL.sh
sh ./MC.sh

# Verify the outputs of AL0 -> train -> MC0 -> AL1 -> train -> MC1
python3 ../check_active_learn.py input.toml 2
echo OK
