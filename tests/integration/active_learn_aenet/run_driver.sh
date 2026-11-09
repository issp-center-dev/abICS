#!/bin/sh
# The whole loop with the abics_al driver (requires libEnsemble):
# two iterations in one invocation, then a third one by raising max_iterations.

export OMP_NUM_THREADS=1

set -e

rm -f ./ALloop.progress ./rundirs.txt ./runtask.log ./al_summary.dat
rm -rf ./AL0 ./AL1 ./AL2 ./MC0 ./MC1 ./MC2 ./aenetXSF ./train0 ./generate0 ./baseinput ./al_driver

echo start abics_al, max_iterations = 2
abics_al input_driver.toml
python3 ../check_active_learn.py input_driver.toml 2
test -f al_summary.dat
test -f al_driver/model_metadata.json

echo continue with max_iterations = 3
abics_al input_driver.toml --max-iterations 3
python3 ../check_active_learn.py input_driver.toml 3
cat al_summary.dat
echo OK
