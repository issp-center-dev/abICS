#!/bin/sh
# Integrated mode: abics_mlref runs the mock solver itself (no mpiexec, no parallel_run.sh)
rm -f active.out
echo start AL integrated
abics_mlref input_integrated.toml

#train
echo start training
abics_train input_integrated.toml
echo 'cat train0/stdout'
cat train0/stdout
echo Done
