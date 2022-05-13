#!/bin/sh
rm -f active.out
echo start AL sample
mpiexec -np 4 --oversubscribe abics_activelearn input.toml
echo start parallel_run 1
sh parallel_run.sh
sleep 5

echo start AL final
mpiexec -np 4 --oversubscribe abics_activelearn input.toml

#train
echo start training
abics_train input.toml > train.out
echo Done
