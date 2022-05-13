#!/bin/sh
rm -f active.out
echo start AL sample
mpiexec -np 4 --oversubscribe abics_activelearn input_aenet.toml >> active.out
echo start parallel_run 1
sh parallel_run.sh
sleep 5

echo start AL final
mpiexec -np 4 --oversubscribe abics_activelearn input_aenet.toml >> active.out

#train
echo start training
abics_train input_aenet.toml > train.out
echo Done
