#!/bin/sh
rm -f active.out
echo start AL sample
mpiexec -np 2 --oversubscribe abics_mlref input.toml
echo start parallel_run 1
sh parallel_run.sh
sleep 5

echo start AL final
mpiexec -np 2 --oversubscribe abics_mlref input.toml
sleep 5

#train
echo start training
abics_train input.toml

# for debug
echo '== cat train0/stdout =='
cat train0/stdout
echo '== end of train0/stdout =='

echo Done
