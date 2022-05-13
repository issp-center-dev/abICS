#!/bin/sh

export OMP_NUM_THREADS=1

sh ./AL.sh
sh ./MC.sh
sh ./AL.sh
sh ./MC.sh
