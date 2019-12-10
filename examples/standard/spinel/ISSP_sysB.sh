#!/bin/sh

#QSUB -queue i18cpu
#QSUB -node 18
#QSUB -over false

#PBS -l walltime=0:30:00
#PBS -N spinel
#PBS -j oe
#PBS -m abe

cd ${PBS_O_WORKDIR}

export MPI_IB_CONGESTED=enabled
mpijob -spawn -np 24 python3 ./spinel_catmix.py issp.toml >> stdout.log
# mpijob -np 17 python3 analyze_result.py 17
