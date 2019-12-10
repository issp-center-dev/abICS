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
mpijob -spawn -np 24 abics issp_qe.toml >> stdout.log
