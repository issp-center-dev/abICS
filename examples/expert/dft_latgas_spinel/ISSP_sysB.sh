#!/bin/sh

#QSUB -queue F144cpu
#QSUB -node 144
##QSUB -mpi 24
##QSUB -omp 1
##QSUB -place distribute
#QSUB -over false

#PBS -l walltime=24:00:00
#PBS -N spinel
#PBS -j oe
#PBS -m abe

cd ${PBS_O_WORKDIR}

#ulimit -c 0
#ulimit -s unlimited

#nmpi=1008
#nomp=1

#export OMP_NUM_THREADS=$nomp
#export KMP_AFFINITY=disabled
#export KMP_STACKSIZE=1g

#export MPI_BUFFER_MAX=65536
#export MPI_BUFS_PER_PROC=256
#export MPI_IB_RAILS=2

#date > stdout 2>&1

#(time mpiexec_mpt -np $nmpi omplace -c 0-23 -nt ${OMP_NUM_THREADS} ../src/vasp.5.3/vasp < /dev/null ) >> stdout 2>&1

#date >> stdout 2>&1

#mpijob -spawn -up 624 -np 48 python3.6 ./dft_spinel_mix_mpi-init.py 12 48 < /dev/null >> stdout.log

export MPI_IB_CONGESTED=enabled
mpijob -spawn -np 48 python3.6 ./spinel_catmix.py 48 48 < /dev/null >> stdout.log
mpijob -np 48 python3.6 analyze_result.py 48





