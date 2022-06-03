#!/bin/sh
#SBATCH -p i8cpu
#SBATCH -N 4
#SBATCH -n 512
#SBATCH -J spinel
#SBATCH -c 1
#SBATCH --time=0:30:00

# Run reference DFT calc.
# module purge
# module load intel_compiler/2019.5.281
# module load openmpi/4.0.4-intel-2019.5.281

# source XXX

echo start AL sample
srun -n 8 abics_mlref input.toml >> active.out
echo start parallel_run 1
sh parallel_run.sh

echo start AL final
srun -n 8 abics_mlref input.toml >> active.out

#train
# module purge
# module load intel_compiler/2020.2.254  
# module load intel_mpi/2020.2.254
echo start training
abics_train input.toml > train.out
echo Done
