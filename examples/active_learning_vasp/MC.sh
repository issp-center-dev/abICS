#!/bin/sh
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=00:30:00

# module purge
# module load intel_compiler/2019.5.281
# module load openmpi/4.0.4-intel-2019.5.281

# source XXX

#sleep 30
srun -n 8 abics_sampling input.toml >> abics_sampling.out

echo Done
