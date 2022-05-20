#!/bin/sh
#SBATCH -p i8cpu #パーティションを指定
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --time=00:30:00

module purge
module load intel_compiler/2019.5.281
module load openmpi/4.0.4-intel-2019.5.281

#sleep 30
srun -n 15 abics_sampling input_aenet.toml >> aenet.out
