#!/bin/sh
#SBATCH -p i8cpu #パーティションを指定
#SBATCH -N 8
#SBATCH -n 512
#SBATCH -J spinel
#SBATCH -c 2
#SBATCH --time=0:30:00

# Run reference DFT calc.
#module purge
#module load intel_compiler/2019.5.281
#module load openmpi/4.0.4-intel-2019.5.281

for i in {1..2} # 2 runs set in baseinput
do
    srun -n 8  abics_mlref input.toml >> abics_mlref.out
    sh parallel_run.sh
done

srun -n 8 abics_mlref input.toml >> abics_mlref.out

#train
module purge
module load intel_compiler/2020.2.254  
module load intel_mpi/2020.2.254
abics_train input.toml >> abics_train.out
