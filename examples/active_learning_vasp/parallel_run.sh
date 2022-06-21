#!/bin/sh
##SBATCH -p F144cpu #パーティションを指定
##SBATCH -N 144
##SBATCH -n 9216
##SBATCH -c 2
##SBATCH --time=00:30:00
module purge
module load intel_compiler/2019.5.281
module load openmpi/4.0.4-intel-2019.5.281

parallel --delay 0.2 -j 16 --joblog runtask.log  \
	 -a rundirs.txt ./run_vasp.sh 
sleep 30
