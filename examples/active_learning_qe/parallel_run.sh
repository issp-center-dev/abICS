#!/bin/sh
#SBATCH -p i8cpu
#SBATCH -N 8
#SBATCH -n 1024
#SBATCH -J spinel
#SBATCH -c 1
#SBATCH --time=0:30:00

RESTART=OFF # ON or OFF
if [ "_$RESTART" = "_ON" ]; then
	RESUME_OPT=--resume-failed
else
	RESUME_OPT=""
fi

parallel --delay 0.2 -j 32 --joblog runtask.log $RESUME_OPT  \
	 -a rundirs.txt "/bin/sh ./run_pw.sh"
sleep 30
