#!/bin/sh
parallel --delay 0.2 -j 8 --joblog runtask.log $RESUME_OPT  \
	 -a rundirs.txt abics_mocksolver
sleep 5
