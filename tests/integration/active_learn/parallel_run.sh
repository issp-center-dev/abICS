#!/bin/sh
parallel --delay 0.2 -j 8 --joblog runtask.log $RESUME_OPT  \
  -a rundirs.txt python3 ./mock_solver.py
sleep 5
