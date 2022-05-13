#!/bin/sh
mpiexec -np 4 abicsAL input_aenet.toml >> aenet.out
echo Done
