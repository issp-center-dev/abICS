#!/bin/sh
mpiexec -np 4 --oversubscribe abicsAL input_aenet.toml >> aenet.out
echo Done
