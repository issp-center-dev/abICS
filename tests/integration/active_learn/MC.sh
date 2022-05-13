#!/bin/sh
mpiexec -np 4 --oversubscribe abicsAL input.toml
echo Done
