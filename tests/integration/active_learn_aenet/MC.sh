#!/bin/sh
mpiexec -np 2 --oversubscribe abics_sampling input.toml
echo Done
