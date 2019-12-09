.. highlight:: none

Overview
-----------

Preparing the input file
========================


Execution
========================

The number of processes specified here must be the same as the number of replicas.

(Various matters derived from MPI_Comm_spawn (the actual process to be allocated, mpiexec options, etc. need to be described somewhere)

::

 $ mpiexec -np 2 abics input.toml

