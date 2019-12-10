.. highlight:: none

Basic usage
-----------

In abICS, physical quantities by solver (VASP, QE, etc.) are calculated while updating the position coordinates for each Monte Carlo step.
Therefore, it is necessary to prepare information other than the position coordinates in advance. This information is obtained by referring to the file according to the solver input format.

Prepare a reference file
===========================

Prepare an input file according to the input format of the solver to be used.
The path of the reference file is specified by ``base_input_dir`` in the ``[solver]`` section in the abICS input file.
The coordinate information does not need to be described because it refers to the abICS input file.
The following is an example of a QE reference file.

.. literalinclude::  ../../../../../examples/standard/spinel/baseinput/scf.in


Make an input file of abICS		    
===========================

The input file of abICS is constructed by the following four sections:

1. [replica] section
   Specify the parameters of the replica exchange Monte Carlo part, such as the number of replicas, the temperature range, and the number of Monte Carlo steps.
  
2. [solver] section
   Specify the parameters for the (first principle calculation) solver, including the type of solver (VASP, QE,...), the path to the solver, and the directory containing immutable input files.
   
3. [observer] section
   Specify the type of physical quantity to be calculated.

4. [config] section
   Specify the configuration of the alloy, etc.

For details, see :doc:`../file_specification/index` .
The following is an example of an input file selecting QE as a solver.

.. literalinclude::  ../../../../../examples/standard/spinel/input_qe.toml

		     
Execution
========================

The number of processes specified here must be the same as the number of replicas.

(Various matters derived from MPI_Comm_spawn (the actual process to be allocated, mpiexec options, etc. need to be described somewhere)

::

 $ mpiexec -np 2 abics input.toml

