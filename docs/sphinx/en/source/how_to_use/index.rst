.. _sec_basic_usage:

***************************
Basic Usage
***************************

.. highlight:: none

In abICS, all solver and sampler parameters other than the position coordinates must be prepared in advance.
To this end, abICS requires two input files: (1) a reference file for specifiying solver parameters that is formatted according to the solver input format, and (2) an input file for specifying parameters that control how abICS performs the configurational sampling.

.. _subsec_basic_reference:

Preparing a reference file
===========================

The user must prepare a reference file according to the input format of the solver to be used.
The path of the reference file is specified by ``base_input_dir`` in the ``[solver]`` section in the abICS input file (see below).
The coordinate information should not be written here because it will obviously change in the course of the simulation. 
The lattice sites are specified in a separate abICS input file (see below), 
and abICS will take care of generating the coordinates section at each sampling step.
The following is an example of a QE reference file.

.. literalinclude::  ../../../../../examples/standard/spinel/baseinput/scf.in

.. _subsec_basic_input:

Preparing an input file of abICS		    
================================

The input file of abICS is constructed by the following four sections:

1. [replica] section specifies the parameters of the replica exchange Monte Carlo part, such as the number of replicas, the temperature range, and the number of Monte Carlo steps.
  
2. [solver] section specifies the parameters for the (first principle calculation) solver, including the type of solver (VASP, QE,...), the path to the solver, and the directory containing reference input files (see :ref:`subsec_basic_reference` ).
   
3. [observer] section specifies the type of physical quantity to be calculated.

4. [config] section specifies the configuration of the alloy, etc.

For details, see :doc:`../inputfiles/index` .
The following is an example of an input file selecting QE as a solver.

.. literalinclude::  ../../../../../examples/standard/spinel/input_qe.toml

		     
Execution
========================

The number of processes specified here must be greater than or equal to the number of replicas.


::

 $ mpiexec -np 2 abics input.toml

This creates a directory named with the replica number under the current directory, and each replica runs the solver in it.
Here, `input.toml` is an input file for abICS (see :ref:`subsec_basic_input`).

Tips for the number of MPI processes
========================================
abICS uses the MPI library function ``MPI_Comm_spawn`` to run the solver.
This function executes another program on new MPI processes.

For example, consider that you have a parallel machine with 4 CPU cores per node
and want to run two replicas and to invoke solvers on 4 CPU cores.
If invoked as ``mpiexec -np2 abics input.toml``, the replica control processes A and B are started on the first two cores of node 0,
each starting a new four-parallel solver a and b.
Then, solver a fills the remaining two cores in node 0 and the first two cores in node 1,
and solver b is placed on the remaining two cores in node 1 and the first two cores in node 2.
This causes inter-node communication within the solver and reduces performance.

By taking more of the initial process, you can align the processes and prevent unnecessary internode communication.
In this example, ``mpiexec -np 4 abics input.toml`` will allow the replica control processes A and B to fill all the cores of the node 0, while solvers a and b fill nodes 1 and 2, respectively.


Comments on MPI implementation
====================================
In the ``MPI_Comm_spawn`` function, 
some MPI implementations require the information of "the number of processes that can be started in total" to be set in variable ``MPI_UNIVERSE_SIZE``.
In this section, we will comment on some MPI implementations including how to set ``MPI_UNIVERSE_SIZE``.

OpenMPI
~~~~~~~~~~~~~
``MPI_UNIVERSE_SIZE`` is automatically set to the number of the CPU cores available.
If you want more processes, you should pass the ``--oversubscribe`` option to the ``mpiexec`` command.

When one of the spawned processes returns a nonzero return code, all the OpenMPI processes will abort.
The ``--mca orte_abort_on_non_zero_status 0`` option allows you to ignore the return code.
Quantum ESPRESSO, for example, may return a nonzero value due to a floating point exception even if the calculation is successfully completed.

MPICH / Intel MPI
~~~~~~~~~~~~~~~~~~~~~
The ``-usize <num>`` option sets ``MPI_UNIVERSE_SIZE``.
However, MPICH and Intel MPI seem not to use this value in ``MPI_Comm_spawn``.

HPE (SGI) MPT
~~~~~~~~~~~~~~~~~~~
The ``-up <num>`` option sets ``MPI_UNIVERSE_SIZE``.
This must be set before the ``-np <num>`` option.

Others
~~~~~~~~~~
On large supercomputers, the vendor may provide a dedicated MPI execution script along with the job scheduler.
In this case, please refer to the manuals provided by the supercomputer site.
On the ISSP supercomputer systems, Sekirei and Enaga, for example, ``mpijob -spawn`` sets ``MPI_UNIVERSE_SIZE`` properly.


.. _solver_specific_notes:

Solver specific Notes
===========================
