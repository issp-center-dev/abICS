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

The number of processes specified here must be greater than or equal to the number of replicas.


::

 $ mpiexec -np 2 abics input.toml

This creates a directory named with the replica number under the current directory, and each replica runs the solver in it.


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

By taking more of the initial process, you can align the processes and prevent unnecessary node stripping.
In this example, ``mpiexec -np 4 abics input.toml`` will allow the replica control processes A and B to fill all the cores of the node 0, while solvers a and b fill nodes 1 and 2, respectively.


Comments on MPI implementation
====================================
In the ``MPI_Comm_spawn`` function, a MPI implementation can use the information "how many process can be started in total" by ``MPI_UNIVERSE_SIZE``.
In this section, we will comment on some MPI implementations including how to set ``MPI_UNIVERSE_SIZE``.

OpenMPI
~~~~~~~~~~~~~
``MPI_UNIVERSE_SIZE`` is automatically set to the number of the CPU cores available.
If you want more processes, you should pass the ``--oversubscribe`` option to the ``mpiexec`` command.

When one of the spawned process returns a nonzero return code, all the OpenMPI processes will abort.
The ``--mca orte_abort_on_non_zero_status 0`` option allows you to ignore the return code.
Quantum ESPRESSO, for example, may return a nonzero value due to a floating point exception even if the calculation is successfully completed.

MPICH / Intel MPI
~~~~~~~~~~~~~~~~~~~~~
The ``-usize <num>`` option set ``MPI_UNIVERSE_SIZE``.
Indeed, MPICH and Intel MPI seem not to use this value in ``MPI_Comm_spawn``.

HPE (SGI) MPT
~~~~~~~~~~~~~~~~~~~
The ``-up <num>`` option set ``MPI_UNIVERSE_SIZE``.
This must be set before the ``-np <num>`` option.

Others
~~~~~~~~~~
On large supercomputers, the vendor may provide a dedicated MPI execution script along with the job scheduler.
In this case, refer to the manuals.
On the ISSP supercomputer systems, Sekirei and Enaga, for example, ``mpijob -spawn`` sets ``MPI_UNIVERSE_SIZE`` properly.


The output of abics
---------------------
The calculation results are output below each replica directory.

``structure.XXX.vasp``
=========================
The atomic position for each step in the POSCAR file format of VASP is saved.
``XXX`` in the filename means the index of the step.

``minE.vasp``
====================
The atomic position minimizing the total energy in the POSCAR file format of VASP is saved.

``obs.dat``
===================
The temperature and the total energy for each step in units of eV is saved.

``obs_save.npy``
==================
The total energy for each step in units of eV in the Numpy binary format is saved.
Users can load it as ``darray`` by using ``numpy.load('obs_save.npy')``.

``kT_hist.npy``
==================
The temperature for each step in units of eV in the Numpy binary format is saved.
Users can load it as ``darray`` by using ``numpy.load('kT_hist.npy')``.

``Trank_hist.npy``
==================
The rank (index) of the temperature for each step in the Numpy binary format is saved.
Users can load it as ``darray`` by using ``numpy.load('Trank_hist.npy')``.

