.. _sec_tutorial:

***************************
Tutorial
***************************

This section contains a tutorial on how to use Quantum ESPRESSO (QE) to perform active learning and structure estimation.
A set of input files used in this tutorial can be found in ``examples/standard/active_learning_qe``.
In the following, we assume that gnu parallel and anet are installed.
We will also use ohtaka, the supercomputer system B of the Institute for Solid State Physics, as the environment for running the calculations.


Preparation of input files
-----------------------------

Preparation of the QE reference file
=========================================

Place the input file to be referenced in the QE scf calculation in ``baseinput_ref``.
The following is a description of the ``scf.in`` file in the sample directory.

.. code-block::

    &CONTROL
    calculation = 'relax'
    tstress = .false.
    tprnfor = .false.
    pseudo_dir = '~/qe/pot'
    disk_io = 'low'
    wf_collect = .false.
    /
    &SYSTEM
      ecutwfc = 60.0
      occupations = "smearing".
      smearing = "gauss"
      degauss = 0.01
    /
    &electrons
      mixing_beta = 0.7
      conv_thr = 1.0d-8
      electron_maxstep = 100
    /
    &ions
    /
    ATOMIC_SPECIES
    Al 26.981 Al.pbe-nl-kjpaw_psl.1.0.0.UPF
    Mg 24.305 Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF
    O 16.000 O.pbe-n-kjpaw_psl.1.0.0.UPF
    ATOMIC_POSITIONS crystal

    K_POINTS gamma

You need to rewrite the directory that contains the pseudopotentials, ``pseudo_dir``, and the pseudopotentials used in ``ATOMIC_SPECIES`` according to your environment. The pseudopotentials used in this sample can be downloaded from the following link.

- https://pseudopotentials.quantum-espresso.org/upf_files/Al.pbe-nl-kjpaw_psl.1.0.0.UPF
- https://pseudopotentials.quantum-espresso.org/upf_files/Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF
- https://pseudopotentials.quantum-espresso.org/upf_files/O.pbe-n-kjpaw_psl.1.0.0.UPF

In this example, ``calculation = 'relax'`` is used for structural optimization during the QE calculation, and ``gammma`` is used for ``K_POINTS`` to speed up the calculation. 

Preparation of the abICS file
===============================

Next, prepare an input file for abICS, referring to :ref:`Input file`.

Running the calculation
-----------------------

The sample scripts ``AL.sh`` and ``MC.sh`` are prepared to simplify the calculation procedure.
By performing these alternately, active learning and structural estimation by Monte Carlo are carried out.

Let's take a look at the contents of ``AL.sh`` first.
Note that before running these shell scripts, you need to change the permissions of ``run_pw.sh`` with ``chmod u+x``.
``run_pw.sh`` is a script to perform QE calculations and is called inside ``parallel_run.sh``, which will be described later.

.. code-block:: shell

    #!/bin/sh
    #SBATCH -p i8cpu
    #SBATCH -N 4
    #SBATCH -n 512
    #SBATCH -J spinel
    #SBATCH -c 1
    #SBATCH --time=0:30:00

    # Run reference DFT calc.
    echo start AL sample
    srun -n 8 abics_mlref input.toml >> active.out
    echo start parallel_run 1
    sh parallel_run.sh

    echo start AL final
    srun -n 8 abics_mlref input.toml >> active.out

    #train
    echo start training
    abics_train input.toml > train.out
    echo Done

The lines starting with ``#SBATCH`` and ``srun`` command are parameters of the job scheduler and the command to invoke parallel program (similar to ``mpiexec``) used on the ISSP supercomputer system B, respectively.
In this example, we are running an MPI parallel with 512 processes.
For more information about the job scheduler, please refer to the manuals of your machine.

.. code-block:: shell

    # Run reference DFT calc.
    echo start AL sample
    srun -n 8 abics_mlref input.toml >> active.out

The above code block generates an input file for ab initio calculation, which is the main source of the training data, using ``abics_mlref``.
At the first execution, the specified number of atomic arrangements are randomly generated, a separate directory is prepared for each atomic arrangement, and an input file is created in the directory.
At the same time, a file ``rundirs.txt`` is generated with the path of those directories.
This directory listing can be used to automate the execution of ab initio computation jobs for individual inputs.
We will then run the ab initio calculation based on the resulting file.

.. code-block:: shell

    echo start parallel_run 1
    sh parallel_run.sh

``parallel_run.sh`` is a script to run the QE exhaustive calculation using gnu parallel.
It will run the QE exhaustive calculation for the directories listed in rundirs.txt.
The results of the QE calculation will be stored in each directory.
Now that we have created the teacher data by the QE coverage calculation, we will move on to create the neural network potential in aenet.
First, we run ``abics_mlref`` again to create a file with the results of the ab initio calculations in a common format that abics_train will read.

.. code-block:: shell

    echo start AL final
    srun -n 8 abics_mlref input.toml >> active.out

Next, we use anet to create a neural network potential based on the training data.
The neural network potential is calculated by ``abics_train``.
The calculation is performed by reading the input file stored in ``base_input_dir`` in the ``[trainer]`` section of the input file.
When the calculation is completed successfully, the trained neural network is output to the baseinput directory.

.. code-block:: shell

    #train
    echo start training
    abics_train input.toml > train.out
    echo Done

The above process completes the AL.sh process for active learning.

Next, we use the trained neural network potential to find the optimization structure by abICS.
This process can be done in ``MC.sh``.
The following is the content of ``MC.sh``.

.. code-block:: shell

    #! /bin/sh
    #SBATCH -p i8cpu
    #SBATCH -N 1
    #SBATCH -n 8
    #SBATCH --time=00:30:00

    srun -n 8 abics_sampling input.toml >> aenet.out
    echo Done

Running abicsAL will create the ``MCxx`` directory (where xx is the number of runs).
With active learning in mind, additional functions have been implemented to obtain information such as the number of calculations by reading ``ALloop.progress``.
Under the ``MCxx`` directory, a folder will be created for the number of replicas.
Then, in these folders, the atomic arrangement (``structure.XXX.vasp``) for each step described in the VASP POSCAR file format, the atomic position given the lowest energy (``minE.vasp``), and each step temperature and energy (``obs.dat``) etc. are output.
For more details, please refer to the `abICS manual output file <https://issp-center-dev.github.io/abICS/docs/sphinx/ja/build/html/outputfiles/index.html>`_.

The results obtained by the above procedure depend on the accuracy of the neural network potential computed by aenet.
In the first step, we trained based on random configurations, thus the accuracy for low temperature structures is expected to be low.
Here, by repeating the step of calculating the energy again by first-principles calculation for the structure estimated by Monte Carlo and relearning it, we expect to improve the accuracy in the whole temperature range.
This process can be calculated by repeating ``AL.sh`` and ``MC.sh`` in turn.
The actual result of the calculation of the inversion rate (DOI) is shown in the figure below.
In this example, the first result is ``MC0``, followed by ``MC1``, ``MC2``, and so on.
The first run is quite different from the others, thus we can expect that it is not accurate.
On the other hand, if we train on the results of one Monte Carlo run, we find that the values are almost identical from the next run.

.. image:: ../../../image/DOI.*
   :width: 800px
   :align: center

In addition, DOI can be calculated by the following procedure.

1. Go to ``MCxxx``.

2. Create ``Tseparate`` directory by ``srun -n 8 abicsRXsepT ../input.toml``. (Align with the number of parallelism when ``abics_sampling`` is executed.
   In this tutorial, the number of parallelism is set to 8, so set it to 8.) 

3. copy ``calc_DOI.py`` and ``MgAl2O4.vasp`` in the sample directory.

4. Calculate the inversion rate for each temperature by ``srun -n 8 python3 calc_DOI.py ../input.toml``. (Align with the number of parallelism when ``abics_sampling`` is executed.
   In this tutorial, the number of parallelism is set to 8, so set it to 8.) 

