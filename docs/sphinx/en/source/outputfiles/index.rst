.. highlight:: none

***************************
Output Files Format
***************************

``RANK`` means the rank of process (replica) (``0, 1, ...``).

``RANK/structure.XXX.vasp``
==============================
The atomic coordinates for each step are saved in the POSCAR file format of VASP.
``XXX`` in the filename corresponds to the index of the step.

Example::

 Mg8 Al16 O32
 1.0
 8.113600 0.000000 0.000000
 0.000000 8.113600 0.000000
 0.000000 0.000000 8.113600
 Al Mg O
 16 8 32
 direct
 0.011208 0.995214 0.998158 Al
 0.758187 0.240787 0.499981 Al
 ... skipped ...
 0.746308 0.744706 0.233021 O
 0.257199 0.255424 0.771040 O

``RANK/minE.vasp``
====================
The lowest-energy structure among the samples in this replica.

``RANK/obs.dat``
===================
The temperature and the total energy for each step in units of eV.

Example::

 0	0.1034076	-41690.28269769395
 1	0.1034076	-41692.06763035158
 2	0.1034076	-41692.06763035158
 3	0.1034076	-41691.98205990787
 4	0.1034076	-41692.74143710456

``RANK/obs_save.npy``
========================
The total energy for each step in units of eV in the Numpy binary format.
Users can load it as ``darray`` by using ``numpy.load('obs_save.npy')``.

Example::

 $ python -c "import numpy; print(numpy.load('obs_save.npy'))"
 [[-41690.28269769]
  [-41692.06763035]
  [-41692.06763035]
  [-41691.98205991]
  [-41692.7414371 ]]

``RANK/kT_hist.npy``
=======================
The temperature for each step in units of eV in the Numpy binary format.
Users can load it as ``darray`` by using ``numpy.load('kT_hist.npy')``.

Example::

 $ python -c "import numpy; print(numpy.load('kT_hist.npy'))"
 [0.1034076 0.1034076 0.1034076 0.1034076 0.1034076]


``RANK/Trank_hist.npy``
=======================
(ONLY for RXMC)
The rank (index) of the temperature for each step in the Numpy binary format.
Users can load it as ``darray`` by using ``numpy.load('Trank_hist.npy')``.

Example::

 $ python -c "import numpy; print(numpy.load('Trank_hist.npy'))"
 [1 1 1 1 1]

``RANK/logweight_hist.npy``
=============================
(ONLY for PAMC)
The logarithm of the Neal-Jarzynski weight for each step in the Numpy binary format.

Example::

 $ python -c "import numpy; print(numpy.load('logweight_hist.npy'))"
 [0 0 0 0 0]


``RANK/acceptance_ratio.dat``
===============================
Acceptance ratio of Monte Carlo steps for each temperature.
The first column is temperature and the second column is acceptance ratio (number of accepted / number of trials).

``logZ.dat``
==============
The logarithm of the partition function, :math:`\log Z_i/Z_0` where :math:`i` is the index of temperature.

- The 1st column is temperature :math:`T_i`.
- The 2nd and 3rd columns are :math:`\log Z_i/Z_0` and its error.
- The 4th and 5th columms are :math:`\log Z_i/Z_{i-1}` and its error.

``<name>.dat``
===============
Canonical expectation value :math:`\langle O \rangle` and statistical error :math:`\sigma[O]` of an observable :math:`O` for each temperature.
``<name>`` is the name of the observable which is specified by ``name`` keyword in ``[[observer.solver]]`` section of the input file.

- The 1st column is temperature :math:`T_i`.
- The 2nd and 3rd columns are :math:`\langle O \rangle` and its error.
- The 4th and 5th columns are :math:`\langle O^2 \rangle` and its error.
- The 6th and 7th columns are fluctuation, :math:`\langle O^2 \rangle - \langle O \rangle^2` and its error.

   - Note that the heat capacity :math:`C` is related to the fluctuation of energy as :math:`k_B T^2 C = \left[ \langle E^2 \rangle - \langle E \rangle^2 \right]`.


Output files of active learning (``abics_mlref``)
==================================================

``abics_mlref`` writes the input and output of the :math:`n`-th active-learning step to the directory ``ALn`` (``AL0``, ``AL1``, ...).
``AL0`` holds randomly generated configurations; ``AL1`` and later hold configurations taken from the preceding Monte Carlo run ``MC(n-1)``.
``RANK`` is the index of the replica given by ``nreplicas`` in the ``[mlref]`` section.

``ALloop.progress``
--------------------------
Progress of the active-learning loop.
``abics_mlref``, ``abics_train`` and ``abics_sampling`` append ``AL0``, ``train``, ``MC0``, ``AL1``, ... when they finish a step,
and use it to decide which program should run next.

``rundirs.txt``
--------------------------
List of the input directories (``ALn/RANK/inputXX/baseinputK``) in which the first-principles solver has to be run.
It is rewritten every time ``abics_mlref`` generates input files.

``ALn/RANK/inputXX/``
--------------------------
Input and output of the solver for configuration ``XX``.
In ``AL0``, ``XX`` is a running number; in ``AL1`` and later it is the index of the sample of the Monte Carlo run
the configuration was taken from (``MC(n-1)/RANK/structure.XX.vasp``).

- ``initial.vasp``: the atomic configuration before the perturbation (``perturb`` of ``[mlref.solver]``) is applied
- ``baseinputK/``: input and output files of the ``K``-th stage of the solver (one directory per stage when several directories are given in ``base_input_dir`` of ``[mlref.solver]``)

``ALn/RANK/structure.XX.vasp``
-----------------------------------
The atomic configuration after the solver run (structure relaxation). These are the training data of ``abics_train``.

``ALn/RANK/energy_corr.dat``
-----------------------------------
For each extracted configuration, the energy predicted by the machine-learning model and the reference energy of the first-principles solver (in eV per supercell).
This file is meant for checking the accuracy of the model: a scatter plot of the predicted energy against the reference energy shows the accuracy as the deviation from the diagonal.

- The 1st column is the energy predicted by the model (the value during the Monte Carlo run). It is ``NaN`` in ``AL0``, where no model exists yet.
- The 2nd column is the reference energy of the solver.
- The 3rd column is the index ``XX`` of the configuration.

Example::

 -93.40812444 -93.1635694601349 0
 -93.48409204 -93.32726602032281 8
 -93.7486019 -93.81712285031085 17

``ALn/RANK/relax_max.dat``
-----------------------------------
(``AL1`` and later) For each configuration, the largest displacement of an atom between ``initial.vasp`` and the structure after the solver run.
The 1st column is the index of the atom with the largest displacement and the 2nd column is the displacement in Angstrom.
The perturbation ``perturb`` is included.

``ALn/convergence.dat``
-----------------------------------
(``AL1`` and later) Error statistics of the machine-learning model in this step:
the differences of the two columns of ``energy_corr.dat`` summarised per replica and over all replicas.
The same table is written to the log of ``abics_mlref``.
It can be used to decide whether to continue the active-learning loop.

- The 1st column is the replica index (the last row, ``all``, is over all replicas).
- The 2nd column is the number of configurations.
- The 3rd column is the root mean square error (RMSE).
- The 4th column is the mean absolute error (MAE).
- The 5th column is the largest absolute error.
- The 6th column is the index ``XX`` of the configuration with the largest error.

The errors are in the units of the energies (eV per supercell).

Example::

 # AL1: errors of the model energies against the reference solver
 # (per supercell, in the units of energy_corr.dat)
 # replica  n  rmse  mae  max_error  max_error_sample
 0  10  0.256721  0.210012  0.465759  61
 1  10  0.213205  0.140375  0.513428  0
 all  20  0.235968  0.175193  0.513428  0

``ALn/parity.png``
-----------------------------------
(``AL1`` and later, only when matplotlib is installed)
Parity plot of ``energy_corr.dat``: the energy predicted by the machine-learning model against the reference energy,
one colour per replica, with the diagonal as a guide.

``ALn/baseinput.progress``
-----------------------------------
(When several ``base_input_dir`` are given, staged mode only) The indices of the solver stages that have been completed.

``ALn/libE/``
-----------------------------------
(``exec_mode = "integrated"`` only) The logs of libEnsemble, which manages the solver runs
(``ensemble.log``, ``libE_stats.txt``, ...). The standard output and error of each solver run are written to
``ALn/RANK/inputXX/baseinputK/stdout`` and ``stderr``, and ``abics_mlref.done`` is written into a stage directory
once the solver exited successfully and its output was parsed (on a rerun such stages are not calculated again).
