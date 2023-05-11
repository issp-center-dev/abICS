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
