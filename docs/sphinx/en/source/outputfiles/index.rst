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


``logZ.dat``
==============
(ONLY for PAMC)
The logarithm of the partition function, :math:`\log Z_i/Z_0` where :math:`i` is the index of temperature.
The 1st column is temperature :math:`T_i`.
The 2nd column is :math:`\log Z_i/Z_0`.
The 3rd colum is :math:`\log Z_i/Z_{i-1}`.

``result.dat``
===============
(ONLY for PAMC)
Canonical expectation value :math:`\langle O \rangle` and statistical error :math:`\sigma[O]` of observables :math:`O` for each temperature.

- The 1st column is temperature :math:`T_i`.
- The 2nd and 3rd columns are energy :math:`\langle E \rangle` and its error.
- The 4th and 5th columns are squared energy :math:`\langle E^2 \rangle` and its error.
- The 6th and 7th columns are fluctuation of energy :math:`\langle E^2 \rangle - \langle E \rangle^2` and its error.

   - Note that the heat capacity :math:`C` is related to the fluctuation of energy as :math:`k_B T^2 C = \left[ \langle E^2 \rangle - \langle E \rangle^2 \right]`.

For Potts solver, total magnetization :math:`\langle M \rangle = \langle \sum_i \delta_{\sigma_i,0} - 1/Q \rangle`, squared magnetization :math:`\langle M^2 \rangle`, and fluctuation :math:`\langle M^2 \rangle - \langle M \rangle^2` will be saved as 8th and the following columns.
