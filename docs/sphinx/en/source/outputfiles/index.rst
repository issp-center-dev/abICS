.. highlight:: none

***************************
Output Files Format
***************************

The calculation results are output in each replica directory.

``structure.XXX.vasp``
=========================
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

``minE.vasp``
====================
The lowest-energy structure among the samples in this replica.

``obs.dat``
===================
The temperature and the total energy for each step in units of eV.

Example::

 0	0.1034076	-41690.28269769395
 1	0.1034076	-41692.06763035158
 2	0.1034076	-41692.06763035158
 3	0.1034076	-41691.98205990787
 4	0.1034076	-41692.74143710456

``obs_save.npy``
==================
The total energy for each step in units of eV in the Numpy binary format.
Users can load it as ``darray`` by using ``numpy.load('obs_save.npy')``.

Example::

 $ python -c "import numpy; print(numpy.load('obs_save.npy'))"
 [[-41690.28269769]
  [-41692.06763035]
  [-41692.06763035]
  [-41691.98205991]
  [-41692.7414371 ]]

``kT_hist.npy``
==================
The temperature for each step in units of eV in the Numpy binary format.
Users can load it as ``darray`` by using ``numpy.load('kT_hist.npy')``.

Example::

 $ python -c "import numpy; print(numpy.load('kT_hist.npy'))"
 [0.1034076 0.1034076 0.1034076 0.1034076 0.1034076]


``Trank_hist.npy``
==================
The rank (index) of the temperature for each step in the Numpy binary format.
Users can load it as ``darray`` by using ``numpy.load('Trank_hist.npy')``.

Example::

 $ python -c "import numpy; print(numpy.load('Trank_hist.npy'))"
 [1 1 1 1 1]

