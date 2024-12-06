.. highlight:: none

abicsRXsepT
-------------------------------

This tool is for reordering the resulting structures and energies at each sampling step of
a RXMC run by temperature. It is used after an abICS RXMC run is finished as::

   $ mpiexec -np NPROCS abicsRXsepT input.toml NSKIP

``NPROCS`` should be equal to or larger than the number of replicas, and ``input.toml`` 
should be replaced by the abICS input file that was used for this run.
``NSKIP`` is an optional parameter and is used for specifying the number of thermalization steps to skip when calculating the energy averages at each temperature (default value is 0).
The results are stored in the ``Tseparate`` directory, and energy averages vs. temperature are stored in ``Tseparate/energies_T.dat``
