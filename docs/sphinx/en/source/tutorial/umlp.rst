Universal machine-learning potentials
=====================================

The example in ``examples/active_learning_umlp/`` trains an on-lattice model
for MgAl2O4 spinel from reference energies calculated with a pretrained Orb
potential. TorchSim relaxes several structures together on one NVIDIA GPU.
An ASE alternative supports sequential CPU or GPU relaxation.

Follow the `step-by-step example guide
<https://github.com/issp-center-dev/abICS/tree/develop/examples/active_learning_umlp>`_
for installation, input preparation, two active learning cycles, and interpretation
of the degree of inversion. The optional GPU packages are installed separately
from abICS and require Python 3.12 or newer.
