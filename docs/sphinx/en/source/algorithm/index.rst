**********
Algorithm
**********
abICS is designed for combining parallel extended ensemble methods with
arbitrary energy calculators. At present, only the replica exchange
Monte Carlo method is implemented.

Replica exchange Monte Carlo method
------------------------------------
A disadvantage of the widely-used Metropolis Monte Carlo algorithm is
that it tends to get stuck in local minima.
The replica exchange approach aims to overcome this problem by
considering multiple copies, or replicas, of the system under study.
The algorithm may be described roughly as follows
(see references below for more accurate descriptions).
Monte Carlo sampling is performed on each replica independently at
varying temperatures. At preset intervals, the temperatures are
exchanged according to a Metropolis criterion that essentially
assigns lower temperatures to replicas that happen to have lower
energies. This allows an efficient sampling of the global configuration
space using replicas at higher temperatures and accurate sampling of
the local energy landscape at lower temperatures.

In abICS, parameters related to the replica exchange Monte Carlo method are specified in the ``[replica]`` section of the input file.
By setting the lower limit of the replica temperature to ``kTstart`` , the upper limit to ``kTend``, and the number of replicas to ``nreplicas``,
``nreplicas`` replica systems with different temperatures ( :math:`T_0, T_1, \cdots T_{\verb|nreplicas|-1}`) are sampled, where

.. math::
   
   T_i = \frac{\bf{kTend}-\bf{kTstart}}{\bf{nreplicas}-1} i + \bf{kTstart}.

In abICS, using ``nprocs_per_replica``, the number of parallel solver processes that performs the calculation on each replica can be specified.
The number of Monte Carlo steps is specified by ``nsteps``, and the exchange transition probability :math:`R` for each ``RXtrial_frequency`` step is defined as

.. math::

   R = \exp\left[-\left(\frac{1}{T_i}-\frac{1}{T_{k}}\right)\left(E(X_i)-E(X_{k})\right)\right],

where  :math:`X_i` is the state for :math:`i` -th replica system. In abICS, the exchange transition is tried between replicas with adjacent temperatures.
The temperature exchange :math:`T_i \leftrightarrow T_{k}` is performed with the exchange transition probability :math:`R`.
Physical quantities such as the total energy is measured at each ``sample_frequency`` step.

- About replica exchange Monte Carlo method

  - `K. Hukushima and K. Nemoto, J. Phys. Soc. Japan, 65, 1604 (1996) <https://journals.jps.jp/doi/abs/10.1143/JPSJ.65.1604>`_.
  - `R. Swendsen and J. Wang, Phys. Rev. Lett. 57, 2607 (1986) <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.57.2607>`_.


About configuration and update
------------------------------------

Here, the outline of the definition of the configuration in abICS and the update by the Monte Carlo method are explained using :numref:`alg_sampling` as an example.

(a)-(c) are schematic figures of ``unitcell``, ``base_structure``, and ``defect_structure``, where blue, green, and black circles are the atomic types defined by ``base_structure``, respectively. The star symbol indicates the location of the defects defined by ``defect_structure``.
(d) is a schematic figure for specifying the atomic species in ``base_structure``. Here, three atomic species of blue, green, and black are defined. How each atom is arranged is defined by ``coords`` for each atom type.
(e) is a schematic figure for specifying the group of atoms to be located at the defect position with ``defect_structure``. orange defines a group consisting of four atoms composed of two types of atoms, and purple forms a group of three atoms composed of three types of atoms. These groups are placed at defect points specified by ``defect_structure.coords``. The arrangement of atoms in each group can be specified by ``coords`` in the ``defect_structure.groups`` section.
``defect_structure`` can be defined multiple times. Groups of each ``defect_structure`` will be placed at points of each one.
(f) is a schematic figure about the update of the Monte Carlo method. In the update, there are two patterns, one that swaps two atom groups of different type
, and the other that changes the orientation within the atom group without changing the arrangement. The type of updates is automatically selected with 1/2 probability. The energy is calculated with the specified solver from the proposed configuration :math:`X_ {trial}` and then the adoption rate :math:`P (X_i \rightarrow X_ {trial})` is calculated.

.. figure:: ../../../image/alg_sampling.png
     :name: alg_sampling
     :scale: 15%
	    
     (a)-(e) Definition of lattice in abICS. (f) A schematic of MonteCarlo method. Details are described in the text.



- Overview of abICS

  - `S. Kasamatsu and O. Sugino, J. Phys. Condens. Matter, 31, 085901 (2019) <https://iopscience.iop.org/article/10.1088/1361-648X/aaf75c/meta>`_.


Active learning
------------------------------------------
abICS was originally developed with the intention of directly combining first-principles calculations with replica-exchange Monte Carlo methods to perform the kind of calculations described above,
but the scale of the models and the number of steps that can be calculated are limited by the large computational cost of first-principles calculations.
In contrast, Ver. 2 implements an active learning method to construct a neural network model that can rapidly predict the energy after structural optimization,
dramatically improving the sampling speed `[preprint] <https://arxiv.org/abs/2008.02572>`_ .

The general flow of the active learning method implemented in abICS is as follows.

1. Perform ab initio calculations on a large number of randomly generated atomic configurations and prepare training data (correspondence between configurations and energies).
2. Build a neural network model that predicts energy from atomic configurations using the prepared training data.
3. Perform statistical thermodynamic sampling of atomic configurations using a replica exchange Monte Carlo method with a neural network model.
4. Evaluate the accuracy of the neural network model by sampling the ion configurations that appear in the Monte Carlo calculations and performing ab initio calculations on each of them.
5. If the accuracy is not sufficient, add the results calculated in 4. to the training data and repeat from 2.

.. image:: ../../../image/schmatic_AR.png
   :width: 800px
   :align: center



