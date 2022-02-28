.. pyMC documentation master file, created by
   sphinx-quickstart on Wed Jul 31 13:13:22 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Overview
------------------------------------------
abICS is a package for statistical thermodynamic calculations of atomic configurations on crystal lattices in multicomponent systems. It aims to enable quantitative prediction of temperature-dependent regular and irregular phase transitions and changes in short- and intermediate-range order in alloy and composite oxide systems, as well as prediction of physical properties taking these into account. abICS was originally developed with the intention of directly combining first-principles calculations with replica-exchange Monte Carlo methods to perform the kind of calculations described above, but the scale of the models and the number of steps that can be calculated are limited by the large computational cost of first-principles calculations. In contrast, Ver. 2 implements an active learning method to construct a neural network model that can rapidly predict the energy after structural optimization, dramatically improving the sampling speed `[preprint] <https://arxiv.org/abs/2008.02572>`_ . In this tutorial, we will explain how to calculate the Mg/Al site inversion of MgAl2O4 spinel crystals using this method. We will use the examples in ``examples/standard/active_learning``. Please note that abICS is basically developed for use on supercomputers or medium-sized computer clusters, and requires a certain number of CPU cores (at least several hundred) to run.

The general flow of the active learning method implemented in abICS is as follows.

1. Perform ab initio calculations on a large number of randomly generated atomic configurations and prepare training data (correspondence between configurations and energies).
2. Build a neural network model that predicts energy from atomic configurations using the prepared training data.
3. Perform statistical thermodynamic sampling of atomic configurations using a replica exchange Monte Carlo method with a neural network model.
4. Evaluate the accuracy of the neural network model by sampling the ion configurations that appear in the Monte Carlo calculations and performing ab initio calculations on each of them.
5. If the accuracy is not sufficient, add the results calculated in 4. to the training data and repeat from 2.

In ::ref:`sec_basic_usage`, we will explain the procedure for using abICS in individual steps.


.. image:: ../../../image/schmatic_AR.png
   :width: 800px
   :align: center


