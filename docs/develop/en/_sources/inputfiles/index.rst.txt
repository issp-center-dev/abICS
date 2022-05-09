.. _input_format:

***************************
Input Files Format
***************************

The input file of abICS is constructed by the following four sections:

1. [replica] section specifies the parameters of the replica exchange Monte Carlo part, such as the number of replicas, the temperature range, and the number of Monte Carlo steps.
  
2. [solver] section specifies the parameters for the (first principle calculation) solver, including the type of solver (VASP, QE,...), the path to the solver, and the directory containing immutable input files.
   
3. [observer] section specifies the type of physical quantity to be calculated.

4. [config] section specifies the configuration of the alloy, etc.

The following sections describe the detail of each section.

.. toctree::
   :maxdepth: 1

   parameter_replica
   parameter_solver
   parameter_observer
   parameter_config
   
   
