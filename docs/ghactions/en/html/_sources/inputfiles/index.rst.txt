.. _input_format:

***************************
Input Files Format
***************************

The input file of abICS is constructed by the following five sections:

1. [sampling] section specifies the parameters of the replica exchange Monte Carlo part, such as the number of replicas, the temperature range, and the number of Monte Carlo steps. In addition, [sampling.solver] subsection specifies the parameters for the (first principle calculation) solver, including the type of solver (VASP, QE,...), the path to the solver, and the directory containing immutable input files.

2. [mlref] section specifies options for extracting only atomic configurations from the sampling results in order to evaluate the accuracy of the neural network model and to expand the training data. In addition, for generating training data, [mlref.solver] subsection specifies the parameters for the (first principle calculation) solver, including the type of solver (VASP, QE,...), the path to the solver, and the directory containing immutable input files. This section is used for ``abics_mlref`` .

3. [train] section specifies optinons for making a trainer to learn a placement energy prediction model from training data.  This section is used for ``abics_train`` .

4. [observer] section specifies the type of physical quantity to be calculated.

5. [config] section specifies the configuration of the alloy, etc.

The following sections describe the detail of each section.

.. toctree::
   :maxdepth: 1

   parameter_sampling
   parameter_mlref
   parameter_train
   parameter_observer
   parameter_config
   
   
