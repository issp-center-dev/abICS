.. highlight:: none

[train] section
-------------------------------

We create neural network potentials by learning the results of ab initio calculations made by abics\_train using aenet.
The input information for abics_train is described in the ``[trainer]`` section. The description of each parameter is as follows.

This section specifies the type of physical quantity to be acquired.
An example is shown as follows:

  ::

     [trainer] # Configure the model trainer.
     type = 'aenet'
     base_input_dir = '. /aenet_train_input'
     exe_command = ['~/git/aenet/bin/generate.x-2.0.4-ifort_serial', 
                  'srun ~/git/aenet/bin/train.x-2.0.4-ifort_intelmpi']
     ignore_species = ["O"]

Input Format
^^^^^^^^^^^^^

Keywords and their values are specified by a keyword and its value in the form ``keyword = value``.
Comments can also be entered by adding # (Subsequent characters are ignored).

Key words
^^^^^^^^^^

- ``type``

  **Format :** str

  **Description :** The trainer to generate the neural network potential (currently only 'aenet').
  
- ``base_input_dir``

  **Format :** str 

  **Description :**
  Path of the directory containing the input files that the learner refers to.

- ``exe_command``

  **Format :** list of str 

  **Description :**
  List of commands to execute; if you use aenet, you need to specify the path to ``generate.x`` and ``train.x``.
  
- ``ignore_species``

  **Format :** list

  **Description :**
  Same as ``ignore_species`` in [sampling.solver] section. Specify atomic species to "ignore" in neural network models such as ``aenet``. For those that always have an occupancy of 1, it is computationally more efficient to ignore their presence when training and evaluating neural network models.
 
