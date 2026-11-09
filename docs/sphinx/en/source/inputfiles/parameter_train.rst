.. highlight:: none

[train] section
-------------------------------

``abics_train`` creates and trains a regression model from configurations to energies.
Indeed, ``abics_train`` uses an external program to train the model.
In the current version, Aenet, Nequip, and MLIP-3 are supported as an external program.
For software-specific notes (such as input file names), see :ref:`trainer_specific_notes`.

The input information for ``abics_train`` is described in the ``[train]`` section. The description of each parameter is as follows.
An example is shown as follows:

  ::

     [train] # Configure the model trainer.
     type = 'aenet'
     base_input_dir = './aenet_train_input'
     ignore_species = ["O"]
     [train.exe_command]
     generate = '~/git/aenet/bin/generate.x-2.0.4-ifort_serial'
     train = 'srun ~/git/aenet/bin/train.x-2.0.4-ifort_intelmpi'

Input Format
^^^^^^^^^^^^^

Keywords and their values are specified by a keyword and its value in the form ``keyword = value``.
Comments can also be entered by adding # (Subsequent characters are ignored).

Key words
^^^^^^^^^^

- ``type``

  **Format :** str

  **Description :** The trainer to generate the neural network potential (currently 'aenet', 'nequip', and 'mlip_3' are available).
  
- ``base_input_dir``

  **Format :** str 

  **Description :**
  Path of the directory containing the input files that the learner refers to.

- ``exe_command``

  **Format :** dict

  **Description :**
  List of commands to execute; if you use aenet, you need to specify the path to ``generate.x`` and ``train.x``.

  - ``type = 'aenet'``

    - ``generate`` and ``train`` keys are required.
    - ``generate``

        - Specify the path to ``generate.x`` of aenet.

    - ``train``

        - Specify the path to ``train.x`` of aenet.
        - The MPI parallel version is available. In that case, set the command to execute MPI (e.g., ``srun``, ``mpirun``) .
    
    - Array format is supported for compatibility with abICS 2.0 and earlier.
      The first element is ``generate``, and the second element is ``train``.

  - ``type = 'nequip'``

    - ``train`` 

        - Specify the path to ``nequip-train``.

  - ``type = 'mlip_3'``

    - ``train``

        - Specify the path to ``mlp``.

  
- ``ignore_species``

  **Format :** list

  **Description :**
  Same as ``ignore_species`` in [sampling.solver] section. Specify atomic species to "ignore" in neural network models such as ``aenet``. For those that always have an occupancy of 1, it is computationally more efficient to ignore their presence when training and evaluating neural network models.

- ``energy_ref``

  **Format :** float (default: 0.0)

  **Description :**
  Reference energy per unit cell. Before training, this value times the supercell multiplicity
  (volume of the structure / volume of the ``[config]`` unit cell) is subtracted from every training energy,
  so that a large configuration-independent constant such as the energy of the host lattice is removed and
  the model learns the configuration-dependent part. ``abics_mlref`` adds the same amount back to the model
  energies it writes to ``energy_corr.dat``, so they remain comparable with the reference energies.

- ``previous_dirs``

  **Format :** list (directories, default: empty)

  **Description :**
  Directories with training data in addition to the ``AL*`` steps. Each directory must contain only
  ``structure.0.xsf``, ``structure.1.xsf``, ... in aenet's XSF format with the energy on the first line
  (the number of structures is taken from the number of files). They are added to the data set before the
  ``AL*`` data, and ``abics_train`` can be run on them alone when ``ALloop.progress`` does not exist yet
  (before any ``abics_mlref`` run).

- ``prev_dirs_energy_ref``

  **Format :** bool (default: false)

  **Description :**
  Whether ``energy_ref`` is applied to the data of ``previous_dirs`` as well. If true, the multiplicity of each
  structure is derived from its volume and ``energy_ref`` is subtracted (the structures must be supercells of
  the ``[config]`` unit cell); if false, the energies of ``previous_dirs`` are used as they are.

