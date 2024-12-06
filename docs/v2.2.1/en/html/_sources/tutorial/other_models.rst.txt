.. _tutorial_nequip:

***********************************************
Sampling using other machine learning models
***********************************************

In abICS, in addition to the aenet, it is possible to perform sampling using 
other machine learning models such as NequIP, Allegro, and MLIP-3.
This section explains how to train and sample using each model.

Sampling with NequIP
----------------------------------------------

Installation of NequIP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use ``nequip``, you need to install NequIP.

Install it with the following command.

.. code-block:: bash

    $ python3 -m pip install wandb
    $ python3 -m pip install nequip

Also, when installing abICS, you can install NequIP by specifying the [nequip] option.

.. code-block:: bash

    $ cd /path/to/abics
    $ python3 -m pip install '.[nequip]'

Preparation of input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, prepare input_nequip.toml and set the parameters required to run NequIP.
Below, we extract [sampling.solver] and [train] with changes from the aenet input.

.. code-block:: toml
    
   [sampling.solver]
   type = 'nequip'
   base_input_dir = './baseinput_nequip'
   perturb = 0.0
   ignore_species = ["O"]

   [train]
   type = 'nequip'
   base_input_dir = './nequip_train_input'
   exe_command = { train = 'nequip-train' }
   ignore_species = ["O"]
   vac_map = []
   restart = false

Also, create the NequIP input file ``input.yaml`` in the ``nequip_train_input/train`` directory.

.. code-block:: yaml

   root: results/spinel
   run_name: run
   seed: 123
   dataset_seed: 456

   # network
   num_basis: 8
   BesselBasis_trainable: true
   PolynomialCutoff_p: 6
   l_max: 1
   r_max: 8.0
   parity: true
   num_layers: 3
   num_features: 16

   nonlinearity_type: gate

   nonlinearity_scalars:
     e: silu
     o: tanh

   nonlinearity_gates:
     e: silu
     o: tanh

   model_builders:
    - SimpleIrrepsConfig
    - EnergyModel
    - PerSpeciesRescale
    - RescaleEnergyEtc


   dataset: ase
   dataset_file_name: structure.xyz
   chemical_symbols:
     - Mg
     - Al

   # logging
   wandb: false
   # verbose: debug

   # training
   n_train: 80%
   n_val: 20%
   batch_size: 5
   train_val_split: random
   #shuffle: true
   metrics_key: validation_loss
   use_ema: true
   ema_decay: 0.99
   ema_use_num_updates: true
   max_epochs: 100
   learning_rate: 0.01
   # loss function
   loss_coeffs: total_energy

The procedure of model learning and sampling is the same as aenet.


Sampling with Allegro
----------------------------------------------

Models implemented as extensions of NequIP can be used as is by installing the extension package and setting the input file of NequIP appropriately. Allegro is one of the extension packages.

Installation of Allegro
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install Allegro with the following command.

.. code-block:: bash

    $ git clone --depth 1 https://github.com/mir-group/allegro.git
    $ cd allegro
    $ python3 -m pip install .


Preparation of input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, prepare input_allegro.toml and set the parameters required to run Allegro.   
Below, we extract ``[sampling.solver]`` and ``[train]`` with changes from the aenet input.

.. code-block:: toml
    
   [sampling.solver]
   type = 'allegro'
   base_input_dir = './baseinput_allegro'
   perturb = 0.0
   ignore_species = ["O"]

   [train]
   type = 'allegro'
   base_input_dir = './allegro_train_input'
   exe_command = {train = 'nequip-train'}
   ignore_species = ["O"]
   vac_map = []
   restart = false

Also, create the Allegro input file ``input.yaml`` in the ``allegro_train_input/train`` directory.

.. code-block:: yaml

   root: results/spinel
   run_name: run
   seed: 123
   dataset_seed: 456

   # network
   num_basis: 8
   BesselBasis_trainable: true
   PolynomialCutoff_p: 6
   l_max: 1
   r_max: 8.0
   parity: o3_full
   num_layers: 2

   env_embed_multiplicity: 16
   embed_initial_edge: true
   two_body_latent_mlp_latent_dimensions: [32, 64]
   two_body_latent_mlp_nonlinearity: silu
   latent_mlp_latent_dimensions: [64, 64]
   latent_mlp_nonlinearity: silu
   latent_mlp_initialization: uniform
   latent_resnet: true
   env_embed_mlp_latent_dimensions: []
   env_embed_mlp_nonlinearity: null
   env_embed_mlp_initialization: uniform
   edge_eng_mlp_latent_dimensions: [16]
   edge_eng_mlp_nonlinearity: null
   edge_eng_mlp_initialization: uniform

   model_builders:
    - allegro.model.Allegro
    - PerSpeciesRescale
    - RescaleEnergyEtc


   dataset: ase
   dataset_file_name: structure.xyz
   chemical_symbols:
     - Mg
     - Al

   # logging
   wandb: false
   # verbose: debug

   # training
   n_train: 80%
   n_val: 20%
   batch_size: 5
   train_val_split: random
   #shuffle: true
   metrics_key: validation_loss
   use_ema: true
   ema_decay: 0.99
   ema_use_num_updates: true
   max_epochs: 100
   learning_rate: 0.01
   # loss function
   loss_coeffs: total_energy

The procedure of model learning and sampling is the same as aenet.


Sampling with MLIP-3
----------------------------------------------

Installation of MLIP-3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use ``mlip-3``, you need to install MLIP-3.

Install it with the following command.

.. code-block:: bash

    $ git clone https://gitlab.com/ashapeev/mlip-3.git
    $ cd mlip-3
    $ ./configure --no-mpi
    $ make mlp


Preparation of input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, prepare ``input_mlip3.toml`` and set the parameters required to run MLIP-3.
Below, we extract ``[sampling.solver]`` and ``[train]`` with changes from the aenet input.

.. code-block:: toml
    
   [sampling.solver]
   type = 'mlip_3'
   path= '~/github/mlip-3/bin/mlp'
   base_input_dir = './baseinput'
   perturb = 0.0
   run_scheme = 'subprocess'
   ignore_species = ["O"]

   [train]
   type = 'mlip_3'
   base_input_dir = './mlip_3_train_input'
   exe_command = { train = '~/github/mlip-3/bin/mlp'}
   ignore_species = ["O"]
   vac_map = []
   restart = false

In the above, the ``path`` in ``[sampling.solver]`` and the ``exe_command`` in ``[train]``
specify the path to the MLIP-3 executable file ``mlp`` .
Please change them according to your environment.

Also, create the MLIP-3 input file ``input.almtp`` in the ``mlip_3_train_input/train`` directory.

.. code-block:: none

   MTP
   version = 1.1.0
   potential_name = MTP1m
   species_count = 3
   potential_tag = 
   radial_basis_type = RBChebyshev
    min_dist = 2.3
   	max_dist = 5
   	radial_basis_size = 8
	radial_funcs_count = 2
   alpha_moments_count = 8
   alpha_index_basic_count = 5
   alpha_index_basic = {{0, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}
   alpha_index_times_count = 5
   alpha_index_times = {{0, 0, 1, 5}, {1, 1, 1, 6}, {2, 2, 1, 6}, {3, 3, 1, 6}, {0, 5, 1, 7}}
   alpha_scalar_moments = 5
   alpha_moment_mapping = {0, 4, 5, 6, 7}


The procedure of model learning and sampling is the same as aenet.