.. _tutorial_nequip:

*************************************
他のモデルを利用したサンプリング
*************************************

abICSでは、機械学習モデルとして、aenet以外にも、
NequIP, Allegro, MLIP-3を利用したサンプリングが可能となっています。
本項では、それぞれのモデルの学習およびサンプリングの方法について説明します。

NequIPを利用したサンプリング
----------------------------------------------

NequIP のインストール
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``nequip`` の利用には、 NequIPのインストールが必要です。

下記コマンドにてインストールします。

.. code-block:: bash

    $ pip3 install wandb
    $ pip3 install nequip

また、abICSインストール時に[nequip]オプションを指定すれば、NequIPもインストールされます。

.. code-block:: bash

    $ cd /path/to/abics
    $ pip3 install .abics[nequip]

インプットファイルの準備
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

まず、input_nequip.tomlを準備し、NequIPの実行に必要なパラメータを設定します。
下では、aenetのインプットから変更のある[sampling.solver]と[train]を抜粋しています。

.. code-block:: toml
    
   [sampling.solver]
   type = 'nequip'
   base_input_dir = './baseinput_nequip'
   perturb = 0.0
   #run_scheme = 'subprocess' #'mpi_spawn_ready'
   ignore_species = ["O"]

   [train]
   type = 'nequip'
   base_input_dir = './nequip_train_input'
   exe_command = ['', 'nequip-train']
   ignore_species = ["O"]
   vac_map = []
   restart = false

また、NequIPのインプットファイルinput.yamlをnequip_train_input/trainディレクトリに作成します。

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
   n_train: 70
   n_val: 10
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

モデル学習、サンプリングの方法に関してはaenetと同様です。


Allegroを利用したサンプリング
----------------------------------------------

Allegro のインストール
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``allegro`` の利用には、Allegroのインストールが必要です。

下記コマンドにてインストールします。

.. code-block:: bash

    $ git clone --depth 1 https://github.com/mir-group/allegro.git
    $ cd allegro
    $ pip3 install .


インプットファイルの準備
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

まず、input_allegro.tomlを準備し、Allegroの実行に必要なパラメータを設定します。
下では、aenetのインプットから変更のある[sampling.solver]と[train]を抜粋しています。

.. code-block:: toml
    
   [sampling.solver]
   type = 'allegro'
   base_input_dir = './baseinput_allegro'
   perturb = 0.0
   #run_scheme = 'subprocess' #'mpi_spawn_ready'
   ignore_species = ["O"]

   [train]
   type = 'allegro'
   base_input_dir = './allegro_train_input'
   exe_command = ['', 'nequip-train']
   ignore_species = ["O"]
   vac_map = []
   restart = false

また、Allegroのインプットファイルinput.yamlをallegro_train_input/trainディレクトリに作成します。

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
   num_features: 16

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
   n_train: 70
   n_val: 10
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

モデル学習、サンプリングの方法に関してはaenetと同様です。


MLIP-3を利用したサンプリング
----------------------------------------------

MLIP-3 のインストール
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``mlip-3`` の利用には、 MLIP-3のインストールが必要です。

下記コマンドにてインストールします。

.. code-block:: bash

    $ git clone https://gitlab.com/ashapeev/mlip-3.git
    $ cd mlip-3
    $ ./configure --no-mpi
    $ make mlp


インプットファイルの準備
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

まず、input_mlip3.tomlを準備し、mlip-3の実行に必要なパラメータを設定します。
下では、aenetのインプットから変更のある[sampling.solver]と[train]を抜粋しています。

.. code-block:: toml
    
   [sampling.solver]
   type = 'mlip_3'
   path= '~/github/mlip-3/bin/mlp'
   base_input_dir = './baseinput'
   perturb = 0.0
   run_scheme = 'subprocess' #'mpi_spawn_ready'
   ignore_species = ["O"]

   [train]
   type = 'mlip_3'
   base_input_dir = './mlip_3_train_input'
   exe_command = ['~/github/mlip-3/bin/mlp','~/github/mlip-3/bin/mlp']
   ignore_species = ["O"]
   vac_map = []
   restart = false

上記の内、[sampling.solver]のpathと[train]のexe_commandのリストでは
MLIP-3の実行ファイルmlpのパスを指定します。お使いの環境に合わせて変更してください。

また、MLIP-3のインプットファイルinput.almtpをmlip_3_train_input/trainディレクトリに作成します。

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

モデル学習、サンプリングの方法に関してはaenetと同様です。