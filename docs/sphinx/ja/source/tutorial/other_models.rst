.. _tutorial_nequip:

*************************************
他のモデルを利用したサンプリング
*************************************

abICSでは、機械学習モデルとして、aenet以外にも、
NequIP, MLIP-3を利用したサンプリングが可能となっています。
本項では、それぞれのモデルの学習およびサンプリングの方法について説明します。

NequIPを利用したサンプリング
----------------------------------------------

NequIP のインストール
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``nequip`` の利用には、 NequIPのインストールが必要です。

下記コマンドにてインストールします。

.. code-block:: bash

    $ python3 -m pip install wandb
    $ python3 -m pip install nequip

また、abICSインストール時に[nequip]オプションを指定すれば、NequIPもインストールされます。

.. code-block:: bash

    $ cd /path/to/abics
    $ python3 -m pip install '.abics[nequip]'

インプットファイルの準備
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

まず、input_nequip.tomlを準備し、NequIPの実行に必要なパラメータを設定します。
下では、aenetのインプットから変更のある[sampling.solver]と[train]を抜粋しています。

.. code-block:: toml
    
   [sampling.solver]
   type = 'nequip'
   base_input_dir = './baseinput_nequip'
   perturb = 0.0
   ignore_species = ["O"]

   [train]
   type = 'nequip'
   base_input_dir = './nequip_train_input'
   exe_command = { train = 'nequip-train'}
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

モデル学習、サンプリングの方法に関してはaenetと同様です。


Allegroを利用したサンプリング
----------------------------------------------

NequIPの拡張として実装されたモデルも、拡張パッケージをインストールし、NequIPの入力ファイルを適切に設定することで、そのまま利用可能です。Allegroは拡張パッケージの一つです。

Allegro のインストール
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

下記コマンドにてインストールします。

.. code-block:: bash

    $ git clone --depth 1 https://github.com/mir-group/allegro.git
    $ cd allegro
    $ python3 -m pip install .


インプットファイルの準備
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

まず、input_allegro.tomlを準備し、Allegroの実行に必要なパラメータを設定します。
下では、aenetのインプットから変更のある[sampling.solver]と[train]を抜粋しています。

.. code-block:: toml
    
   [sampling.solver]
   type = 'allegro'
   base_input_dir = './baseinput_allegro'
   perturb = 0.0
   ignore_species = ["O"]

   [train]
   type = 'allegro'
   base_input_dir = './allegro_train_input'
   exe_command = { train =  'nequip-train' }
   ignore_species = ["O"]
   vac_map = []
   restart = false

また、Allegroのインプットファイル ``input.yaml`` を ``allegro_train_input/train`` ディレクトリに作成します。

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
   # num_features: 16

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
   path= '~/git/mlip-3/bin/mlp'
   base_input_dir = './baseinput'
   perturb = 0.0
   run_scheme = 'subprocess'
   ignore_species = ["O"]

   [train]
   type = 'mlip_3'
   base_input_dir = './mlip_3_train_input'
   exe_command = { train = '~/git/mlip-3/bin/mlp'}
   ignore_species = ["O"]
   vac_map = []
   restart = false

上記の内、 ``[sampling.solver]`` の ``path`` と ``[train]`` の ``exe_command`` では
MLIP-3の実行ファイル ``mlp`` のパスを指定します。お使いの環境に合わせて変更してください。

また、MLIP-3のインプットファイル ``input.almtp`` を ``mlip_3_train_input/train`` ディレクトリに作成します。

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

SevenNetを利用したサンプリング
----------------------------------------------

SevenNetのインストール
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``sevennet`` の利用には、 SevenNetのインストールが必要です。

下記コマンドにてインストールします。

.. code-block:: bash

    $ python3 -m pip install sevenn

学習済みモデルの利用
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SevenNetでは、モデルを学習してからサンプリングを行う事以外に、
学習済みモデルを利用してサンプリングを行う事も可能です。

学習済みモデルを用いる場合は、 ``[sanmping.solver]`` セクションを下記のように設定します。

.. code-block:: toml

   [sampling.solver]
   type = 'sevennet'
   perturb = 0.0

サンプリングの方法については、aenetと同様です。

モデル学習から実行する場合
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

モデルの学習から実行する場合は、 ``[train]`` セクションを適切に設定の上で、
``[sanmping.solver]`` セクションに、
``use_pretrained = false`` を追加する必要があります。
``relax = false`` とする事で、構造の最適化を行わずにサンプリングを行う事も可能です。

.. code-block:: toml

   [sampling.solver]
   type = 'sevennet'
   perturb = 0.0
   base_input_dir = './baseinput_sevennet'
   use_pretrained = false

   [train]
   type = 'sevennet'
   base_input_dir = './sevennet_train_input'
   exe_command = ['', 'sevenn']
   vac_map = []
   restart = false

また、SevenNetのインプットファイル ``input.yaml`` を ``sevennet_train_input/train`` ディレクトリに作成します。
ここではコマンドsevennのインプットファイルを作成しています。
各パラメータの詳しい説明はSevenNetのドキュメントを参照してください。

.. code-block:: yaml

   model:  # model keys should be consistent except for train_* keys
       chemical_species: 'Auto'
       cutoff: 5.0
       channel: 128
       is_parity: False
       lmax: 2
       num_convolution_layer: 5
       irreps_manual:
           - "128x0e"
           - "128x0e+64x1e+32x2e"
           - "128x0e+64x1e+32x2e"
           - "128x0e+64x1e+32x2e"
           - "128x0e+64x1e+32x2e"
           - "128x0e"

       weight_nn_hidden_neurons: [64, 64]
       radial_basis:
           radial_basis_name: 'bessel'
           bessel_basis_num: 8
       cutoff_function:
           cutoff_function_name: 'XPLOR'
           cutoff_on: 4.5
       self_connection_type: 'linear'

       train_shift_scale: False   # customizable (True | False)
       train_denominator: False   # customizable (True | False)

   train:  # Customizable
       random_seed: 1
       is_train_stress: False
       epoch: 5

       optimizer: 'adam'
       optim_param:
           lr: 0.004
       scheduler: 'exponentiallr'
       scheduler_param:
           gamma: 0.99

       force_loss_weight: 0.1
       stress_loss_weight: 1e-06

       per_epoch: 1  # Generate checkpoints every this epoch

       # ['target y', 'metric']
       # Target y: TotalEnergy, Energy, Force, Stress, Stress_GPa, TotalLoss
       # Metric  : RMSE, MAE, or Loss
       error_record:
           - ['Energy', 'RMSE']
           - ['TotalLoss', 'None']

       continue:
           reset_optimizer: True
           reset_scheduler: True
           reset_epoch: True
           checkpoint: 'SevenNet-0_11July2024'

   data:  # Customizable
       batch_size: 4
       data_divide_ratio: 0.1

       # SevenNet automatically matches data format from its filename.
       # For those not `structure_list` or `.pt` files, assumes it is ASE readable
       # In this case, below arguments are directly passed to `ase.io.read`
       data_format_args:
           index: ':'                                # see `https://wiki.fysik.dtu.dk/ase/ase/io/io.html` for more valid arguments

       # validset is needed if you want '_best.pth' during training. If not, both validset and testset is optional.
       load_trainset_path: ['./structure.xyz']  # Example of using ase as data_format, support multiple files and expansion(*)

モデル学習、サンプリングの方法に関してはaenetと同様です。

Maceを利用したサンプリング
----------------------------------------------

Maceのインストール
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``mace`` の利用には、 Maceのインストールが必要です。

下記コマンドにてインストールします。

.. code-block:: bash

    $ python3 -m pip install mace

モデル学習から実行する場合
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Maceでは、モデルを学習してからサンプリングを行う事以外に、
学習済みモデルを利用してサンプリングを行う事も可能です。

学習済みモデルを用いる場合は、 ``[sanmping.solver]`` セクションを下記のように設定します。

.. code-block:: toml

   [sampling.solver]
   type = 'mace'
   perturb = 0.0

サンプリングの方法については、aenetと同様です。

モデル学習から実行する場合
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

モデルの学習から実行する場合は、 ``[train]`` セクションを適切に設定の上で、
``[sanmping.solver]`` セクションに、
``use_pretrained = false`` を追加する必要があります。
``relax = false`` とする事で、構造の最適化を行わずにサンプリングを行う事も可能です。

.. code-block:: toml

   [sampling.solver]
   type = 'mace'
   perturb = 0.0
   base_input_dir = './baseinput_mace'
   use_pretrained = false

   [train]
   type = 'mace'
   base_input_dir = './mace_train_input'
   exe_command = ['', 'mace_run_train']
   vac_map = []
   restart = false

また、Maceのインプットファイル ``input.yaml`` を ``mace_train_input/train`` ディレクトリに作成します。
ここではコマンドmace_run_trainのインプットファイルを作成しています。
各パラメータの詳しい説明はMaceのドキュメントを参照してください。

.. code-block:: yaml

   name: spinel
   foundation_model: "small"
   seed: 2024
   train_file: structure.xyz
   swa: yes
   start_swa: 1200
   max_num_epochs: 5
   device: cpu
   E0s:
     8: -2042.0
     12: -1750.0
     13: -1750.0
   energy_weight: 1.0
   forces_weight: 0.0
   stress_weight: 0.0

モデル学習、サンプリングの方法に関してはaenetと同様です。

CHGNetを利用したサンプリング
----------------------------------------------

CHGNetのインストール
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``chgnet`` の利用には、 CHGNetのインストールが必要です。

下記コマンドにてインストールします。

.. code-block:: bash

    $ python3 -m pip install chgnet

学習済みモデルの利用
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CHGNetでは、モデルを学習してからサンプリングを行う事以外に、
学習済みモデルを利用してサンプリングを行う事も可能です。

学習済みモデルを用いる場合は、 ``[sanmping.solver]`` セクションを下記のように設定します。

.. code-block:: toml

   [sampling.solver]
   type = 'chgnet'
   perturb = 0.0

サンプリングの方法については、aenetと同様です。

モデル学習から実行する場合
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

モデルの学習から実行する場合は、 ``[train]`` セクションを適切に設定の上で、
``[sanmping.solver]`` セクションに、
``use_pretrained = false`` を追加する必要があります。
``relax = false`` とする事で、構造の最適化を行わずにサンプリングを行う事も可能です。

.. code-block:: toml

   [sampling.solver]
   type = 'chgnet'
   perturb = 0.0
   base_input_dir = './baseinput_chgnet'
   use_pretrained = false

   [train]
   type = 'chgnet'
   base_input_dir = './chgnet_train_input'
   exe_command = ['', 'chgnet']
   vac_map = []
   restart = false

また、CHGNetのインプットファイル ``input.yaml`` を ``chgnet_train_input/train`` ディレクトリに作成します。

.. code-block:: yaml

   finetuning : False
   batch_size : 4
   train_ratio : 0.9
   val_ratio : 0.05
   learning_rate : 0.004
   epochs : 100
   model_params:
     atom_fea_dim : 8
     bond_fea_dim : 8
     angle_fea_dim : 8
     num_radial : 9
     num_angular : 9
     num_conv : 2
     atom_conv_hidden_dim : 4
     bond_conv_hidden_dim : 4
     mlp_hidden_dims :
       - 16
       - 16
     atom_graph_cutoff : 7.5
     bond_graph_cutoff : 6.0

このインプットファイルは、CHGNetの学習に必要なパラメータを設定するものであり、
abICS側で定義されているものとなります。
各ファイルのパラメータは下記の通りとなります。

- finetuning : ファインチューニングを行うかどうか。デフォルト値: True
- batch_size : バッチサイズ。デフォルト値: 4
- train_ratio : 学習データの割合。デフォルト値: 0.9
- val_ratio : 検証データの割合。デフォルト値: 0.05
- learning_rate : 学習率。デフォルト値: 0.01
- epochs : エポック数。デフォルト値: 5
- model_params: finetuningがFalseの時に、CHGNetのパラメータとして使用される。パラメータに関しては https://chgnet.lbl.gov/api#class-chgnet を参照の事。

モデル学習、サンプリングの方法に関してはaenetと同様です。

