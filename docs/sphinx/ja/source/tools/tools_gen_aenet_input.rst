.. highlight:: none

``abics_gen_aenet_input``
-------------------------------

学習器に aenet を用いる場合、 ``[train]`` セクションの ``base_input_dir`` で指定したディレクトリに
``generate`` 、 ``train`` 、 ``predict`` の 3 つのディレクトリを作り、それぞれに aenet の入力ファイルを
用意する必要があります（内容は :ref:`チュートリアル<sec_tutorial>` の「aenet用の入力ファイル」を参照）。
``abics_gen_aenet_input`` は、abICS 入力ファイルの ``[config]`` セクションと ``[train]`` セクションから
これらの入力ファイルの雛形を生成します。
生成される設定は ``examples/active_learning_qe/aenet_train_input`` と同じもので、
出発点として用い、必要に応じて手で調整してください。

- 記述子: Chebyshev 基底（ ``radial_Rc = 8.0`` , ``radial_N = 16`` , ``angular_Rc = 6.5`` , ``angular_N = 4`` ）, ``RMIN 0.55``
- ネットワーク: 隠れ層 2 層、各 15 ノード（活性化関数 tanh）。ポテンシャルファイル名は ``<元素>.15t-15t.nn``
- 訓練: BFGS 法、500 反復、テストデータ 10%

``abics_gen_aenet_input`` は以下のように使用します::

    $ abics_gen_aenet_input -h
    usage: abics_gen_aenet_input [-h] [-o OUTPUT_DIR] [--force] [tomlfile]

    Generate default aenet training input files (generate/train/predict) from the
    [config] and [train] sections of an abICS input.toml file.

    positional arguments:
      tomlfile              abICS input file (default: input.toml)

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                            Output directory (default: train.base_input_dir from
                            the toml file, or ./aenet_train_input)
      --force               Overwrite existing non-empty output directories

- 対象となる元素は、 ``[config]`` セクションに現れる元素（ ``base_structure`` と ``defect_structure`` の
  すべての元素）から ``[train]`` セクションの ``ignore_species`` を除いたものです。
- 出力先は ``-o`` で指定しない場合、 ``[train]`` セクションの ``base_input_dir``
  （未指定の場合は ``./aenet_train_input`` ）です。
  出力先の ``generate`` 、 ``train`` 、 ``predict`` のいずれかが既に存在して空でない場合はエラーになります。
  ``--force`` を付けると上書きします。
- ``[sampling.solver]`` セクションの ``type`` が ``aenetPyLammps`` の場合は、 ``predict/predict.in`` の代わりに
  LAMMPS 用の ``predict/in.lammps`` を生成します（ :ref:`tutorial_aenet_lammps` 参照）。
  このとき、 ``[sampling.solver]`` の ``ignore_species`` を除いた元素がすべて訓練対象に含まれている必要があり、
  不足している元素があればエラーになります。
- ``[train]`` セクションの ``type`` が ``aenet`` 以外の場合は警告を出します（生成されるのは aenet 用の入力のみです）。

生成されるファイルは以下の通りです（MgAl2O4 スピネルの例。酸素は ``ignore_species`` で除いています）::

    aenet_train_input/
    ├── generate/
    │   ├── generate.in.head
    │   ├── Al.fingerprint.stp
    │   └── Mg.fingerprint.stp
    ├── train/
    │   └── train.in
    └── predict/
        └── predict.in        （aenetPyLammps の場合は in.lammps）

例えば、チュートリアルの入力ファイルから雛形を生成するには次のようにします::

    $ cd examples/active_learning_qe
    $ abics_gen_aenet_input input-aenet.toml -o aenet_train_input_generated
    Generated default aenet training input files for species ['Al', 'Mg'] in aenet_train_input_generated
