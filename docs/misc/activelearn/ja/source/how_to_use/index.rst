.. _sec_basic_usage:

***************************
基本的な使用方法
***************************

.. highlight:: none

aenetのインストール
-----------------------

abICSでは、ニューラルネットワークモデルの構築のためにaenetを利用します。
aenetは http://ann.atomistic.net からダウンロードすることが可能です。
DocumentationのInstallationに従い、インストールを実施します。
なお、abICSでは、ニューラルネットワークの学習と評価にaenetの ``train.x`` と ``predict.x`` を使います。
``train.x`` についてはMPI並列版が利用可能ですが、 ``predict.x`` についてはMPIを使用しない実行ファイル(serial)を使用する必要があります。
そのため、makefilesの下にあるserial版もインストールするようにしてください。

GNU parallelのインストール(オプション)
-----------------------------------------
チュートリアルではGNU parallelを用いて、Quantum Espressoによる第一原理計算を並列実行します。
そこで、GNU parallelを最初にインストールします。
GNU parallelは https://www.gnu.org/software/parallel/ からダウンロードすることができます(Macの場合はhomebrewにより直接インストールすることも可能です)。
また、オンラインで無償で利用できるチュートリアルを本サイトで配布しています。
インストールは基本的には、ダウンロードして解凍したディレクトリ に移動した後、

::

  $ ./configure && make && make install

でインストールできます。詳細な設定は公式マニュアルを参考にしてください。


入力ファイルの準備
-----------------------

abICSの利用にあたっては、4通りの入力ファイルを準備する必要があります。

abICS制御ファイル (input.toml)
++++++++++++++++++++++++++++++++++++++++++++++++++++
計算対象とする格子構造の定義と、abICSによる能動学習のループ全体の制御、および
レプリカ交換モンテカルロ法に関するパラメータを設定します。
st2abicsツールを使うことで、結晶構造ファイルからinput.tomlのひな形を自動で生成することができます。

::

  $ cd [example_dir]
  $ st2abics st2abics_MgAl2O4.toml MgAl2O4.vasp > input.toml


今回の例題では、このようにして生成したinput.tomlの中の
``[solver]`` セクションのpathをご自身の環境におけるaenetの ``predict.x`` のパスに設定し、
``[trainer]`` セクションのexe_commandをaenetの ``generate.x`` 、 ``train.x`` 実行時のコマンドに
置き換えます。さらに、 ``[solver]`` と ``[trainer]`` で ``ignore_species = ["O"]`` と設定することで動作します。

ここでは、input.tomlのセクションごとの設定内容をもう少し詳しく解説します。例題をとりあえず
実行したい場合は、読み飛ばしても大丈夫です。

(i)  [replica]セクション
****************************************************
.. code-block:: toml

    [replica] 
    nreplicas = 15            
    nprocs_per_replica = 1    
    kTstart = 600.0           
    kTend = 2000.0            
    nsteps = 6400 
    RXtrial_frequency = 4
    sample_frequency = 16
    print_frequency = 1
    reload = false

レプリカ交換モンテカルロ(RXMC)法のレプリカの数や温度範囲などに関する設定を行います（マニュアル参照リンク）。
今回は、RXMC計算のエネルギーソルバーとしてaenetの ``predict.x`` を用います。
現状、mpi版の ``predict.x`` はサポートしていないため、nprocs_per_replicaは
1を指定してください。

(ii)  [replicaRef]セクション
****************************************************
.. code-block:: toml

    [replicaRef] 
    nreplicas = 15
    nprocs_per_replica = 1
    nsteps = 400
    sample_frequency = 20

RXMC計算の結果から、ニューラルネットワークモデルの精度評価と訓練データの拡張のために原子配置を取り出す際の
オプションが設定できます。基本的に、nreplicasやnprocs_per_replicaは[replica]セクションと同じ値にしてください。
nstepsは、RXMC計算で出力される配置の数（[replica]セクションのnsteps/sample_frequencyの値）のうち、
最初の何ステップまでを取り出すかを指定します。従って、RXMC計算で出力される配置の数以下の値に設定してください。
より小さい値に設定することで、RXMC計算の初期の部分を重点的に取ってくることができます（RXMC計算が平衡化しきる前
を重点的に抜き出したい場合など）。
また、[replicaRef]セクションのsample_frequencyは、配置を抜き出してくる間隔を指定します。
上記の場合は、20ステップ間隔でstep 0, 20, 40, ... 380の20通りの構造が1つのレプリカごとに抜き出されます。
全部で20×15=300通りの配置に対する第一原理計算の入力ファイルが生成されることになります。

(iii)  [solver]セクション
****************************************************
.. code-block:: toml

    [solver] # RXMC計算に使うソルバーの設定
    type = 'aenet'
    path= '~/git/aenet/bin/predict.x-2.0.4-ifort_serial'
    base_input_dir = './baseinput'
    perturb = 0.0
    run_scheme = 'subprocess' 
    ignore_species = ["O"]

RXMC計算に使うエネルギーソルバーの設定を行います。今回は、aenetを使ってニューラルネットワークモデルの評価を行います。
type, perturb, run_schemeに関しては、能動学習スキームを用いる場合は上の例のまま変更しないでください。
pathには、ご自身の環境におけるaenetの ``predict.x`` のパスを指定してください。base_input_dirは自由に設定して構いません。
設定したディレクトリの中に ``predict.x`` に対応した入力ファイルを設置します（後述）。

また、ignore_speciesでは、
ニューラルネットワークモデルで「無視」する原子種を指定できます。今回の例題では、Oの副格子は常に占有率1なので、Oの
配置はエネルギーに影響を及ぼしません。こういった場合は、ニューラルネットワークモデルの訓練および評価時に存在を無視した方が、
計算効率が高くなります。

.. code-block:: toml

    [solverRef] # 参照第一原理ソルバーの設定
    type = 'qe'
    path = '' # active learning では無視される
    base_input_dir = ['./baseinput_ref', './baseinput_ref', './baseinput_ref'] #, './baseinput_ref']
    perturb = 0.05
    run_scheme = 'subprocess'
    only_input = true
    ignore_species = []
    vac_convert = []

    [trainer] # モデル学習器の設定
    type = 'aenet'
    base_input_dir = './aenet_train_input'
    exe_command = ['~/git/aenet/bin/generate.x-2.0.4-ifort_serial', 
                  'srun ~/git/aenet/bin/train.x-2.0.4-ifort_intelmpi']
    ignore_species = ["O"]
    vac_map = []
    restart = false

    [config] # 以下、結晶格子の情報と、格子上に配置される原子や空孔の情報が続く
    unitcell = [[8.1135997772, 0.0000000000, 0.0000000000],
                [0.0000000000, 8.1135997772, 0.0000000000],
                [0.0000000000, 0.0000000000, 8.1135997772]]
    supercell = [1,1,1]

    [[config.base_structure]]
    type = "O"
    coords = [
        [0.237399980, 0.237399980, 0.237399980],
        [0.762599945, 0.762599945, 0.762599945],
        [0.512599945, 0.012600004, 0.737399936],
        [0.487399966, 0.987399936, 0.262599975],
        ... 
    


第一原理ソルバーの入力ファイル
++++++++++++++++++++++++++++++++++++++++++++++++++++

aenetを使った訓練用の入力ファイル
++++++++++++++++++++++++++++++++++++++++++++++++++++

aenetを使った配置エネルギー計算用の入力ファイル
++++++++++++++++++++++++++++++++++++++++++++++++++++


Active learningの実施
-----------------------



訓練データの生成
++++++++++++++++++++++++++++++++++++++++++++++++++++


(i)  第一原理計算用入力ファイルの生成
****************************************************

abics_activelearnを用いて、訓練データの大元となる第一原理計算用の入力ファイルを生成します。初回実行時は、
指定した数だけ原子配置をランダムに生成し、それぞれの原子配置に対して個別のディレクトリを用意し、入力ファイルを
設置します。同時に、それらのディレクトリのpathが記載されたファイルrundirs.txtも生成します。
このディレクトリリストを使って、個々の入力に対する第一原理計算ジョブの実行を自動化することができます。
本チュートリアルでは、スケジューラとして
slurmがインストールされている共用計算機を念頭に、gnu parallelを利用した一括実行方法を紹介します。

abics_activelearnの入力ファイルの情報は以下の通りで、[solverRef]セクションにある情報を読み取り、第一原理計算用の入力ファイルを生成します。

- type : 第一原理計算ソルバーを表します。 'vasp', 'qe', 'openmx'が選択できます。

- path: aenetの実行ファイル(predict.x)へのパスを指定します。

- base_input_dir: 第一原理ソルバーの参照する入力ファイルが格納されたディレクトリのリストを表します。

- perturb, run_scheme, ignore\_species:  [solver]セクションと同様。

abics\_activelearnでは、baseinput\_dirにあるフォルダの入力ファイルを利用し、計算用入力ファイルの生成が行われます
(厳しい条件で計算させたい場合などに、複数入力ファイルを準備するなどの応用も可能)。
実行すると中間ファイルを出力し、abics\_activelearnを実行した回数を記録します。
それを読み取ることで、baseinput\_dirの対応する入力フォルダに格納された入力ファイルを読み込みます。
実行回数がbaseinput\_dirの要素よりも多い場合には、エラーを吐きます。

(ii)  第一原理計算の実行
****************************************************

(i)で作成した入力ファイルをもとに第一原理計算を実行します。
gnu parallelを用いて計算する場合には、rundirs.txtを-aオプションで指定することで簡単に並列計算を実施することができます。

baseinput\_dirで複数回実行する場合には、abics\_activelearnの実行と(ii)の計算を回数分実行する必要があります。


aenetを用いたニューラルネットワークポテンシャルの作成
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

abics_trainで作成した第一原理計算の計算結果をaenetを用いて学習し、ニューラルネットワークポテンシャルを作成します。
abics_trainの入力情報は[trainer]セクションで記載します。各パラメータの説明は以下の通りです。

- type: ニューラルネットワークポテンシャルを生成するための学習器 (現状では 'aenet' のみ)
- base_input_dir:  学習器が参照する入力ファイルが格納されたディレクトリのパス。
- exe_command:  実行コマンドのリスト。aenetを利用する場合は、generate.xとtrain.xへのパスを指定する必要があります。
- ignore_species: [solver]セクションと同様。

aenetをソルバーとして利用しモンテカルロ法を利用した構造推定
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

abICSを走らせてモンテカルロ法を利用してします。
[solver]セクションのtypeをaenetに変更し、pathに ``predict.x`` へのパスを通します。
また、base_input_dirには、 ``predict.x`` の実行用に、(b)で作成されたニューラルネットワークと入力ファイルを置かれたディレクトリへのパスを指定します。
入力ファイル作成後、abICSを実行させることで、構造推定が行われます。

