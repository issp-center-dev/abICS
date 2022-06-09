.. _sec_tutorial:

***************************
aenetを用いた例
***************************

ここではQuantum ESPRESSO (QE) を用い、実際に能動学習を行い、構造推定を行うまでのチュートリアルを記載します。
なお、本チュートリアルで使用する入力ファイル一式は、 ``examples/standard/active_learning_qe`` にあります。
以下ではgnu parallelおよびaenetはインストールしてあるものとします。
また、計算実行の環境は物性研究所スーパーコンピュータシステムBのohtakaを利用します。

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

abICSファイルの準備
============================

abICSによる能動学習の利用にあたっては、abICS、利用する第一原理ソルバー、aenetの
3つに対応する入力ファイルを準備する必要があります。

abICS制御ファイル (``input.toml``)
++++++++++++++++++++++++++++++++++++++++++++++++++++
計算対象とする格子構造の定義と、abICSによる能動学習のループ全体の制御、および
レプリカ交換モンテカルロ法に関するパラメータを設定します。
``st2abics`` ツールを使うことで、結晶構造ファイルから ``input.toml`` のひな形を自動で生成することができます。

::

  $ cd [example_dir]
  $ st2abics st2abics_MgAl2O4.toml MgAl2O4.vasp > input.toml


今回の例題では、このようにして生成した ``input.toml`` の中の
``[sampling.solver]`` セクションのpathをご自身の環境におけるaenetの ``predict.x`` のパスに設定し、
``[train]`` セクションの ``exe_command`` をaenetの ``generate.x`` 、 ``train.x`` 実行時のコマンドに
置き換えます。さらに、 ``[sampling.solver]`` と ``[train]`` で ``ignore_species = ["O"]`` と設定することで動作します。

ここでは、 ``input.toml`` のセクションごとの設定内容をもう少し詳しく解説します。例題をとりあえず
実行したい場合は、読み飛ばしても大丈夫です。

(i)  ``[sampling]`` セクション
****************************************************
.. code-block:: toml

    [sampling]
    nreplicas = 8
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
現状、mpi版の ``predict.x`` はサポートしていないため、 ``nprocs_per_replica`` は
1を指定してください。

(ii)  ``[mlref]`` セクション
****************************************************
.. code-block:: toml

    [mlref]
    nreplicas = 8
    nprocs_per_replica = 1
    ndata = 5

RXMC計算の結果から、ニューラルネットワークモデルの精度評価と訓練データの拡張のために原子配置を取り出す際の
オプションが設定できます。基本的に、 ``nreplicas`` は ``[sampling]`` セクションと同じ値にしてください。
``ndata`` は、RXMC計算で出力される配置の数（ ``[sampling]`` セクションの ``nsteps/sample_frequency`` の値）から、
何個のデータを機械学習ように取り出すかを指定します。従って、RXMC計算で出力される配置の数以下の値に設定してください。

(iii)  ``[sampling.solver]`` セクション
****************************************************
.. code-block:: toml

    [sampling.solver] # RXMC計算に使うソルバーの設定
    type = 'aenet'
    path= 'predict.x-2.0.4-ifort_serial'
    base_input_dir = './baseinput'
    perturb = 0.0
    run_scheme = 'subprocess'
    ignore_species = ["O"]

RXMC計算に使うエネルギーソルバーの設定を行います。今回は、aenetを使ってニューラルネットワークモデルの評価を行います。
``type`` , ``perturb`` , ``run_scheme`` に関しては、能動学習スキームを用いる場合は上の例のまま変更しないでください。
``path`` には、ご自身の環境におけるaenetの ``predict.x`` のパスを指定してください。 ``base_input_dir``
は自由に設定して構いません。
設定したディレクトリの中に ``predict.x`` に対応した入力ファイルが自動で設置されます（後述）。

また、 ``ignore_species`` では、
ニューラルネットワークモデルで「無視」する原子種を指定できます。今回の例題では、Oの副格子は常に占有率1なので、Oの
配置はエネルギーに影響を及ぼしません。こういった場合は、ニューラルネットワークモデルの訓練および評価時に存在を無視した方が、
計算効率が高くなります。

(iv)  ``[mlref.solver]`` セクション
****************************************************
.. code-block:: toml

    [mlref.solver] # 参照第一原理ソルバーの設定
    type = 'qe'
    base_input_dir = ['./baseinput_ref', './baseinput_ref', './baseinput_ref'] #, './baseinput_ref']
    perturb = 0.05
    ignore_species = []

訓練データ（配置エネルギー）の計算に用いるソルバーの設定を行います。この例ではQuantum Espressoで配置エネルギーを求めます。
``base_input_dir`` は自由に設定して構いません。設定したディレクトリの中に、ソルバーの入力ファイルを設置します（後述）。
この例のように、リスト形式で複数設定した場合は、各々の入力を使った計算が順番に実行されます。このときに、2番目以降の計算では
前の計算の最終ステップでの構造が初期座標として用いられます。そして、最後の計算のエネルギーが
学習に使われます。例えば、1つ目の入力ファイルでで精度を犠牲にして高速な構造最適化を行い、2番目以降の入力ファイルで
高精度な設定で構造最適化を行うといった
ことが可能になります。あるいは、格子ベクトルの緩和を行う場合に、設定した平面波カットオフに基づいて計算メッシュをリセット
するために同じ入力の計算を複数回実行するといったことも可能です。

``perturb`` は、ランダムに各原子を変位させることで、対称性を崩した構造から構造最適化を開始するための
設定です。この場合は、構造緩和を行う原子を全て0.05 Å、ランダムな方向に変位させた構造から1番目の計算が開始されます。

また、 ``ignore_species`` は、第一原理ソルバーを訓練データ生成に用いる場合は空リストを指定しますが、一部の元素を無視するような
モデルを使って訓練データを生成する場合は、無視する元素を指定します。



(v)  ``[train]`` セクション
****************************************************
.. code-block:: toml

    [train] # モデル学習器の設定
    type = 'aenet'
    base_input_dir = './aenet_train_input'
    exe_command = ['generate.x-2.0.4-ifort_serial',
                  'srun train.x-2.0.4-ifort_intelmpi']
    ignore_species = ["O"]
    vac_map = []
    restart = false

訓練データから配置エネルギー予測モデルを学習する学習器の設定を行います。現在のところ、abICSではaenetのみに
対応しています。 ``base_input_dir`` は自由に設定して構いません。設定したディレクトリの中に、学習器の設定ファイルを
設置します（後述）。 ``exe_command`` にはaenetの ``generate.x`` と ``train.x`` へのパスを指定します。
``train.x`` についてはMPI並列版が利用可能で、その場合は、上の例で示すように、MPI実行するためのコマンド
（ ``srun`` 、 ``mpirun`` など）を合わせて設定してください。

また、 ``ignore_species`` は、第一原理ソルバーを訓練データ生成に用いる場合は空リストを指定しますが、
一部の元素を無視するような
モデルを使って訓練データを生成する場合は、無視する元素を指定します。 ``vac_map`` 、 ``restart`` については現状対応していないので、
例のように設定してください。

(vi)  ``[config]`` セクション
****************************************************
.. code-block:: toml

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

モンテカルロサンプリングを行う原子配置の情報を設定します。基本的に ``st2abics`` ツールで生成されたものを
そのまま利用できます。

QE参照ファイルの準備
============================

``baseinput_ref`` にQEのscf計算で参照する入力ファイルをおきます。
以下、サンプルディレクトリにある ``scf.in`` ファイルを記載します。

.. code-block::

    &CONTROL
    calculation = 'relax'
    tstress = .false.
    tprnfor = .false.
    pseudo_dir = './pseudo'
    disk_io = 'low'
    wf_collect = .false.
    /
    &SYSTEM
      ecutwfc      =  60.0
      occupations  = "smearing"
      smearing     = "gauss"
      degauss      = 0.01
    /
    &electrons
      mixing_beta = 0.7
      conv_thr = 1.0d-8
      electron_maxstep = 100
    /
    &ions
    /
    ATOMIC_SPECIES
    Al 26.981 Al.pbe-nl-kjpaw_psl.1.0.0.UPF
    Mg 24.305 Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF
    O  16.000 O.pbe-n-kjpaw_psl.1.0.0.UPF
    ATOMIC_POSITIONS crystal

    K_POINTS gamma

なお、擬ポテンシャルを格納したディレクトリ ``pseudo_dir`` や
``ATOMIC_SPECIES`` で使用する擬ポテンシャルについて、自分の環境に従い書き換える必要があります。
なお、本サンプルで使用している擬ポテンシャルは以下のリンクからダウンロードできます。

- https://pseudopotentials.quantum-espresso.org/upf_files/Al.pbe-nl-kjpaw_psl.1.0.0.UPF
- https://pseudopotentials.quantum-espresso.org/upf_files/Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF
- https://pseudopotentials.quantum-espresso.org/upf_files/O.pbe-n-kjpaw_psl.1.0.0.UPF

このサンプルでは、QE計算時に構造最適化を行うため ``calculation = 'relax'`` を、
計算高速化のため、 ``K_POINTS`` は ``gammma`` を選択しています。

aenetを使った訓練および配置エネルギ－計算用の入力ファイル
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

aenet用の入力ファイルを ``[train]`` セクションの ``base_input_dir`` で
設定したディレクトリ内の ``generate`` 、 ``train`` 、および ``predict``
ディレクトリに設置します。

generate
********

aenetでは、訓練用の原子配置とエネルギーのデータを、原子環境記述子とエネルギーの
関係に変換した中間バイナリフォーマットにまとめてから訓練を行います。この変換を
行う ``generate.x`` 用の
入力ファイルを ``generate`` ディレクトリに設置します。

まず、元素種ごとの
記述子設定ファイルを用意します。ファイル名は任意ですが、チュートリアルでは
``Al.fingerprint.stp`` , ``Mg.fingerprint.stp`` のような名前にしています。
例として ``Al.fingerprint.stp`` の内容を示します：

.. code-block ::

  DESCR
   N. Artrith and A. Urban, Comput. Mater. Sci. 114 (2016) 135-150.
   N. Artrith, A. Urban, and G. Ceder, Phys. Rev. B 96 (2017) 014112.
  END DESCR

  ATOM Al # 元素を指定

  ENV 2 # ATOMで指定した元素と相互作用する元素種の数と元素名を指定
  Al
  Mg

  RMIN 0.55d0 # 原子間の最隣接距離

  BASIS type=Chebyshev # チェビシェフ記述子の設定
  radial_Rc = 8.0  radial_N = 16 angular_Rc = 6.5  angular_N = 4

記述子設定の詳細についてはaenetのドキュメントをご参照ください。

次に、
``generate.in.head`` という名前の以下のようなファイルを準備します：

.. code-block ::

    OUTPUT aenet.train

    TYPES
    2
    Al -0.0  ! eV
    Mg -0.0  ! eV

    SETUPS
    Al   Al.fingerprint.stp
    Mg    Mg.fingerprint.stp


``OUTPUT`` には必ず ``aenet.train`` を指定してください。
``TYPES`` 以下には訓練データ中の元素種とその数を指定します。
元素種ごとにエネルギーの基準を指定することもできますが、基本的には
0に設定しておくのが無難です。
``SETUPS`` 以下には元素種ごとの記述子設定ファイルを指定します。
ファイルの末尾には必ず改行が入っていることを確認してください。
abICSは ``generate.in.head`` の末尾に座標ファイルのリストを
追加して、 ``generate.in`` を生成し、 ``generate.x`` を実行します。

train
*****

``generate`` で生成された訓練データを読み込み、訓練を行う
``train.x`` 用の入力ファイルを ``train`` ディレクトリに設置します。
ファイル名は ``train.in`` としてください：

.. code-block ::

    TRAININGSET aenet.train
    TESTPERCENT 10
    ITERATIONS  500

    MAXENERGY 10000

    TIMING

    !SAVE_ENERGIES

    METHOD
    bfgs

    NETWORKS
    ! atom   network         hidden
    ! types  file-name       layers  nodes:activation
      Al     Al.15t-15t.nn    2      15:tanh 15:tanh
      Mg       Mg.15t-15t.nn    2      15:tanh 15:tanh

基本的には、 ``NETWORKS`` セクション以外は変更の必要はありません。
``NETWORKS`` セクションでは、生成する元素種ごとのポテンシャル
ファイル名と、ニューラルネットワーク構造、および活性化関数を指定します。

predict
*******

訓練したポテンシャルモデルを使って入力座標に対してエネルギーを
評価するための ``predict.x`` 用の入力ファイル ``predict.in`` を、 ``predict``
ディレクトリに設置します：

.. code-block ::

    TYPES
    2
    Mg
    Al

    NETWORKS
    Mg  Mg.15t-15t.nn
    Al  Al.15t-15t.nn

    VERBOSITY low

``TYPES`` セクションには元素種の数と元素名を、 ``NETWORKS``
セクションには元素種ごとのポテンシャルファイル名（ ``train.in`` で
設定したもの）を入力してください。

また、 ``VERBOSITY`` は必ず ``low`` に設定してください。

計算実行
-----------------------

サンプルスクリプトには、計算手順を簡略化するためのスクリプト  ``AL.sh`` と ``MC.sh`` が準備されています。
これらを交互に実行することで能動学習およびモンテカルロによる構造推定が実施されれます。

それでは、最初に ``AL.sh`` の中身をみてます。
なお、これらのシェルスクリプトの実行前に、 ``run_pw.sh`` を ``chmod u+x`` で権限を変更する必要があります。
``run_pw.sh`` はQEの計算を実行するためのスクリプトで、後述する ``parallel_run.sh`` 内部で呼び出されます。

.. code-block:: shell

    #!/bin/sh
    #SBATCH -p i8cpu
    #SBATCH -N 4
    #SBATCH -n 512
    #SBATCH -J spinel
    #SBATCH -c 1
    #SBATCH --time=0:30:00

    # Run reference DFT calc.
    echo start AL sample
    srun -n 8 abics_mlref input.toml >> active.out

    echo start parallel_run 1
    sh parallel_run.sh

    echo start AL final
    srun -n 8 abics_mlref input.toml >> active.out

    #train
    echo start training
    abics_train input.toml > train.out

    echo Done

最初の ``#SBATCH`` で始まる数行は物性研スパコンでのジョブスケジューラに関するコマンドです。
ここでは、プロセス数512のMPI並列を実行しています。
また、 ``srun`` は並列環境でプログラムを実行するためのコマンドです（ ``mpiexec`` に相当します）。
ジョブスケジューラに関する詳細は、実際に利用する計算機のマニュアルを参照してください。

.. code-block:: shell

    # Run reference DFT calc.
    echo start AL sample
    srun -n 8 abics_mlref input.toml >> active.out

で、 ``abics_mlref`` を用いて、訓練データの大元となる第一原理計算用の入力ファイルを生成します。
初回実行時は、指定した数だけ原子配置をランダムに生成し、
それぞれの原子配置に対して個別のディレクトリを用意した上で、ディレクトリ内に入力ファイルを作成します。
同時に、それらのディレクトリのpathが記載されたファイル ``rundirs.txt`` も生成します。
このディレクトリリストを使って、個々の入力に対する第一原理計算ジョブの実行を自動化することができます。
次に得られたファイルをもとに、第一原理計算を実行します。

.. code-block:: shell

    echo start parallel_run 1
    sh parallel_run.sh

``parallel_run.sh`` は、gnu parallelを用いてQEの網羅計算を行うためのスクリプトで、
これによりrundirs.txtに記載されたディレクトリを対象にQEの網羅計算が行われます。
QEの計算結果は、それぞれのディレクトリに格納されます。
QEの網羅計算により、教師データを作成したので、次はaenetでのニューラルネットワークポテンシャルの作成に移ります。
最初に、 ``abics_mlref`` を再度実行し、第一原理計算の結果をabics_trainが読み込む共通フォーマットにしたファイルを作成します。

.. code-block:: shell

    echo start AL final
    srun -n 8 abics_mlref input.toml >> active.out

次に、学習データをもとにaenetによりニューラルネットワークポテンシャルの作成を行います。
ニューラルネットワークポテンシャルは ``abics_train`` により計算されます。
入力ファイルの ``[train]`` セクションにある ``base_input_dir`` に格納された入力ファイルを読み込むことで、計算が実施されます。
計算が無事終了すると、 ``baseinput`` ディレクトリに学習済みのニューラルネットワークが出力されます。

.. code-block:: shell

    #train
    echo start training
    abics_train input.toml > train.out

以上のプロセスで、能動学習を行うための ``AL.sh`` のプロセスが終了となります。

次に、学習したニューラルネットワークポテンシャルを用い、abICSにより最適化構造を求めます。
このプロセスは ``MC.sh`` で行うことができます。
以下が、 ``MC.sh`` の中身です。

.. code-block:: shell

    #!/bin/sh
    #SBATCH -p i8cpu
    #SBATCH -N 1
    #SBATCH -n 8
    #SBATCH --time=00:30:00

    srun -n 8 abics_sampling input.toml >> aenet.out

    echo Done

``abics_sampling`` を実行することで ``MCxx`` ディレクトリが作成されます(xxは実行回数)。
``active learning`` を念頭にしており、ALloop.progressを読むことで計算回数などの情報を取得する機能が追加実装されています。
``MCxx`` ディレクトリの下には、レプリカ数分だけのフォルダが作成され、
VASPのPOSCARファイル形式で記載された各ステップごとの原子配置(``structure.XXX.vasp``)、
最低エネルギーを与えた原子位置(``minE.vasp``)や、各ステップごとの温度とエネルギー(``obs.dat``)などが出力されます。
詳細については `abICSマニュアルの出力ファイル <https://issp-center-dev.github.io/abICS/docs/sphinx/ja/build/html/outputfiles/index.html>`_ を参考にしてください。

上の手続きで得られた結果は、aenetにより求められたニューラルネットワークポテンシャルの精度に依存します。
はじめのステップではランダムな配置をもとに学習を行ったので、低温の構造については精度が低いことが予想されます。
そこで、モンテカルロで推定された構造に対して、
再度第一原理計算でエネルギーを計算し再学習させるステップを繰り返すことで、
全温度領域での精度を高めることが期待されます。
このプロセスは、AL.shとMC.shを順番に繰り返すことで計算できます。
実際に下図に反転率(DOI)を計算した結果を掲載します。
この例では最初の一回目の結果がMC0、その後MC1, MC2, ..., MC5と5回実行させています。
最初の一回目が、他のものとかなりずれていることから精度が出ていないことが予想されます。
一方で、一度モンテカルロを行った結果を元に学習させると、その次からはほぼ同じような値が得られていることがわかります。

.. image:: ../../../image/doi_aenet.*
   :width: 800px
   :align: center

なお、DOIについては以下の手順で計算が可能です。

1. MCxxxに移動する。

2. ``srun -n 8 abicsRXsepT ../input.toml`` で ``Tseparate`` ディレクトリを作成する
(abics_samplingを実行した際の並列数に揃える。本チュートリアルでは並列数を8にしているので8に設定)。

3. sampleディレクトリにある ``calc_DOI.py`` と ``MgAl2O4.vasp`` をコピーする。

4. ``srun -n 8 python3 calc_DOI.py ../input.toml`` で温度ごとの反転率を計算する。
(abics_samplingを実行した際の並列数に揃える。本チュートリアルでは並列数を8にしているので8に設定)。

以上が、チュートリアルになります。
