.. _sec_tutorial:

***************************
チュートリアル
***************************

ここではQuantum ESPRESSO (QE) を用い、実際に能動学習を行い、構造推定を行うまでのチュートリアルを記載します。
なお、本チュートリアルで使用する入力ファイル一式は、 ``examples/standard/active_learning_qe`` にあります。
以下ではgnu parallelおよびaenetはインストールしてあるものとします。
また、計算実行の環境は物性研究所スーパーコンピュータシステムBのohtakaを利用します。

入力ファイルの準備
-----------------------

QE参照ファイルの準備
============================

``baseinput_ref`` にQEのscf計算で参照する入力ファイルをおきます。
以下、サンプルディレクトリにある ``scf.in`` ファイルを記載します。

.. code-block:: toml

    &CONTROL
    calculation = 'relax'
    tstress = .false.
    tprnfor = .false.
    pseudo_dir = '~/qe/pot'
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

- https://www.quantum-espresso.org/upf_files/Al.pbe-nl-kjpaw_psl.1.0.0.UPF
- https://www.quantum-espresso.org/upf_files/Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF
- https://www.quantum-espresso.org/upf_files/O.pbe-n-kjpaw_psl.1.0.0.UPF

このサンプルでは、QE計算時に構造最適化を行うため ``calculation = 'relax'`` を、
計算高速化のため、 ``K_POINTS`` は ``gammma`` を選択しています。


abICSファイルの準備
============================

次に、 :ref:`Input file` を参考にabICS用の入力ファイルを準備します。


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
    srun -n 8 abics_activelearn input_aenet.toml >> active.out
    echo start parallel_run 1
    sh parallel_run.sh

    echo start AL final
    srun -n 8 abics_activelearn input_aenet.toml >> active.out

    #train
    echo start training
    abics_train input_aenet.toml > train.out
    echo Done

最初のSBATCHは物性研スパコンでのジョブスケジューラに関するコマンドです。
ここでは、プロセス数512のMPI並列を実行しています。
ジョブスケジューラに関する詳細は、物性研スパコンのマニュアルを参照してください。

.. code-block:: shell

    # Run reference DFT calc.
    echo start AL sample
    srun -n 8 abics_activelearn input_aenet.toml >> active.out

で、abics_activelearnを用いて、訓練データの大元となる第一原理計算用の入力ファイルを生成します。
初回実行時は、指定した数だけ原子配置をランダムに生成し、
それぞれの原子配置に対して個別のディレクトリを用意した上で、ディレクトリ内に入力ファイルを作成します。
同時に、それらのディレクトリのpathが記載されたファイルrundirs.txtも生成します。
このディレクトリリストを使って、個々の入力に対する第一原理計算ジョブの実行を自動化することができます。
次に得られたファイルをもとに、第一原理計算を実行します。

.. code-block:: shell

    echo start parallel_run 1
    sh parallel_run.sh

``parallel_run.sh`` は、gnu parallelを用いてQEの網羅計算を行うためのスクリプトで、
これによりrundirs.txtに記載されたディレクトリを対象にQEの網羅計算が行われます。
QEの計算結果は、それぞれのディレクトリに格納されます。
QEの網羅計算により、教師データを作成したので、次はaenetでのニューラルネットワークポテンシャルの作成に移ります。
最初に、``abics_activelearn`` を再度実行し、第一原理計算の結果をabics_trainが読み込む共通フォーマットにしたファイルを作成します。

.. code-block:: shell

    echo start AL final
    srun -n 8 abics_activelearn input_aenet.toml >> active.out

次に、学習データをもとにaenetによりニューラルネットワークポテンシャルの作成を行います。
ニューラルネットワークポテンシャルは ``abics_train`` により計算されます。
入力ファイルの[trainer]セクションにあるbase_input_dirに格納された入力ファイルを読み込むことで、計算が実施されます。
計算が無事終了すると、baseinputディレクトリに学習済みのニューラルネットワークが出力されます。

.. code-block:: shell

    #train
    echo start training
    abics_train input_aenet.toml > train.out
    echo Done

以上のプロセスで、能動学習を行うためのAL.shのプロセスが終了となります。

次に、学習したニューラルネットワークポテンシャルを用い、abICSにより最適化構造を求めます。
このプロセスはMC.shで行うことができます。
以下が、MC.shの中身です。

.. code-block:: shell

    #!/bin/sh
    #SBATCH -p i8cpu
    #SBATCH -N 1
    #SBATCH -n 8
    #SBATCH --time=00:30:00

    srun -n 8 abicsAL input_aenet.toml >> aenet.out
    echo Done

abicsALを実行することで ``MCxx`` ディレクトリが作成されます(xxは実行回数)。
active learningを念頭にしており、ALloop.progressを読むことで計算回数などの情報を取得する機能が追加実装されています。
``MCxx`` ディレクトリの下には、レプリカ数分だけのフォルダが作成され、
VASPのPOSCARファイル形式で記載された各ステップごとの原子配置(structure.XXX.vasp)、
最低エネルギーを与えた原子位置(minE.vasp)や、各ステップごとの温度とエネルギー(obs.dat)などが出力されます。
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

.. image:: ../../../image/DOI.pdf
   :width: 800px
   :align: center

なお、DOIについては以下の手順で計算が可能です。

1. MCxxxに移動する。

2. srun -n 8 abicsRXsepT ../input_aenet.toml で Tseparate ディレクトリを作成する
(abicsALを実行した際の並列数に揃える。本チュートリアルでは並列数を8にしているので8に設定)。

3. sampleディレクトリにあるcalc_DOI.py と MgAl2O4.vasp をコピーする。

4. srun -n 8 python3 calc_DOI.py ../input_aenet.toml で温度ごとの反転率を計算する。
(abicsALを実行した際の並列数に揃える。本チュートリアルでは並列数を8にしているので8に設定)。

以上が、チュートリアルになります。