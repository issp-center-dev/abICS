.. highlight:: none

[mlref] セクション
-------------------------------

MC計算の結果から原子配置のみを取り出す際のオプションを設定します.
例えばニューラルネットワークモデルの精度評価と訓練データの拡張などに利用します.

以下のようなファイルフォーマットをしています.

  ::
  
        [mlref]
        nreplicas = 3
        ndata = 50


入力形式
^^^^^^^^^^^^
``keyword = value`` の形式でキーワードとその値を指定します.
また, #をつけることでコメントを入力することができます(それ以降の文字は無視されます).

キーワード
^^^^^^^^^^

- レプリカに関する指定

    -  ``nreplicas``

       **形式 :** int型 (自然数)

       **説明 :** レプリカ数を指定します.

    -  ``ndata``

       **形式 :** int型 (自然数)

       **説明 :** 取り出すデータ（原子配位）の数
    
    -  ``sampler``

       **形式 :** 文字列 ("linspace" or "random", デフォルトは "linspace")

       **説明 :** :math:`N` 個生成されている MC サンプルから :math:`N_\text{data}` のデータをどのようにして取り出すか. 

          - "linspace"

            ``numpy.linspace(0, N-1, num=ndata, dtype=int)`` を用いて等間隔に取り出す

          - "random"

            ``numpy.random.choice(range(N), size=ndata, replace=False)`` を用いたランダムサンプリング

    -  ``seed``

       **形式 :** int型 (デフォルトは 0)

       **説明 :** 乱数のシード. 正の値を指定すると, 初回の原子配置の生成, 原子座標の摂動(``perturb``), ``sampler = "random"`` のサンプリングが再現可能になります(レプリカごとに ``seed + 137 * (レプリカ番号)`` が用いられます). 0 の場合はシードを固定しません.

- 実行方式に関する指定

    -  ``exec_mode``

       **形式 :** 文字列 ("staged" or "integrated", デフォルトは "staged")

       **説明 :** ``abics_mlref`` の実行方式.

          - "staged"

            従来の方式. ``abics_mlref`` は第一原理計算ソルバーの入力ファイルを書き出して終了し,
            ユーザーがソルバーを実行したあと, 再度 ``abics_mlref`` を実行して結果を回収します
            (``abics_mlref`` は ``[mlref]`` の ``nreplicas`` 個の MPI プロセスで実行します).

          - "integrated"

            1 回の ``abics_mlref`` 実行で, 配置の生成, ソルバーの実行, 結果の回収までを行います.
            ソルバーは libEnsemble(``pip install libensemble``)を通じて ``abics_mlref`` 自身が
            起動するので, ``abics_mlref`` は **mpiexec / srun を付けずに 1 プロセスで** 実行してください
            (ジョブスクリプト内で実行する場合は, 確保したノード上でこの 1 プロセスがソルバーを
            ``mpirun`` / ``srun`` で起動します). ``[mlref.solver]`` セクションには
            ソルバーの実行コマンド ``path`` と ``run_scheme = "subprocess"`` が必要です.
            ソルバーの実行は ``[mlref.integrated]`` セクションで設定します.
            出力ファイルは staged と同じで, 中断後に再実行すると出力の揃っている配置は再計算されません.

- ``[mlref.integrated]`` サブセクション(``exec_mode = "integrated"`` のとき)

    -  ``nworkers``

       **形式 :** int型 (自然数, デフォルトは 1)

       **説明 :** 同時に実行するソルバーの本数(libEnsemble のワーカー数).

    -  ``nprocs``

       **形式 :** int型 (自然数, デフォルトは 1)

       **説明 :** ソルバー 1 本あたりの MPI プロセス数.

    -  ``nthreads``

       **形式 :** int型 (自然数, デフォルトは 1)

       **説明 :** MPI プロセスあたりのスレッド数(ソルバーのコマンドライン生成に渡されます).

    -  ``timeout``

       **形式 :** float型 (秒, デフォルトは 0)

       **説明 :** ソルバー 1 本の実行時間の上限. 超えると強制終了して失敗扱いにします. 0 は無制限.

    -  ``retries``

       **形式 :** int型 (0 以上, デフォルトは 1)

       **説明 :** ソルバーが失敗(異常終了, タイムアウト, 出力が読めない)したときに再実行する回数.
       それでも失敗した配置があると ``abics_mlref`` は原因を表示して終了コード 1 で終わり,
       ``ALloop.progress`` は更新しません. 原因を取り除いて再実行すると失敗した配置だけが計算されます.

    -  ``scratch_dir``

       **形式 :** 文字列 (デフォルトは空)

       **説明 :** 指定すると, ソルバーの入力をこのディレクトリ(ノードローカルなスクラッチなど)に
       コピーしてそこで実行し, 終了後に出力を ``ALn/RANK/inputXX/baseinputK`` へコピーして
       スクラッチ側を削除します. 空の場合はその場で実行します.

    -  ``copy_back_exclude``

       **形式 :** 文字列のリスト (デフォルトは空)

       **説明 :** ``scratch_dir`` からコピーして戻さないファイルのパターン(例: ``["*.wfc*", "*.save"]``).

    -  ``cores``

       **形式 :** int型 (デフォルトは 0)

       **説明 :** ソルバーの実行に使えるコア数. 0 のときは ``nworkers * nprocs`` とみなします.

  例::

        [mlref]
        nreplicas = 8
        ndata = 50
        exec_mode = "integrated"

        [mlref.integrated]
        nworkers = 16     # 16 本の DFT を同時に実行
        nprocs = 8        # 1 本あたり 8 MPI プロセス
        timeout = 3600
        retries = 1
        scratch_dir = "/scratch/$USER"
        copy_back_exclude = ["*.wfc*"]

        [mlref.solver]
        type = "qe"
        path = "pw.x"
        run_scheme = "subprocess"
        base_input_dir = "./baseinput_ref"
        perturb = 0.05

.. - その他
..
..     -  ``output_frequency``
..
..        **形式 :** list型 (自然数)
..
..        **説明 :**  ``nsteps`` をRXMC計算で出力される配置の数（ ``[replica]`` セクションの ``nsteps/sample_frequency`` の値）のうち何ステップまでを取り出すかのステップ数、 ``sample_frequency`` を配置を抜き出してくる間隔として、 [ ``nsteps`` , ``sample_frequency`` ] のlist形式で配置を抜き出す間隔を指定します(抜き出す最初のステップ数は0に固定しています)。その配置に対応する第一原理計算ソルバーの入力ファイルが各フォルダ内に作成されます。

