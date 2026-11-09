.. highlight:: none

[al] セクション
-------------------------------

能動学習ループ全体を 1 つのプログラムで実行する ``abics_al`` の設定です.
``abics_al`` は, ``[mlref] exec_mode = "integrated"`` の ``abics_mlref``(第一原理計算),
``abics_train``(訓練), ``abics_sampling``(モンテカルロ)を ``AL0 → train → MC0 → AL1 → ...`` の順に,
``max_iterations`` 回のモンテカルロ計算が揃うまで繰り返します.
第一原理計算はこのプロセス自身が libEnsemble を通じて実行し, 訓練とサンプリングは
``train_command`` と ``sampling_command`` で指定したコマンドを子プロセスとして起動します.
``abics_mlref`` の integrated モードと同様に, ``abics_al`` は **mpiexec / srun を付けずに 1 プロセスで**
実行してください. libEnsemble が必要です(``pip install libensemble``).

``[mlref]`` セクションは integrated モードの設定(``[mlref.solver]`` の ``path`` と ``run_scheme = "subprocess"``,
``[mlref.integrated]``)が必要です(``exec_mode`` は ``abics_al`` が integrated に上書きします).

ループの状態は ``ALloop.progress`` と ``MC*`` ディレクトリから読み取るので, 中断後に ``abics_al`` を
再実行すると完了済みのステージは繰り返さず続きから進みます. サンプリングの途中で止まった場合
(``MC<n>`` はあるが ``ALloop.progress`` に記録されていない)は, その ``MC<n>`` を ``al_driver/`` 以下に
退避してサンプリングをやり直します. ``max_iterations`` を増やして再実行すれば, 完了したループを延長できます.

以下のようなファイルフォーマットをしています.

  ::

        [al]
        max_iterations = 5
        sampling_command = "srun -n {nprocs} abics_sampling {input}"
        walltime = 86400
        walltime_margin = 1800
        resubmit_command = "sbatch job.sh"

        [al.stage_estimates]
        mlref = 7200
        train = 1800
        sampling = 3600

キーワード
^^^^^^^^^^

-  ``max_iterations``

   **形式 :** int型 (自然数, デフォルトは 1)

   **説明 :** 完了させるモンテカルロ計算の数(``MC0`` から ``MC(max_iterations-1)`` まで).
   コマンドラインの ``--max-iterations N`` でも指定できます.

-  ``train_command``

   **形式 :** 文字列 (デフォルトは ``"abics_train {input}"``)

   **説明 :** 訓練ステージのコマンド. ``{input}`` は入力ファイル名に置き換えられます
   (``{nreplicas}`` と ``{nprocs}`` も使えます).
   GPU での訓練など資源の指定が必要な場合は ``"srun -n 1 --gpus 1 abics_train {input}"`` のように
   ランチャを含めて指定します.

-  ``sampling_command``

   **形式 :** 文字列 (デフォルトは ``"mpiexec -n {nprocs} abics_sampling {input}"``)

   **説明 :** サンプリングステージのコマンド. ``{nprocs}`` は ``abics_sampling`` に必要な MPI プロセス数
   (``[sampling]`` の ``nreplicas`` × アンサンブルのメンバー数, 通常は ``nreplicas``)に,
   ``{nreplicas}`` は ``nreplicas`` に置き換えられます. Slurm では
   ``"srun -n {nprocs} abics_sampling {input}"`` のようにします.

   コマンドはシェルを介さずに空白で分割して実行されます. 置き換えられるのは上記のプレースホルダだけで,
   それ以外の波括弧(``--cpu-bind=map_cpu:{0,1}`` など)はそのまま残ります.
   各引数の ``~`` と環境変数は展開されます.

-  ``stage_timeout``

   **形式 :** float型 (秒, デフォルトは 0)

   **説明 :** 訓練・サンプリングの子プロセスを終了させるまでの時間. 0 は無制限.
   超過すると SIGTERM を送り(``mpirun`` / ``srun`` が MPI プロセスを片付けられるように), 30 秒後にまだ
   終わらなければ SIGKILL します. ステージは失敗(終了コード 1)扱いになります.

-  ``walltime``

   **形式 :** float型 (秒, デフォルトは 0)

   **説明 :** この ``abics_al`` の実行に使える壁時間. 0 は無制限.
   各ステージの開始前に, そのステージの所要時間の見積もりと ``walltime_margin`` の和が残り時間に
   収まらないと判断すると, そのステージを始めずに終了コード 2 で終了します(完了済みのステージは保存されています).
   見積もりには, これまでの実行でそのステージにかかった時間の最大値(``al_driver/durations.json`` に記録され,
   次回の実行に引き継がれます)を使い, まだ一度も実行していないステージでは ``stage_estimates`` の値を使います.
   どちらもない場合は警告を出して開始します.

-  ``walltime_margin``

   **形式 :** float型 (秒, デフォルトは 600)

   **説明 :** 壁時間の余裕. ステージ終了後の後処理やジョブの終了処理のために確保します.

-  ``stage_estimates``

   **形式 :** テーブル (``mlref``, ``train``, ``sampling`` をキーとする秒数, デフォルトは空)

   **説明 :** 各ステージの所要時間の初期見積もり. 最初のジョブで壁時間の判断に使われます
   (例: ``[al.stage_estimates]`` に ``sampling = 3600``).

-  ``progress_wait``

   **形式 :** float型 (秒, デフォルトは 30)

   **説明 :** 訓練・サンプリングの子プロセスが終了したあと, 書き込まれた ``ALloop.progress`` が
   見えるようになるまで待つ最大時間(共有ファイルシステムの遅延対策).

-  ``resubmit_command``

   **形式 :** 文字列 (デフォルトは空)

   **説明 :** 壁時間のために終了するときに実行するコマンド(例: ``"sbatch job.sh"``).
   ジョブスクリプトから ``abics_al`` を呼び, ここに同じジョブスクリプトの投入コマンドを書いておくと,
   ループが自動的に継続します.

-  ``generate_train_input``

   **形式 :** bool型 (デフォルトは true)

   **説明 :** ``[train] type = "aenet"`` で ``base_input_dir`` が存在しない場合に,
   最初の訓練の前に ``abics_gen_aenet_input`` で訓練入力の雛形を生成します.

出力
^^^^^^^^^^

- ``al_summary.dat``: 各ステップ(``AL1`` 以降)の全レプリカをまとめた誤差統計(配置数, RMSE, MAE, 最大誤差).
  ``AL<n>/convergence.dat`` から作られ, ループを続けるかの判断に使います.
- ``al_driver/``: 訓練・サンプリングの子プロセスの出力(``train_<n>.log``, ``sampling_<n>.log``. 再実行時は追記),
  ステージの所要時間(``durations.json``), 途中で止まったサンプリングの退避先(``MC<n>.interrupted-<日時>``),
  および目標回数に達したときの最終モデルの記述 ``model_metadata.json`` (abICS のバージョン,
  入力ファイルのハッシュ, 完了した反復回数と参照計算の数, 参照ソルバーと学習器の設定,
  元素種・温度範囲・スーパーセル, 各ステップの誤差統計).

終了コードは, 目標回数に達したとき 0, 壁時間のために途中で止まったとき 2, ステージが失敗したとき 1 です.
