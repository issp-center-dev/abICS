.. highlight:: none

abicsRXsepT
-------------------------------

このツールは、RXMC実行の各サンプリングステップで得られる構造とエネルギーを
温度ごとに並べ替えるためのツールです。abICSのRXMC実行が終了した後に使用します::

   $ mpiexec -np NPROCS abicsRXsepT input.toml NSKIP

``NPROCS`` はレプリカの数と同じかそれ以上でなければならず、``input.toml`` はabICSの実行に使用された
abICS入力ファイルに置き換える必要があります。``NSKIP`` はオプションのパラメータで、
各温度でのエネルギー平均を計算する際にスキップするステップ数を指定するために使用します。
結果は ``Tseparate`` ディレクトリに保存され、温度ごとのエネルギー平均は ``Tseparate/energies_T.dat`` に保存されます。

