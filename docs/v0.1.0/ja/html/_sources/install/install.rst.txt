.. highlight:: none

ダウンロード
~~~~~~~~~~~~~~~~~~~~~~

abICS のソースコードは `GitHub page <https://github.com/issp-center-dev/abICS>`_ からダウンロードできます。

``$ git clone https://github.com/issp-center-dev/abICS``

必要なライブラリ・環境
~~~~~~~~~~~~~~~~~~~~~~

- python3
- numpy
- scipy
- toml (for parsing input files)
- mpi4py (for parallel tempering)
- pymatgen (for parsing vasp I/O)
- qe-tools (for parsing QE I/O)

VASPをソルバーとして利用する際には、MPI_COMM_SPAWNを利用するためのパッチをあてる必要があります。利用されたい場合には、:doc:`../contact/index` のその他に記載された連絡先までご連絡ください。

ディレクトリ構成
~~~~~~~~~~~~~~~~~~~~~~

abICSのディレクトリ構成は以下のようになっています.
``examples/standard`` には簡易ファイルで実行可能なサンプルが, 
``examples/expert`` にはpythonモジュールを直接用いて作成されたサンプルがあります.
pythonモジュールは ``abics`` ディレクトリ以下に一式格納されています.

:: 

 .
 |-- COPYING
 |-- README.md
 |-- abics/
 |   |-- __init__.py
 |   |-- applications/
 |   |-- mc.py
 |   |-- mc_mpi.py
 |   |-- scripts/
 |   `-- util.py
 |-- docs/
 |   `-- sphinx/
 |-- examples/
 |   |-- expert/
 |   `-- standard/
 |-- make_wheel.sh
 |-- setup.cfg
 `-- setup.py



インストール
~~~~~~~~~~~~~~~~~~~~~~~~~~

1. wheelファイルを作成します.

``$ ./make_wheel.sh``

2. 作成されたファイルを使用して以下のようにインストールします.

``$ pip install dist/abics-*.whl``

インストールディレクトリを変更したい場合には, ``--user`` オプションもしくは ``--prefix = DIRECTORY`` ( ``DIRECTORY`` にインストールしたいディレクトリを指定) オプションを指定してください:

``$ pip install --user dist/abics-*.whl``
