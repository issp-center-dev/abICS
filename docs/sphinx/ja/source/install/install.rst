.. highlight:: none

必要なライブラリ・環境
~~~~~~~~~~~~~~~~~~~~~~

abICS をインストール・実行するには、 バージョン3.7 以上の Python が必要です。
また、以下の Python パッケージが必要です。

- numpy
- scipy
- toml
- mpi4py
- pymatgen (>=2019.12.3)
- qe-tools

これらのライブラリは自動でインストールされますが、 mpi4py と pymatgen はあらかじめ関連ソフトウェアが必要です。

- mpi4py をインストールするには、なんらかのMPI 環境をあらかじめインストールしておく必要があります。
- pymatgen をインストールするには、 Cython をインストールしておく必要があります。::

   $ pip3 install cython

VASPをソルバーとして利用する際には、MPI_COMM_SPAWNを利用するためのパッチをあてる必要があります。利用されたい場合には、:doc:`../contact/index` のその他に記載された連絡先までご連絡ください。


PyPI からインストールする
~~~~~~~~~~~~~~~~~~~~~~~~~~

abICS は PyPI に登録されているため、 ``pip`` コマンドで簡単にインストールできます。::

   $ pip3 install abics

書き込み権限がないなどで、ユーザローカルのディレクトリにインストールする場合には ``--user`` オプションを追加してください。
この場合、 ``~/.local/`` 以下に実行可能スクリプトやライブラリがインストールされます。
また、インストールディレクトリを指定したい場合には、 ``--prefix=DIRECTORY`` ( ``DIRECTORY`` はインストールしたいディレクトリ) オプションを指定してください:

   ``$ pip3 install --user abics``


ソースからインストールする
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

多くの場合には PyPI からインストールすれば良いですが、機能追加する場合などはソースからインストールしてください。

ダウンロード
..................

abICS のソースコードは `GitHub page <https://github.com/issp-center-dev/abICS>`_ からダウンロードできます。

``$ git clone https://github.com/issp-center-dev/abICS``


ディレクトリ構成
.......................

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
.................

- ``pip3 install`` の引数に abICS のルートディレクトリを渡すことでインストール可能です ::

   $ pip3 install ./abICS


アンインストール
~~~~~~~~~~~~~~~~~

- ``pip3 uninstall abics`` でアンインストールできます.
