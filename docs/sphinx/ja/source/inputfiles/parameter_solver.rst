.. highlight:: none

[solver] セクション
-------------------------------

ソルバーの種類 (VASP, QE, ...)、ソルバーへのパス、不変な入力ファイルのあるディレクトリなど（第一原理計算）ソルバーのパラメータを指定します.
以下のようなファイルフォーマットをしています.

  :: 
  
    [solver]
    type = 'vasp'
    path = './vasp'
    base_input_dir = './baseinput'
    perturb = 0.1
    run_scheme = 'mpi_spawn_ready'

入力形式
^^^^^^^^^^^^
``keyword = value`` の形式でキーワードとその値を指定します.
また, #をつけることでコメントを入力することができます(それ以降の文字は無視されます).

キーワード
^^^^^^^^^^

    -  ``type``

       **形式 :** str型

       **説明 :**
       ソルバーの種類(``OpenMX, QE, VASP`` )を指定します.

    -  ``path``

       **形式 :** str型

       **説明 :**
       ソルバーへのパスを指定します.

    -  ``base_input_dir``

       **形式 :** str型

       **説明 :** 
       ベースとなる入力ファイルへのパスを指定します.

    -  ``run_scheme``

       **形式 :** str型

       **説明 :**
       ソルバーを起動する方法を指定します.
       詳細は :ref:`solver_specific_notes` を参照してください。

    -  ``perturb``

       **形式 :** float型

       **説明 :**
       対称性が良い構造を入力にしてしまうと、構造最適化が鞍点で止まってしまいがちである。これを避けるため、各原子をこのパラメータに比例するようにランダムに変位させたものを初期構造とする。0.0あるいはfalseに設定することも可能. デフォルト値 = 0.0.

    -  ``parallel_level`` (QuantumESPRESSO のみ)

       **形式 :** 辞書型

       **説明 :** 
       `Parallelization levels <https://www.quantum-espresso.org/Doc/user_guide/node18.html>`_ について、各level の並列数を指定します。
       長い形式のコマンドラインオプションから ``-`` を抜いたもの、すなわち、
       ``nimage``, ``npools``, ``nband``, ``ntg``, ``ndiag`` をキーとして、各level の分割数を値とする辞書として指定します。
       指定した要素のみ、実際のコマンドラインオプションとして ``pw.x`` に渡されます。

