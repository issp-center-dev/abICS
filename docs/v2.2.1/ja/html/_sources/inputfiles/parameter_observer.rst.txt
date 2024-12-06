.. highlight:: none

[observer] セクション
-------------------------------

計算する物理量を指定します.
以下のようなファイルフォーマットをしています.

  :: 

    [observer]
    [observer.similarity]
    reference_structure = './MgAl2O4.vasp'
    ignored_species = ["O"]

    [[observer.solver]]
    name = "magnetization"
    type = 'aenet'
    path= '~/opt/aenet/bin/predict.x_serial'
    base_input_dir = './baseinput_mag'
    perturb = 0.0
    run_scheme = 'subprocess'
    ignore_species = ["O"]


入力形式
^^^^^^^^^^^^
``keyword = value`` の形式でキーワードとその値を指定します.
また, #をつけることでコメントを入力することができます(それ以降の文字は無視されます).

キーワード
^^^^^^^^^^

- ``[[observer.solver]]``

  物理量を計算するためのオプションを指定します.
  本セクションは複数指定することができます.
  ``name`` を除き, ``sampling.solver`` セクションと同様の形式で指定します.

  なお、 ``name = "energy"`` という物理量については、 ``sampling.solver`` セクションで指定したものが自動的に適用されます.

  - ``name``

    **形式 :** str型 

    **説明 :**
    物理量の名前を指定します.
    計算終了後、温度ごとの期待値として, ``<name>.dat`` というファイルが出力されます.

- ``[observer.similarity]``

  原子配置の「類似度(similarity)」を計算するためのオプションを指定します.
  類似度は元素種ごとに、参照状態と同じ場所にある原子の割合として計算されます.
  計算終了後、温度ごとの期待値として ``similarity_X.dat`` というファイルが出力されます (``X`` は元素記号).
  本サブセクションが指定されていない場合は、類似度は計算されません.

  - ``reference_structure``

    **形式 :** str型 

    **説明 :**
    参照状態の構造ファイルを指定します.

  - ``ignored_species``

    **形式 :** list型 

    **説明 :**
    類似度を計算する際に無視する原子種を指定します.
    例えば, 酸素を固定した計算の場合には、このキーワードに ``["O"]`` を指定します.
