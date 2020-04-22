.. highlight:: none

st2abics
-------------------------------

abICS :ref:`入力ファイル<input-format>` の :ref:`[config]セクション<config-section>` を比較的簡単に作成するために、 
``st2abics`` ツールを使うことができます。これはpymatgenで読める原子構造ファイルを読み取り、
[config]セクションを埋めたabICS入力テンプレートファイルに変換します。元の構造ファイルからどのようにして
config.base_structureとconfig.defect_structureを構成するか指定するため、 ``st2abics`` 専用の制御ファイルが必要です。
``st2abics`` は以下のように使用します::


    $ st2abics -h
    usage: st2abics [-h] inputfi structurefi [outfi]

    Prepare abICS config from structure file

    positional arguments:
      inputfi      toml input file for st2abics
      structurefi  Structure file that can be read by pymatgen Structure.from_file() method
      outfi        Output file to be used as abics input. Defaults to standard output

    optional arguments:
      -h, --help   show this help message and exit


例がいくつか ``examples/standard/st2abics`` 用意されています::

    $ cd examples/standard/st2abics
    $ st2abics st2abics_MgAl2O4.toml MgAl2O4.vasp abics_MgAl2O4.toml # spinel
    $ st2abics st2abics_CuZn.toml CuZn.vasp abics_CuZn.toml # brass
    $ st2abics st2abics_BZY.toml BaZrO3.vasp abics_BZY.toml # Y-doped BaZrO3

結果として得られたファイル(上記の例ではabics_MgAl2O4.toml, abics_CuZn.toml, abics_BZY.toml)は、
``[replica]`` 、``[solver]`` 、 ``[observer]`` セクションのキーワードが空になっており、
これらを記入した後、abICSの入力として使用することができます。

入力フォーマット
^^^^^^^^^^^^^^^^^
``st2abics`` の入力ファイルの例は、``examples/standard/st2abics`` にあります(上の例ではst2abics_CuZn.toml、st2abics_MgAl2O4.toml、st2abics_BZY.toml)。

フォーマットはabICS入力ファイルの[config]セクションに似ています。

キーワード
^^^^^^^^^^
-  ``supercell`` 

   **形式 :** list型 
   
   **説明 :** スーパーセルの大きさをリスト形式 [ :math:`\bf{a}, \bf{b}, \bf{c}` ] で指定します.

-  ``[[config.base_structure]]`` セクション

   ここでは、モンテカルロ計算中に格子サイト間で原子を交換しないbase_structureを指定します。

   -  ``species`` 

      **形式 :** 文字列のlist
      
      **説明 :** base_structureの原子種。対応する座標は入力構造ファイルから自動的に抽出されます。

   -  ``fix``
   
      **形式 :** bool型 ("true" or "false")
      
      **説明 :** base_structureの局所的な緩和を行わない場合true、行う場合はfalseに設定する。

-  ``[[config.defect_structure]]`` セクション

   このセクションは、配置サンプリングを行う副格子を指定します。
   複数の[[config.defect_structure]]セクションが存在しても構いません。
   例えば、陽イオンの配置サンプリングを行う副格子と、陰イオンの配置サンプリングを行う副格子を指定することが
   できます。
  
   -  ``site_center_species``

      **形式 :** 文字列のlist
      
      **説明 :** 元の構造ファイルに含まれる元素種のうち、配置サンプリングを行う格子サイトと対応するもの。

   -  ``[[config.defect_structure.groups]]`` サブセクション
   
      このセクションでは、配置サンプリングを行う格子サイト上に配置される原子グループを指定します。もしこのセクション
      がない場合は、 ``site_center_species`` を基に、入力された構造ファイルから自動的に構築されます。
      
      -  ``name``

         **形式 :** 文字列
         
         **説明 :** 原子グループの名前。

      -  ``species``
         
         **形式 :** 文字列のlist
         
         **説明 :** 原子グループに属する元素種のリスト。デフォルト値は、[``name``]です。

      -  ``coords``
      
         **形式 :** listのlistのlist あるいは文字列
         
         **説明 :** 原子グループがとることのできる配向ごとの原子グループ内の各原子の座標（入力ファイルのcoords定義 :ref:`参照 <coords-orr>`）。
         デフォルト値は[[[0.0, 0.0, 0.0]]です。

      -  ``num``
      
         **形式 :** int
         
         **説明 :** このセクションで指定した原子グループの数。スーパーセル内のサイト数に応じた数を指定してください。
