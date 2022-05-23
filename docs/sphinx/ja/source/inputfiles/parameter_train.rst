.. highlight:: none

[train] セクション
-------------------------------

訓練データから配置エネルギー予測モデルを学習する学習器の設定を行います。
現在のところ、abICSではaenetのみに対応しています。
以下のようなファイルフォーマットをしています.

.. code-block:: toml

    [trainer] # モデル学習器の設定
    type = 'aenet'
    base_input_dir = './aenet_train_input'
    exe_command = ['~/git/aenet/bin/generate.x-2.0.4-ifort_serial', 
                  'srun ~/git/aenet/bin/train.x-2.0.4-ifort_intelmpi']
    ignore_species = ["O"]

    
入力形式
^^^^^^^^^^^^
``keyword = value`` の形式でキーワードとその値を指定します.
また, #をつけることでコメントを入力することができます(それ以降の文字は無視されます).

キーワード
^^^^^^^^^^

 -  ``type``

    **形式 :** str型 
    
    **説明 :**
    訓練データから配置エネルギー予測モデルを学習する学習器の設定を行います.現在のところ、abICSではaenetのみに対応しています.


 -  ``base_input_dir``

    **形式 :** str型 

    **説明 :**
    設定したディレクトリの中に、学習器の設定ファイルを設置します.


 -  ``exe_command``

    **形式 :** strのlist型 
    
    **説明 :**
    aenetの ``generate.x`` と ``train.x`` へのパスを指定します。 ``train.x`` についてはMPI並列版が利用可能で、その場合は、上の例で示すように、MPI実行するためのコマンド（ ``srun`` 、 ``mpirun`` など）を合わせて設定してください。


 -  ``ignore_species``
   
    **形式 :** list型

    **説明 :**
    ``aenet`` などのニューラルネットワークモデルで「無視」する原子種を指定します. 常に占有率が1のものについては、ニューラルネットワークモデルの訓練および評価時に存在を無視した方が、計算効率が高くなります.
