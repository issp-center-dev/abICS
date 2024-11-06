.. highlight:: none

[train] セクション
-------------------------------

訓練データから配置エネルギー予測モデルを学習する学習器の設定を行います。

予測モデルの作成・学習には外部のプログラムを利用します。
現在はaenet, NequIP, MLIP-3に対応しています。
ソフトウェア固有の注意事項（入力ファイル名など）は :ref:`trainer_specific_notes` を参照してください.

本セクションは以下のようなファイルフォーマットをしています.

.. code-block:: toml

    [train] # モデル学習器の設定
    type = 'aenet'
    base_input_dir = './aenet_train_input'
    ignore_species = ["O"]
    [train.exe_command]
    generate = '~/git/aenet/bin/generate.x-2.0.4-ifort_serial'
    train = 'srun ~/git/aenet/bin/train.x-2.0.4-ifort_intelmpi'

    
入力形式
^^^^^^^^^^^^
``keyword = value`` の形式でキーワードとその値を指定します.
また, #をつけることでコメントを入力することができます(それ以降の文字は無視されます).

キーワード
^^^^^^^^^^

 -  ``type``

    **形式 :** str型 
    
    **説明 :**
    訓練データから配置エネルギー予測モデルを学習する学習器の設定を行います.
    aenet, nequip, mlip_3を利用できます.


 -  ``base_input_dir``

    **形式 :** str型 

    **説明 :**
    設定したディレクトリの中に、学習器の設定ファイルを設置します.


 -  ``exe_command``

    **形式 :** 辞書型
    
    **説明 :**
    学習器で使う実行コマンドを指定します.
    コマンドライン引数も指定できますが, それぞれの学習機の入力ファイル (``input.yaml`` など)は含めないようにしてください.
    
    - ``type = 'aenet'``

      - ``generate`` と ``train`` の2つのキーを持ちます.
      - ``generate``

         - aenetの ``generate.x`` へのパスを指定します.

      - ``train``

         - aenetの ``train.x`` へのパスを指定します.
         - MPI並列版が利用可能です. その場合、上の例で示すように、MPI実行するためのコマンド（ ``srun`` 、 ``mpirun`` など）を合わせて設定してください。
      
      - abICS 2.0 以前との互換性のために、配列形式もサポートしています.
        最初の要素が ``generate``, 2番目の要素が ``train`` です.

    - ``type = 'nequip'``

      - ``train`` 

         - ``nequip-train`` へのパスを指定します.

    - ``type = 'mlip_3'``

      - ``train``

         - ``mlp`` へのパスを指定します.


 -  ``ignore_species``
   
    **形式 :** list型

    **説明 :**
    ``aenet`` などのニューラルネットワークモデルで「無視」する原子種を指定します. 常に占有率が1のものについては、ニューラルネットワークモデルの訓練および評価時に存在を無視した方が、計算効率が高くなります.
