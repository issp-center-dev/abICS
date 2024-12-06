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

.. - その他
..
..     -  ``output_frequency``
..
..        **形式 :** list型 (自然数)
..
..        **説明 :**  ``nsteps`` をRXMC計算で出力される配置の数（ ``[replica]`` セクションの ``nsteps/sample_frequency`` の値）のうち何ステップまでを取り出すかのステップ数、 ``sample_frequency`` を配置を抜き出してくる間隔として、 [ ``nsteps`` , ``sample_frequency`` ] のlist形式で配置を抜き出す間隔を指定します(抜き出す最初のステップ数は0に固定しています)。その配置に対応する第一原理計算ソルバーの入力ファイルが各フォルダ内に作成されます。

