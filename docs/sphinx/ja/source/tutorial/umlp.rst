汎用機械学習ポテンシャルを使う
=================================

``examples/active_learning_umlp/`` では、学習済みのOrbで計算した参照エネルギーから、
MgAl2O4スピネルのon-latticeモデルを訓練します。TorchSimで複数構造をGPU 1枚に
まとめて処理し、ASE版ではCPUまたはGPUで逐次緩和もできます。

インストール、入力の準備、2回のactive learningサイクル、反転度の読み方は
`exampleの入門ガイド
<https://github.com/issp-center-dev/abICS/blob/develop/examples/active_learning_umlp/README.ja.md>`_
に沿って進めてください。GPU用の追加パッケージはabICSと別にインストールし、
Python 3.12以上を使います。
