.. pyMC documentation master file, created by
   sphinx-quickstart on Wed Jul 31 13:13:22 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

概要
------------------------------------------
abICSは、第一原理計算を再現する機械学習モデルを訓練し、
不規則系での統計熱力学サンプリングを高速に
実行するためのソフトウェアフレームワークです。
金属や酸化物合金などの多成分固体系に特に重点を置いています。
現在はaenetで実装されているニューラルネットワークポテンシャルを
機械学習モデルとして利用することができます。
機械学習の基となる第一原理計算用入力ファイルの自動生成にも対応しており、
Quantum Espresso, VASP, OpenMXを利用することができます。
サンプリングアルゴリズムは、拡張モンテカルロ法であるレプリカ交換モンテカルロ法(RXMC)とポピュレーションアニーリング・モンテカルロ法(PAMC)を実装しています。また、β版としてグランドカノニカルサンプリングに対応しています。


開発者
------------------------------------------
abICSは以下のメンバーで開発しています.

- ver. 2.0-
   - 笠松 秀輔 (山形大学 学術研究院(理学部主担当))
   - 本山 裕一 (東京大学 物性研究所)
   - 青山 龍美 (東京大学 物性研究所)
   - 吉見 一慶 (東京大学 物性研究所)
   - 杉野 修 (東京大学 物性研究所)

- ver. 1.0
   - 笠松 秀輔 (山形大学 学術研究院(理学部主担当))
   - 本山 裕一 (東京大学 物性研究所)
   - 吉見 一慶 (東京大学 物性研究所)
   - 山本 良幸 (東京大学 物性研究所)
   - 杉野 修 (東京大学 物性研究所)
   - 尾崎 泰助 (東京大学 物性研究所)

   
バージョン履歴
------------------------------------------

- ver.2.2.1    : 2024/12/06.
- ver.2.2.0    : 2024/11/07.
- ver.2.1.0    : 2023/06/12.
- ver.2.0.1    : 2022/11/04.
- ver.2.0      : 2022/06/24.
- ver.1.0      : 2020/05/01.
- ver.1.0-beta : 2020/03/31.
- ver.0.1      : 2019/12/09.


ライセンス
--------------
本ソフトウェアのプログラムパッケージおよびソースコード一式はGNU General Public License version 3 (GPL v3) に準じて配布されています。

abICSを引用する際は, 以下の文献を引用してください。

Shusuke Kasamatsu, Yuichi Motoyama, Kazuyoshi Yoshimi, Tatsumi Aoyama, “Configuration sampling in multi-component multi-sublattice systems enabled by ab Initio Configuration Sampling Toolkit (abICS)”, `accepted in STAM: Methods <https://doi.org/10.1080/27660400.2023.2284128>`_ (`arXiv:2309.04769 <https://arxiv.org/abs/2309.04769>`_).

Bibtex::

   @article{kasamatsu2023configuration,
   author = {Shusuke Kasamatsu, Yuichi Motoyama, Kazuyoshi Yoshimi and Tatsumi Aoyama},
   title = {Configuration sampling in multi-component multi-sublattice systems enabled by ab initio Configuration sampling toolkit ({abICS})},
   journal = {Science and Technology of Advanced Materials: Methods},
   volume = {0},
   number = {ja},
   pages = {2284128},
   year = {2023},
   publisher = {Taylor & Francis},
   doi = {10.1080/27660400.2023.2284128},
   URL = {https://doi.org/10.1080/27660400.2023.2284128},
   eprint = {https://doi.org/10.1080/27660400.2023.2284128}

コピーライト
------------------

*(c) 2019- The University of Tokyo. All rights reserved.*

本ソフトウェアは2019, 2022年度 東京大学物性研究所 ソフトウェア高度化プロジェクトの支援を受け開発されており、その著作権は東京大学が所持しています。
