.. highlight:: none

***************************
出力ファイルフォーマット
***************************

``RANK`` はプロセス番号＝レプリカ番号を表します。


``RANK/structure.XXX.vasp``
=================================
各ステップごとの原子配置が VASP の POSCAR ファイル形式で出力されます。
ステップ番号がファイル名の ``XXX`` に入ります。

例::

 Mg8 Al16 O32
 1.0
 8.113600 0.000000 0.000000
 0.000000 8.113600 0.000000
 0.000000 0.000000 8.113600
 Al Mg O
 16 8 32
 direct
 0.011208 0.995214 0.998158 Al
 0.758187 0.240787 0.499981 Al
 ... skipped ...
 0.746308 0.744706 0.233021 O
 0.257199 0.255424 0.771040 O

``RANK/minE.vasp``
==========================
最低エネルギーを与えた原子位置が VASP の POSCAR ファイル形式で出力されます。

``RANK/obs.dat``
=========================
各ステップごとの温度とエネルギーが電子ボルト単位で出力されます。

例::

 0	0.1034076	-41690.28269769395
 1	0.1034076	-41692.06763035158
 2	0.1034076	-41692.06763035158
 3	0.1034076	-41691.98205990787
 4	0.1034076	-41692.74143710456

``RANK/obs_save.npy``
==========================
各ステップごとのエネルギーが電子ボルト単位で出力されます。
``numpy.load('obs_save.npy')`` で、 ``darray`` として読み取ることができます。

例::

 $ python -c "import numpy; print(numpy.load('obs_save.npy'))"
 [[-41690.28269769]
  [-41692.06763035]
  [-41692.06763035]
  [-41691.98205991]
  [-41692.7414371 ]]

``RANK/kT_hist.npy``
========================
各ステップごとの温度（電子ボルト単位）が Numpy バイナリ形式で出力されます。
``numpy.load('kT_hist.npy')`` で、 ``darray`` として読み取ることができます。

例::

 $ python -c "import numpy; print(numpy.load('kT_hist.npy'))"
 [0.1034076 0.1034076 0.1034076 0.1034076 0.1034076]


``RANK/Trank_hist.npy``
===========================
（RXMC のみ）
各ステップごとの温度インデックスが Numpy バイナリ形式で出力されます。
``numpy.load('Trank_hist.npy')`` で、 ``darray`` として読み取ることができます。

例::

 $ python -c "import numpy; print(numpy.load('Trank_hist.npy'))"
 [1 1 1 1 1]


``RANK/logweight_hist.npy``
=============================
(PAMC のみ)
各ステップにおけるNeal-Jarzynski 重みの対数が Numpy バイナリ形式で出力されます。

Example::

 $ python -c "import numpy; print(numpy.load('logweight_hist.npy'))"
 [0 0 0 0 0]


``RANK/acceptance_ratio.dat``
================================
(PAMC のみ)
各温度におけるモンテカルロ更新の採択率。
1列目に温度, 2列目に採択率（採択回数/更新回数）が出力されます。


``logZ.dat``
==============
(PAMC のみ)
分配関数の対数 :math:`\log Z_i/Z_0` (:math:`i` は温度点の番号)。

- 1列目は温度 :math:`T_i`
- 2列目、3列目は分配関数 :math:`\log Z_i/Z_0` とその誤差
- 4列目、5列目は直前の温度との比 :math:`\log Z_i/Z_{i-1}` とその誤差

``result.dat``
=================
(PAMC のみ)
温度ごとの物理量 :math:`O` のカノニカル平均 :math:`\langle O \rangle` とその統計誤差 :math:`\sigma[O]` 。

- 1列目 は温度 :math:`T_i`
- 2列目、3列目はエネルギーの期待値 :math:`\langle E \rangle` と統計誤差
- 4列目、5列目はエネルギーの2乗の期待値 :math:`\langle E^2 \rangle` と統計誤差
- 6列目、7列目はエネルギーのゆらぎ :math:`\langle E^2 \rangle - \langle E \rangle^2` と統計誤差

   - エネルギーのゆらぎは熱容量 :math:`C` と次のようにして結びついています: :math:`k_B T^2 C = \left[ \langle E^2 \rangle - \langle E \rangle^2 \right]`


Potts ソルバーの場合、全磁化 :math:`\langle M \rangle = \langle \sum_i \delta_{\sigma_i,0} - 1/Q \rangle`,
全磁化の2乗 :math:`\langle M^2 \rangle`, 磁化のゆらぎ :math:`\langle M^2 \rangle - \langle M \rangle^2`
も8列目以降に出力されます。
