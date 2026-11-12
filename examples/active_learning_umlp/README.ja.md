# 機械学習ポテンシャルではじめるabICS

[English](README.md)

MgAl<sub>2</sub>O<sub>4</sub>スピネルのMgとAlの並び方は、温度によってどう変わるでしょうか。
このチュートリアルでは、用意された入力ファイルから計算を始め、その変化を表す曲線を求めます。

使うモデルは2つです。学習済みの汎用機械学習ポテンシャル **Orb** で、原子位置とセルを
動かして安定な構造に近づけ（構造緩和）、エネルギーを計算します。そのデータから
**ænet** で、固定した格子上の原子配置を入力として、**構造緩和後のエネルギー**を
すばやく予測する **on-latticeモデル**を訓練します。abICSはこのモデルを使い、
試行ごとに構造緩和を繰り返すことなく、Monte Carlo計算で各温度の原子配置を
サンプリングします。ここでのMonte Carlo計算では、格子上のMgとAlを入れ替える
試行を繰り返し、温度によってどのような配置が現れるかを調べます。

ここでは検証済みの **Linux・NVIDIA GPU 1枚・TorchSim** の手順に沿って進めます。
TorchSimは複数の構造をまとめてGPUで緩和するためのソフトウェアです。
CPUは8プロセス分、Pythonは3.12を使います。ASE（Atomic Simulation Environment）を
使って構造を1つずつ緩和する方法、CPUだけで実行する方法、別のポテンシャルを使う方法は
[詳細ガイド](docs/reference.ja.md)にまとめています。

## 1. 計算機を準備する

Git、MPI（abICSを8プロセスで動かすため）、CMake、Fortranコンパイラ、BLAS/LAPACK
（ænetをビルドするため）が必要です。Ubuntuでは次でインストールできます。

```sh
sudo apt install git cmake build-essential gfortran libblas-dev liblapack-dev openmpi-bin libopenmpi-dev
```

共同利用計算機では、これらがすでに用意されている場合があります。その場合は計算機の
案内に従ってMPIとコンパイラの環境をロードしてください。計算は **CPU 8タスク・GPU 1枚**を
確保したジョブ内で実行します。[クラスタでの実行例](docs/reference.ja.md#slurm)
にSlurmの指定例があります。

以下の手順ではcondaを使います。未導入なら、先に
[Miniforge](https://github.com/conda-forge/miniforge#install)をインストールしてください。

## 2. abICSをインストールし、作業ディレクトリにコピーする

以下は検証に使った **PyTorch 2.8.0のCUDA 12.6版（`cu126`）**をインストールする例です。
[PyTorch公式のインストール手順](https://pytorch.org/get-started/previous-versions/)の**v2.8.0**の項目と
[NVIDIAのCUDA互換性ガイド](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html)で、
GPUとドライバーに対応したCUDA版を確認し、必要に応じて`--index-url`を変更してください。
同梱のrequirementsに合わせて`torch==2.8.0`は維持します。PyTorchの版を変更する場合は、
他のパッケージとの互換性も確認する必要があります。

以下ではabICSを取得してインストールし、exampleを`~/abics-work/spinel`にコピーします。
コピー先には、まだ存在しないディレクトリを指定してください。場所は変更してもかまいません。
このexampleを含むabICSのローカル
リポジトリがすでにある場合は、`git clone`を省略し、`cd abICS`をその場所に読み替えてください。

```sh
git clone --branch develop https://github.com/issp-center-dev/abICS.git
cd abICS
abics_src="$PWD"
conda create -n abics_torchsim python=3.12
conda activate abics_torchsim
python -m pip install "$abics_src"
mkdir -p ~/abics-work
cp -r "$abics_src/examples/active_learning_umlp" ~/abics-work/spinel
cd ~/abics-work/spinel
python -m pip install 'torch==2.8.0' --index-url https://download.pytorch.org/whl/cu126
python -m pip install -r requirements-torchsim.txt
```

以降もこの環境を有効にし、コピー先の作業ディレクトリ（`~/abics-work/spinel`）から
コマンドを実行してください。`abics_src`にはローカルのabICSリポジトリのフルパスを保存し、
`python -m pip install "$abics_src"`でそこからabICSをインストールします。
入力ファイルの編集と計算結果の保存は、コピー先で行います。
requirementsファイルには、このチュートリアルに追加で必要なパッケージの版を指定しています。
版を固定している理由は[詳細ガイド](docs/reference.ja.md#追加パッケージ)で説明しています。

PythonからGPUが見えるか確認します。

```sh
python -c 'import torch; print(torch.cuda.is_available())'
```

`True`と表示されれば準備できています。クラスタではGPUを確保したジョブ内で確認します。
`False`の場合は、GPUを利用できる状態にしてから先へ進んでください。

## 3. ænetをビルドする

ænetの3つのプログラムのうち、`generate`は訓練用データの準備、`train`はモデルの訓練、
`predict`はサンプリング中のエネルギー予測に使います。
作業ディレクトリから、`../aenet`（上記の手順では`~/abics-work/aenet`）にænetを取得し、
CMakeでこれらのプログラムをコンパイルします。

```sh
git clone https://github.com/atomisticnet/aenet.git ../aenet
cmake -S ../aenet -B ../aenet/build -DBUILD_AENET=ON -DUSE_MPI=OFF -DCMAKE_Fortran_COMPILER=gfortran
cmake --build ../aenet/build --target main --parallel 8
```

すでにænetをインストールしている場合は、その実行ファイルを使うこともできます。
次のコマンドで、ビルドしたプログラムがある`build/bin`のフルパスを表示します。

```sh
realpath ../aenet/build/bin
```

表示されたパスを、次の設定で使います。別のコンパイラを使う場合やビルドでエラーが
出た場合は、[ænetのインストール手順](https://github.com/atomisticnet/aenet#installation)を参照してください。

## 4. ænetのパスを設定し、入力ファイルを生成する

テキストエディタで`input.toml`を開き、下に示す3つの実行ファイルのパスを変更します。
`/your/path/aenet/build/bin`の部分は、先ほど表示したフルパスに置き換えてください。

```toml
[sampling.solver]
path = '/your/path/aenet/build/bin/predict.x_serial'

[train.exe_command]
generate = '/your/path/aenet/build/bin/generate.x_serial'
train = '/your/path/aenet/build/bin/train.x_serial'
```

**すでにあるセクションの中の値**を変更し、他の設定は残してください。
温度、組成、並列に動かすシミュレーションの数（**レプリカ数**）は、用意された設定のままで
実行できます。サンプリング中は、各レプリカが1つの温度を担当します。
このGPUを使う手順では、Pythonスクリプトを編集する必要はありません。

`input.toml`からænetの入力ファイルを生成します。

```sh
abics_gen_aenet_input input.toml
```

`aenet_train_input/`が作られ、AlとMgの記述子・訓練・予測の設定が保存されます。
用意された入力の`train.ignore_species`により、酸素はon-latticeモデルから除外されます。
Orbによる参照構造の緩和では酸素も含めて計算します。
このコマンドは最初の学習サイクルの前に1回だけ実行し、生成されたディレクトリは
以降のサイクルでも使います。既存のファイルがある場合は上書きせずに停止するため、
後から設定を調整した場合もその内容を保持できます。

## 5. 最初の学習サイクルを実行する

ここでいう「学習サイクル」は、構造の用意、Orbによる構造緩和とエネルギー計算、
on-latticeモデルの訓練、そのモデルを使ったサンプリング、という一連の流れです。
Orbで求めたエネルギーを、学習の正解となる**参照エネルギー**として使います。
このサイクルを繰り返して訓練データを追加する方法を**能動学習（active learning）**と呼びます。

以下は上から順に、各コマンドが正常に終了するのを待って実行してください。
`mpirun -n 8`は、協調して動く8プロセスでabICSを実行する指定です。
構造緩和と訓練のコマンドは、記載どおり1回ずつ起動します。
末尾の`>> ファイル名 2>&1`は、メッセージとエラーをログファイルに保存する指定なので、
端末には何も表示されない場合があります。同じコマンドを再実行すると、ログに追記されます。

実行前に、同じ端末またはジョブスクリプトで次を設定します。

```sh
export OMP_NUM_THREADS=1
export TORCHDYNAMO_DISABLE=1
```

### 40個の構造を用意する

```sh
mpirun -n 8 abics_mlref input.toml >> mlref.out 2>&1
```

8レプリカそれぞれに5構造、合計40構造を用意します。
最初の参照計算を保存する`AL0/`が作られ、計算ディレクトリの一覧が`rundirs.txt`に
保存されます。番号は0から始まるため、次のサイクルでは`AL1/`を使います。
次のGPU計算スクリプトは、`rundirs.txt`を自動で読み込みます。

### GPUで構造を緩和し、エネルギーを求める

```sh
python run_umlp_torchsim.py --batch-size 8 --summary torchsim-AL0.json >> torchsim.out 2>&1
```

このコマンドは **Python 1プロセス**で実行します。GPU上で8構造ずつまとめて処理します。
学習済みOrbモデルのパラメータは初回にダウンロードされる場合があり、
その際はインターネット接続が必要です。正常に終了すると、各計算ディレクトリに
緩和後の`structure.vasp`と`energy.dat`が保存されます。
別の端末で`tail -f torchsim.out`を実行すると進行状況を確認できます。

失敗した場合は、先へ進む前に`torchsim.out`を確認してください。
エネルギーを訓練に使うには、構造緩和で力が許容値より小さくなる必要があります。
1構造でもこの収束条件を満たさなければ、スクリプトはエラー終了します。

### エネルギーを回収する

```sh
mpirun -n 8 abics_mlref input.toml >> mlref.out 2>&1
```

最初と同じコマンドですが、今度は参照計算が終わっているため、abICSが結果を読み込み、
訓練データを用意します。

### On-latticeモデルを訓練する

```sh
abics_train input.toml >> train.out 2>&1
```

訓練済みモデルが`baseinput/`に用意されます。訓練の進行状況は`train.out`と
`train0/stdout`に記録されます。

### 温度依存性を計算する

```sh
mpirun -n 8 abics_sampling input.toml >> sampling.out 2>&1
```

600〜2,000 Kの8温度でMonte Carlo計算を行います。
最初のサイクルのMonte Carlo結果は`MC0/`に保存されます。
`MC0/DOI.dat`には、MgとAlの混ざり方を表す反転度（DOI）が記録されます。
反転度の読み方は第7節で説明します。

## 6. 2回目のサイクルでモデルを改善する

最初のモデルはランダムなMgとAlの配置で学習しました。次は、最初のMonte Carlo計算で
現れた配置を選び、サンプリング中に訪れる配置のデータを追加します。
同じディレクトリで、第5節の5つの手順を繰り返してください。

```sh
mpirun -n 8 abics_mlref input.toml >> mlref.out 2>&1
python run_umlp_torchsim.py --batch-size 8 --summary torchsim-AL1.json >> torchsim.out 2>&1
mpirun -n 8 abics_mlref input.toml >> mlref.out 2>&1
abics_train input.toml >> train.out 2>&1
mpirun -n 8 abics_sampling input.toml >> sampling.out 2>&1
```

追加の参照計算は`AL1/`、新しいMonte Carloの結果は`MC1/`に保存されます。
各コマンドが正常に終わってから、次のコマンドを実行してください。

## 7. 結果を読む

```sh
cat MC1/DOI.dat
```

先頭の3列は次の内容です。

| 列 | 内容 |
| --- | --- |
| 1 | 温度（K） |
| 2 | 反転度（DOI）の平均 |
| 3 | 平均の推定誤差 |

反転度は、正常スピネルでMgが占める8サイトに、どれだけAlが入っているかを表します。
0は反転なし、1は完全な反転です。`MC0/DOI.dat`と`MC1/DOI.dat`の1列目・2列目を
比較し、訓練データを追加したことで曲線がどう変わったかを見てみましょう。

**2サイクルは手順を学ぶための例で、収束を確認したことにはなりません。**
信頼できる曲線を得るには、学習サイクルを追加してもDOIの曲線がほとんど変わらないことを
確認します。また、乱数の種を変えて長いMonte Carlo計算を行い、結果が統計誤差の範囲で
一致することも確認します。別途実施した7サイクルの計算について、これらの確認方法と
収束を判断した数値基準を示しています。詳しくは[収束確認の例](docs/validation.ja.md#doiの収束確認)を参照してください。

## うまくいかないとき

| 症状 | 確認すること |
| --- | --- |
| GPU確認で`False`になる | 環境を有効にし、GPUを確保したジョブ内で確認する。 |
| MPIのスロット不足と言われる | `input.toml`と対応する8 MPIタスクを確保する。 |
| ænetのプログラムが見つからない | `input.toml`の3つのフルパスを確認する。 |
| 緩和が失敗する | `torchsim.out`と該当ディレクトリの`torchsim_relaxation.json`を確認する。 |
| 参照計算の出力がすでにあると言われる | 結果を保存し、別の計算は入力ファイルの新しいコピーで始める。 |

`AL*`、`MC*`、進行管理ファイルはまとめて残してください。abICSはこれらを使って
計算の段階を判断します。学習を続けるときは、第5節を最初からやり直すのではなく、
第6節に従ってください。独立した別の計算を始めるときは、ローカルのabICSリポジトリ内の
`examples/active_learning_umlp`を新しい作業ディレクトリにコピーし、
その`input.toml`にænetのパスを設定します。
以前の計算ディレクトリはそのまま残してください。

オフラインでのモデル読み込み、別のポテンシャル、最適化の設定、詳しい検証結果は
[詳細ガイド](docs/reference.ja.md)にまとめています。
