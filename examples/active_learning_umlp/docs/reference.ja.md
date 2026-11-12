# 詳しい設定と使い方

[入門チュートリアル](../README.ja.md) · [English](reference.md)

コマンドは[入門チュートリアル](../README.ja.md)の第2節で用意した、
コピー先の作業ディレクトリ（例：`~/abics-work/spinel`）から実行します。

## 追加パッケージ

TorchSimとOrbのため、このexampleにはPython 3.12以上が必要です。abICS本体の
Python要件や依存関係は変更しません。入門チュートリアルのとおり専用環境を作り、
ローカルリポジトリからabICSを入れた後、このexampleの追加requirementsをインストールします。

TorchSim 0.6.0、Orb 0.7.0とGPU演算用ライブラリを固定しています。TorchSim 0.6.1/0.6.2の
`nvalchemi-toolkit-ops>=0.4`という要件は、Orb 0.7.0の`<0.4`と競合します。
PyPIのWarp 1.11.1はCUDA 12に対応し、新しい版では新しいGPUドライバが必要な場合が
あります。[該当版のWarpインストール案内](https://nvidia.github.io/warp/v1.11/user_guide/installation.html)
も参照してください。READMEのとおり、CUDA対応PyTorchを先に入れます。

## Slurm

実行したいコマンドを、例えば`run.sh`に保存し、その中でPythonとMPIの環境を
有効にしてください。MPI 8タスクとGPU 1枚を確保します。

```sh
sbatch --nodes=1 --ntasks=8 --cpus-per-task=1 --gres=gpu:1 \
    --time=04:00:00 --wrap='bash run.sh'
```

GPUはスケジューラの割り当てに従ってください。GPUが複数あるワークステーションでは、
例えば`export CUDA_VISIBLE_DEVICES=0`で空いているGPUを選んでから緩和を実行します。
MPSやcompute modeの変更は不要です。`mpirun`はabICSのコマンドに使い、
`run_umlp_torchsim.py`はGPUを使う1プロセスとして起動します。

## オフラインでのモデル読み込み

Orbは重みが保存済みでも、キャッシュ確認のためHTTP HEADを送る場合があります。
外部HTTP通信なしで実行するには、ローカルの重みを明示します。

```sh
export ORB_WEIGHTS_PATH=/absolute/path/to/orb-v3-conservative-20-omat-20250404.ckpt
```

`orb_v3_conservative_20_omat`に対応する重みを使い、チェックサムを結果と一緒に
記録してください。この指定でHTTP通信を禁止しても読み込めることを確認しています。

## 最適化と出力ファイル

TorchSim版は初期設定で最大8構造をまとめて処理します。GPUメモリが不足する場合は
`--batch-size`を減らしてください。Orb v3 conservative（最大20近傍、OMAT）を
`float32-highest`で使い、モデルのコンパイルは無効にしています。
`TORCHDYNAMO_DISABLE=1`は検証環境で起きた動的shapeの内部コンパイルエラーを回避する設定です。

Frechet cell filterで原子位置とセルを緩和します。FIREは`fmax=0.2`で最大200ステップ、
続いてBFGSを`fmax=0.04`で最大1,000ステップ実行します。原子とセル両方の力を収束判定に
含めます。FIREの状態はいったん通常のsimulation stateへ変換してからBFGSを初期化します。

各計算ディレクトリには次のファイルが保存されます。

| ファイル | 内容 |
| --- | --- |
| `structure_norel.vasp` | 緩和に渡した構造 |
| `structure.vasp` | 収束に成功した場合の緩和後の構造 |
| `energy.dat` | 収束に成功した場合の緩和後の全エネルギー（eV） |
| `torchsim_relaxation.json` | エネルギー、最大の力、収束状態 |
| `structure_unconverged.vasp` | 未収束だった場合の最後の構造 |

全体のsummaryにはバッチサイズ、時間、PyTorchが割り当てたGPUメモリの最大値を記録します。
最大ステップ数への到達は収束を意味しません。未収束の構造に`energy.dat`は書かず、
コマンドは非ゼロで終了します。参照出力がすでにある場合も停止します。失敗した計算は
診断用に残し、条件を変えてやり直す場合は新しい入力ディレクトリを使ってください。

abICSは回収した緩和後の構造を基準格子へマッピングして訓練します。マッピング失敗が
ないか`train.out`を確認してください。Oは参照構造の緩和に含めますが、on-latticeモデルは
Mg/Alの占有状態を入力とするため、Oをその記述子から除外します。

## ASE版とCPUでの実行

`run_umlp_ase.py`はTorchSimの代わりにASEを使い、1プロセスで1構造を緩和します。
Orb 0.7.0のインタフェース、精度、Frechet filterは同じです。初期設定はCPUです。
割り当てたNVIDIA GPUを使う場合は`--device cuda`を追加してください。
新しい入力コピーで最初の`abics_mlref`を実行し、TorchSimの緩和コマンドを次に置き換えます。

```sh
while IFS= read -r directory; do
    python run_umlp_ase.py "$directory" || exit 1
done < rundirs.txt
```

結果回収・訓練・サンプリングのコマンドは同じです。ASE版も収束を確認してから
`energy.dat`を書き出します。逐次実行にGNU Parallelは不要です。
同じ収束基準でも、ASEとTorchSimの緩和軌道は異なる場合があります。

## 別のポテンシャル

参照ソルバーは`user`です。abICSが`structure.vasp`を用意し、後で`energy.dat`と
緩和後の`structure.vasp`を読みます。`run_umlp_ase.py`のOrbをASE対応calculatorへ、
バッチ版のモデルをTorchSim対応モデルへ置き換えられます。使うポテンシャルの依存関係を
導入し、単位、応力への対応、元素の対応範囲、収束を確認してください。

## DOIの収束確認

2サイクルは手順を学ぶための例です。モデル更新によるDOI曲線の変化が選んだ許容幅に
収まるまで参照データを追加してください。その後はモデルを固定し、独立な乱数seedで
Monte Carloを長く実行して確認します。`DOI.dat`の誤差列は時間相関を十分に考慮している
とは限らないため、温度ごとの時系列をブロック解析してください。
[検証例](validation.ja.md#doiの収束確認)も参照してください。

このexampleのスピネル入力とon-lattice学習手順は、
[元のチュートリアル](https://github.com/skasamatsu/abics_umlp)を出発点にしています。
