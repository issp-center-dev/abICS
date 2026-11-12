# Additional setup and usage

[Beginner tutorial](../README.md) · [日本語](reference.ja.md)

Run these commands from the working copy prepared in Section 2 of the
[beginner tutorial](../README.md), for example `~/abics-work/spinel`.

## Optional packages

This example needs Python 3.12 or newer because of TorchSim and Orb. It does
not change the Python requirement or dependencies of abICS itself. Install
abICS from your local repository as shown in the beginner tutorial, then install this
example's optional requirements in a separate environment.

The requirements pin TorchSim 0.6.0, Orb 0.7.0, and their GPU-operation libraries.
TorchSim 0.6.1/0.6.2 require `nvalchemi-toolkit-ops>=0.4`, which conflicts with
Orb 0.7.0's `<0.4` requirement. The PyPI Warp 1.11.1 wheel supports CUDA 12;
newer Warp releases may require a newer GPU driver. See the
[versioned Warp installation guide](https://nvidia.github.io/warp/v1.11/user_guide/installation.html).
Install CUDA-enabled PyTorch before the requirements file, as in the README.

## Slurm

Save the commands you want to run in a shell script, for example `run.sh`, and
activate the Python and MPI environment in that script. Request eight MPI tasks
and one GPU:

```sh
sbatch --nodes=1 --ntasks=8 --cpus-per-task=1 --gres=gpu:1 \
    --time=04:00:00 --wrap='bash run.sh'
```

Use the GPU assigned by the scheduler. On a workstation with several GPUs,
select an available one before running the Python relaxation, for example
`export CUDA_VISIBLE_DEVICES=0`. MPS and compute-mode changes are unnecessary.
Only the abICS commands use `mpirun`; `run_umlp_torchsim.py` is one GPU process.

## Offline model loading

Orb may make an HTTP HEAD request to check its cache even if the weights were
already downloaded. To run without external HTTP requests, point to the local
checkpoint explicitly:

```sh
export ORB_WEIGHTS_PATH=/absolute/path/to/orb-v3-conservative-20-omat-20250404.ckpt
```

Use the checkpoint for `orb_v3_conservative_20_omat`. Record its checksum with
your results. Model loading with this option was verified with HTTP disabled.

## Optimization and output files

The TorchSim script processes up to eight structures per batch by default.
Reduce `--batch-size` if GPU memory is insufficient. It uses the Orb v3
conservative 20-neighbor OMAT model at `float32-highest` precision with model
compilation disabled. `TORCHDYNAMO_DISABLE=1` avoids an internal dynamic-shape
compilation error encountered in the tested environment.

Atomic positions and the full cell are relaxed with a Frechet cell filter.
FIRE runs for at most 200 steps at `fmax=0.2`, followed by BFGS for at most
1,000 steps at `fmax=0.04`. Both atomic and cell forces enter the final
convergence check. The FIRE state is converted to a plain simulation state
before BFGS initialization.

Each calculation directory contains:

| File | Meaning |
| --- | --- |
| `structure_norel.vasp` | Structure supplied to the relaxer |
| `structure.vasp` | Relaxed structure, after successful convergence |
| `energy.dat` | Relaxed total energy in eV, after successful convergence |
| `torchsim_relaxation.json` | Energy, force maxima, and convergence status |
| `structure_unconverged.vasp` | Last structure if the relaxation did not converge |

The overall summary records batch sizes, timings, and peak PyTorch-allocated
GPU memory. Reaching the optimizer's step limit does not count as convergence:
failed structures receive no `energy.dat`, and the command exits nonzero.
Existing reference outputs are rejected. Keep failed calculations for diagnosis
and use fresh input directories when retrying with different settings.

abICS maps the collected relaxed structures onto the reference lattice for
training. Check `train.out` for failed mappings. Oxygen participates in the
reference relaxations, but is omitted from the on-lattice descriptors because
this model uses the Mg/Al occupations as its input.

## ASE alternative and CPU execution

`run_umlp_ase.py` uses ASE instead of TorchSim and relaxes one structure per
process. It uses the same Orb 0.7.0 interface, precision, and Frechet filter.
The default device is CPU; add `--device cuda` to use the assigned NVIDIA GPU.
Use a fresh input copy, run the first `abics_mlref` command, and replace the
TorchSim relaxation command with:

```sh
while IFS= read -r directory; do
    python run_umlp_ase.py "$directory" || exit 1
done < rundirs.txt
```

Keep the collection, training, and sampling commands unchanged. The ASE script
also checks convergence before writing `energy.dat`. No GNU Parallel is needed
for sequential execution. Relaxation trajectories can differ between ASE and
TorchSim even with matching force tolerances.

## Other potentials

The reference solver is `user`: abICS prepares `structure.vasp` and later reads
`energy.dat` plus the relaxed `structure.vasp`. An ASE-compatible calculator
can replace Orb in `run_umlp_ase.py`; a TorchSim-compatible model can replace
it in the batched script. Install the selected potential's dependencies and
check units, stress support, species coverage, and convergence before use.

## Checking DOI convergence

Two cycles are an introduction to the workflow. Continue adding reference data
until model updates no longer change the DOI curve within your chosen tolerance.
Then hold the model fixed and test longer Monte Carlo runs with independent
random seeds. The supplied `DOI.dat` error column does not account for all time
correlations; use a blocking analysis of the temperature-resolved time series
for a convergence study. See the [validation example](validation.md#doi-convergence).

The [original standalone tutorial](https://github.com/skasamatsu/abics_umlp)
provided the starting spinel inputs and on-lattice workflow for this example.
