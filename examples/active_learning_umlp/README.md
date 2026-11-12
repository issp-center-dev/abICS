# Learn abICS with a machine-learning potential

[日本語](README.ja.md)

How does the arrangement of Mg and Al atoms in MgAl<sub>2</sub>O<sub>4</sub>
spinel change with temperature? In this tutorial, you will calculate a curve
showing that change, starting from the supplied input files.

You will use two models. **Orb**, a pretrained universal machine-learning
potential, lets atoms and the cell move toward a low-energy structure
("relaxation") and calculates its energy. Using these data, **ænet** learns an
**on-lattice model** that predicts the **relaxed energy** from the atomic
arrangement on a fixed lattice. abICS uses that model to sample arrangements
at different temperatures using Monte Carlo, without repeating structure
relaxation at every trial move. Here, Monte Carlo sampling means repeatedly
trying swaps of Mg and Al on lattice sites to find how their arrangements
vary with temperature.

This guide follows one tested route: **Linux, one NVIDIA GPU, and TorchSim**,
which relaxes several structures together on the GPU. You need eight CPU slots
and Python 3.12. For an alternative that relaxes one structure at a time using
ASE (Atomic Simulation Environment), instructions for CPU execution, or other
potentials, see the [detailed reference](docs/reference.md).

## 1. Prepare your computer

You need Git, MPI (to run eight abICS processes), CMake, a Fortran compiler,
and BLAS/LAPACK libraries (to build ænet). On Ubuntu, install them with:

```sh
sudo apt install git cmake build-essential gfortran libblas-dev liblapack-dev openmpi-bin libopenmpi-dev
```

On a shared cluster, these tools may already be provided. Load the site's MPI
and compiler environment instead. Run the calculations inside an allocation
with **eight CPU tasks and one GPU**; the [cluster notes](docs/reference.md#slurm)
give a Slurm example.

The commands below use conda. If you do not have it, install
[Miniforge](https://github.com/conda-forge/miniforge#install) first.

## 2. Install abICS and copy the tutorial to a working directory

The commands below use **PyTorch 2.8.0 with CUDA 12.6 (`cu126`)**, the build
used for this tutorial's validation. Check the **v2.8.0** section of the
[official PyTorch installation instructions](https://pytorch.org/get-started/previous-versions/)
and [NVIDIA's CUDA compatibility guide](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html)
to choose a build compatible with your GPU and driver, and adjust `--index-url`
if needed. Keep `torch==2.8.0` to match the supplied requirements; changing the
PyTorch version also requires checking compatibility with the other packages.

The commands below download abICS, install it, and copy the example to
`~/abics-work/spinel`. Use a destination directory that does not already exist.
You can choose another working directory. If you already
have a local abICS repository containing this example, skip `git clone` and
replace `cd abICS` with the path to that repository.

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

Keep this environment active. Run all remaining commands from the copied
working directory (`~/abics-work/spinel`). The `abics_src` variable holds the
full path to your local abICS repository; `python -m pip install "$abics_src"`
installs abICS from there. Input edits and calculation results go into the
working copy. The requirements file selects only the
additional packages for this example; the [reference](docs/reference.md#optional-packages)
explains the version constraints.

Check that Python can see your GPU:

```sh
python -c 'import torch; print(torch.cuda.is_available())'
```

The expected output is `True`. On a cluster, run this check inside the GPU job.
If it prints `False`, resolve GPU access before continuing.

## 3. Build ænet

ænet provides three programs: `generate` prepares data for training, `train`
fits the model, and `predict` evaluates its energies during sampling.
From the working directory, download ænet into `../aenet`
(`~/abics-work/aenet` with the paths above)
and compile these programs with CMake:

```sh
git clone https://github.com/atomisticnet/aenet.git ../aenet
cmake -S ../aenet -B ../aenet/build -DBUILD_AENET=ON -DUSE_MPI=OFF -DCMAKE_Fortran_COMPILER=gfortran
cmake --build ../aenet/build --target main --parallel 8
```

If ænet is already installed, you can use that installation instead.
Print the full path to the programs built in `build/bin`:

```sh
realpath ../aenet/build/bin
```

Copy the printed path for the next step. For other compilers or build errors,
see the [ænet installation instructions](https://github.com/atomisticnet/aenet#installation).

## 4. Set the ænet paths and generate its input files

Open `input.toml` in a text editor. Replace the three existing executable paths
shown below. Use the full path printed above in place of `/your/path/aenet/build/bin`:

```toml
[sampling.solver]
path = '/your/path/aenet/build/bin/predict.x_serial'

[train.exe_command]
generate = '/your/path/aenet/build/bin/generate.x_serial'
train = '/your/path/aenet/build/bin/train.x_serial'
```

Change these entries **inside the existing sections** and keep the other
settings. The supplied temperatures, composition, and number of parallel
simulations (called **replicas**) are ready to use. During sampling, each
replica runs at one temperature. No changes to the Python scripts are needed
for these GPU instructions.

Generate the ænet input files from `input.toml`:

```sh
abics_gen_aenet_input input.toml
```

This creates `aenet_train_input/` with the descriptor, training, and prediction
settings for Al and Mg. Oxygen is excluded from the on-lattice model by
`train.ignore_species` in the supplied input; Orb still includes oxygen during
reference relaxation. Run this once before the first learning cycle and keep
the generated directory for subsequent cycles. The command refuses to overwrite
existing files, so you can safely keep any settings you later customize.

## 5. Run your first learning cycle

A *learning cycle* prepares structures, calculates their relaxed energies with
Orb (the **reference energies** used as training targets), trains the
on-lattice model, and samples arrangements with that model. Repeating this
cycle to add training data is called **active learning**.

Run the commands below in order, waiting for each to finish successfully.
The `mpirun -n 8` prefix starts eight cooperating processes for abICS.
The relaxation and training commands are launched once, as shown.
The `>> filename 2>&1` suffix saves messages and errors to a log file, so
these commands may show no output in the terminal. It appends to the file
when the same command is used again.

Before starting, set the following in the same terminal or job script:

```sh
export OMP_NUM_THREADS=1
export TORCHDYNAMO_DISABLE=1
```

### Prepare 40 structures

```sh
mpirun -n 8 abics_mlref input.toml >> mlref.out 2>&1
```

This prepares five structures for each of eight replicas, giving 40 in total.
It creates `AL0/` for the first set of reference calculations and writes their
directory paths to `rundirs.txt`. Numbering starts at 0: the next cycle uses
`AL1/`. The GPU script in the next step reads `rundirs.txt` automatically.

### Calculate their relaxed energies on the GPU

```sh
python run_umlp_torchsim.py --batch-size 8 --summary torchsim-AL0.json >> torchsim.out 2>&1
```

Run this command as **one Python process**. It handles eight structures together
on the GPU. The pretrained Orb model parameters may be downloaded on first
use, which requires an internet connection. A successful run
writes a relaxed `structure.vasp` and an `energy.dat` in each calculation
directory. You can watch progress in another terminal with `tail -f torchsim.out`.

If this command fails, check `torchsim.out` before proceeding. The script stops
with an error if any structure has not converged: the forces must fall below
the relaxation tolerance before its energy can be used for training.

### Collect the energies

```sh
mpirun -n 8 abics_mlref input.toml >> mlref.out 2>&1
```

This is the same command as the first step. Now that the reference calculations
are complete, abICS reads their results and prepares the training data.

### Train the on-lattice model

```sh
abics_train input.toml >> train.out 2>&1
```

This creates the trained model in `baseinput/`. Training progress is recorded
in `train.out` and `train0/stdout`.

### Calculate the temperature dependence

```sh
mpirun -n 8 abics_sampling input.toml >> sampling.out 2>&1
```

This runs Monte Carlo at eight temperatures from 600 to 2,000 K. The result is
written to `MC0/DOI.dat`; `MC0/` holds the first cycle's Monte Carlo results.
The degree of inversion (DOI), explained in Section 7, measures Mg/Al mixing.

## 6. Improve the model with a second cycle

The first model learned from random Mg/Al arrangements. abICS now selects
arrangements visited during the first Monte Carlo run, so the added training
data cover configurations encountered during sampling. In the same directory,
repeat the five steps from Section 5:

```sh
mpirun -n 8 abics_mlref input.toml >> mlref.out 2>&1
python run_umlp_torchsim.py --batch-size 8 --summary torchsim-AL1.json >> torchsim.out 2>&1
mpirun -n 8 abics_mlref input.toml >> mlref.out 2>&1
abics_train input.toml >> train.out 2>&1
mpirun -n 8 abics_sampling input.toml >> sampling.out 2>&1
```

The new reference calculations go into `AL1/`, and the new Monte Carlo results
into `MC1/`. Wait for each command to finish successfully before running the next.

## 7. Read your result

```sh
cat MC1/DOI.dat
```

The first three columns are:

| Column | Meaning |
| --- | --- |
| 1 | Temperature in K |
| 2 | Mean degree of inversion (DOI) |
| 3 | Estimated error of the mean |

DOI measures how often Al occupies the eight sites assigned to Mg in the normal
spinel structure. A value of 0 means no inversion; 1 means full inversion.
Compare columns 1 and 2 in `MC0/DOI.dat` and `MC1/DOI.dat` to see how adding
training data changes the curve.

**Two cycles demonstrate the workflow; they do not establish convergence.**
For a reliable curve, check that adding learning cycles no longer changes the
DOI curve appreciably, and that longer Monte Carlo runs started with different
random seeds give consistent results within their statistical errors.
A separate extended calculation with seven learning cycles illustrates these
checks and the numerical criteria used to judge convergence. See the
[convergence example](docs/validation.md#doi-convergence).

## If something goes wrong

| Symptom | What to check |
| --- | --- |
| GPU check prints `False` | Activate the environment and check inside your allocated GPU job. |
| MPI reports too few slots | Allocate eight MPI tasks, matching `input.toml`. |
| ænet program not found | Check the three full executable paths in `input.toml`. |
| Relaxation fails | Read `torchsim.out` and the affected directory's `torchsim_relaxation.json`. |
| Script reports existing reference outputs | Keep those results; start a separate calculation in a fresh copy of the input files. |

Keep the `AL*`, `MC*`, and progress files together. abICS uses them to track the
workflow. To continue learning, follow Section 6 rather than starting again
at Section 5. To start an independent calculation, copy
`examples/active_learning_umlp` from your local abICS repository to a new working
directory and set the ænet paths in its `input.toml`. Keep the previous
calculation directory intact.

For offline model loading, other potentials, optimization settings, and complete
verification results, see the [detailed reference](docs/reference.md).
