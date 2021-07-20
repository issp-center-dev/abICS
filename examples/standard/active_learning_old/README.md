1. Modify the following paths and commands to the aenet executables in the `input_aenet.toml`
  - `path` in `[solver]`
    - Path to aenet predictor
  - `path` in `[solverRef]`
    - Path to energy calculator for making training dataset
  - `exe_command` in `[trainer]`
    - aenet generator and aenet trainer
2. Run `mpiexec -np 4 abics_activelearn input_aenet.toml` to calculate energies for existing configurations (samples)
  - At first time, `replicaRef.nsteps` configurations will be generated randomly
  - In this example, energies are calculated by using the aenet solver with pretrained potential.
3. Run `abics_train input_aenet.toml` to train aenet potential
4. Run `mpiexec -np 4 abics input_aenet.toml` to perform RXMC calculation by using aenet energy solver
5. Repeat 2-4 until the calculation converges
