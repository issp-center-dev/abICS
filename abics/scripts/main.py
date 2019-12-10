import copy
import sys

from mpi4py import MPI
import numpy as np

from abics.mc import CanonicalMonteCarlo
from abics.mc_mpi import RX_MPI_init, TemperatureRX_MPI, RXParams
from abics.applications.latgas_abinitio_interface import default_observer
from abics.applications.latgas_abinitio_interface.model_setup import (
    dft_latgas,
    ObserverParams,
)
from abics.applications.latgas_abinitio_interface.defect import (
    defect_config,
    DFTConfigParams,
)
from abics.applications.latgas_abinitio_interface.run_base_mpi import runner
from abics.applications.latgas_abinitio_interface.vasp import VASPSolver
from abics.applications.latgas_abinitio_interface.qe import QESolver
from abics.applications.latgas_abinitio_interface.params import DFTParams


def main():
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    rxparams = RXParams.from_toml(tomlfile)
    nreplicas = rxparams.nreplicas
    nprocs_per_replica = rxparams.nprocs_per_replica

    kB = 8.6173e-5

    comm = RX_MPI_init(rxparams)

    # RXMC parameters
    # specify temperatures for each replica, number of steps, etc.
    kTstart = rxparams.kTstart
    kTend = rxparams.kTend
    kTs = kB * np.linspace(kTstart, kTend, nreplicas)

    # Set Lreload to True when restarting
    Lreload = rxparams.reload

    nsteps = rxparams.nsteps
    RXtrial_frequency = rxparams.RXtrial_frequency
    sample_frequency = rxparams.sample_frequency
    print_frequency = rxparams.print_frequency

    dftparams = DFTParams.from_toml(tomlfile)

    if dftparams.solver == 'vasp':
        solver = VASPSolver(dftparams.path)
    elif dftparams.solver == 'qe':
        solver = QESolver(dftparams.path)
    else:
        print('unknown solver: {}'.format(dftparams.solver))
        sys.exit(1)

    # model setup
    # we first choose a "model" defining how to perform energy calculations and trial steps
    # on the "configuration" defined below
    energy_calculator = runner(
        base_input_dir=dftparams.base_input_dir,
        Solver=solver,
        nprocs_per_solver=nprocs_per_replica,
        comm=MPI.COMM_SELF,
        perturb=dftparams.perturb,
        solver_run_scheme=dftparams.solver_run_scheme
    )
    model = dft_latgas(energy_calculator, save_history=False)


    # defect sublattice setup

    configparams = DFTConfigParams.from_toml(tomlfile)

    spinel_config = defect_config(configparams)

    configs = []
    for i in range(nreplicas):
        configs.append(copy.deepcopy(spinel_config))


    obsparams = ObserverParams.from_toml(tomlfile)

    # RXMC calculation
    RXcalc = TemperatureRX_MPI(comm, CanonicalMonteCarlo, model, configs, kTs)
    if Lreload:
        RXcalc.reload()
    obs = RXcalc.run(
        nsteps,
        RXtrial_frequency,
        sample_frequency,
        print_frequency,
        observer=default_observer(comm, Lreload),
        subdirs=True,
    )


    if comm.Get_rank() == 0:
        print(obs)
