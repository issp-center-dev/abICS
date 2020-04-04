# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

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
from abics.applications.latgas_abinitio_interface.aenet import aenetSolver
from abics.applications.latgas_abinitio_interface.openmx import OpenMXSolver
from abics.applications.latgas_abinitio_interface.params import DFTParams


def main_impl(tomlfile):
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
    elif dftparams.solver == 'aenet':
        solver = aenetSolver(dftparams.path)
    elif dftparams.solver == 'openmx':
        solver = OpenMXSolver(dftparams.path)
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


def main():
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    main_impl(tomlfile)
