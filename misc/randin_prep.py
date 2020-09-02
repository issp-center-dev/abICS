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
import sys,os

import numpy as np
import scipy.constants as constants

# from abics.mc import CanonicalMonteCarlo, RandomSampling
# from abics.mc_mpi import (
#     RX_MPI_init,
#     TemperatureRX_MPI,
#     RXParams,
#     SamplerParams,
#     ParallelRandomParams,
#     EmbarrassinglyParallelSampling,
# )
# from abics.applications.latgas_abinitio_interface import default_observer
from abics.applications.latgas_abinitio_interface.model_setup import (
    dft_latgas,
    ObserverParams,
)
from abics.applications.latgas_abinitio_interface.defect import (
    defect_config,
    DFTConfigParams,
)
#from abics.applications.latgas_abinitio_interface.run_base_mpi import runner, runner_multistep
from abics.applications.latgas_abinitio_interface.vasp import VASPSolver
from abics.applications.latgas_abinitio_interface.qe import QESolver
from abics.applications.latgas_abinitio_interface.aenet import aenetSolver
from abics.applications.latgas_abinitio_interface.openmx import OpenMXSolver
from abics.applications.latgas_abinitio_interface.params import DFTParams

nsamples = int(sys.argv[2])
ndigits = len(sys.argv[2])
def main_impl(tomlfile):
    dftparams = DFTParams.from_toml(tomlfile)

    if dftparams.solver == 'vasp':
        solver = VASPSolver(dftparams.path)
    elif dftparams.solver == 'qe':
        parallel_level = dftparams.properties.get('parallel_level', {})
        solver = QESolver(dftparams.path, parallel_level=parallel_level)
    elif dftparams.solver == 'aenet':
        solver = aenetSolver(dftparams.path, dftparams.ignore_species, dftparams.solver_run_scheme)
    elif dftparams.solver == 'openmx':
        solver = OpenMXSolver(dftparams.path)
    else:
        print('unknown solver: {}'.format(dftparams.solver))
        sys.exit(1)

    # defect sublattice setup

    configparams = DFTConfigParams.from_toml(tomlfile)
    spinel_config = defect_config(configparams)
    solver_input = solver.input
    solver_input.from_directory(dftparams.base_input_dir[0])

    for i in range(nsamples):
        spinel_config.shuffle()
        solver_input.update_info_by_structure(spinel_config.structure)
        fmtstr = "{:0>"+str(ndigits)+"d}"
        solver_input.write_input(fmtstr.format(i))



def main():
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    main_impl(tomlfile)

if __name__ == "__main__":
    main()
