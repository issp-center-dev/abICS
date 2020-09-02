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

from mpi4py import MPI
import numpy as np
import scipy.constants as constants

from abics.mc import CanonicalMonteCarlo, RandomSampling
from abics.mc_mpi import (
    RX_MPI_init,
    TemperatureRX_MPI,
    RXParams,
    SamplerParams,
    ParallelRandomParams,
    EmbarrassinglyParallelSampling,
)
from abics.applications.latgas_abinitio_interface import default_observer
from abics.applications.latgas_abinitio_interface.model_setup import (
    dft_latgas,
    ObserverParams,
)
from abics.applications.latgas_abinitio_interface.defect import (
    defect_config,
    DFTConfigParams,
)
from abics.applications.latgas_abinitio_interface.run_base_mpi import runner, runner_multistep
from abics.applications.latgas_abinitio_interface.vasp import VASPSolver
from abics.applications.latgas_abinitio_interface.qe import QESolver
from abics.applications.latgas_abinitio_interface.aenet import aenetSolver
from abics.applications.latgas_abinitio_interface.openmx import OpenMXSolver
from abics.applications.latgas_abinitio_interface.params import DFTParams


def main_impl(tomlfile):
    samplerparams = SamplerParams.from_toml(tomlfile)

    rxparams = ParallelRandomParams.from_toml(tomlfile)
    nreplicas = rxparams.nreplicas
    nprocs_per_replica = rxparams.nprocs_per_replica
    comm = RX_MPI_init(rxparams)

    # Set Lreload to True when restarting
    Lreload = rxparams.reload

    nsteps = rxparams.nsteps
    sample_frequency = rxparams.sample_frequency
    print_frequency = rxparams.print_frequency

        

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

    # model setup
    # we first choose a "model" defining how to perform energy calculations and trial steps
    # on the "configuration" defined below
    if len(dftparams.base_input_dir) == 1:
        energy_calculator = runner(
            base_input_dir=dftparams.base_input_dir[0],
            Solver=solver,
            nprocs_per_solver=nprocs_per_replica,
            comm=MPI.COMM_SELF,
            perturb=dftparams.perturb,
            solver_run_scheme=dftparams.solver_run_scheme
        )
    else:
        energy_calculator = runner_multistep(
            base_input_dirs=dftparams.base_input_dir,
            Solver=solver,
            runner=runner,
            nprocs_per_solver=nprocs_per_replica,
            comm=MPI.COMM_SELF,
            perturb=dftparams.perturb,
            solver_run_scheme=dftparams.solver_run_scheme
        )

    myreplica = comm.Get_rank()
    energy_ref = np.loadtxt(os.path.join(str(myreplica), "energies.dat"))
    rootdir = os.getcwd()
    os.makedirs(os.path.join("activelearn",str(myreplica)))
    os.chdir(os.path.join("activelearn",str(myreplica)))
    nskip = 0
    nsteps = 400
    energy_corrlist = []
    relax_max = []
    e_gs =  0 #-.31774118E+04
    from pymatgen import Structure
    from map2perflat import map2perflat
    perf_st = Structure.from_file(os.path.join(rootdir, "MgAl2O4.vasp"))
    perf_st.make_supercell([2,2,2])
    for i in range(nskip, nsteps, 20):
        st_in = Structure.from_file(os.path.join(rootdir, str(myreplica), 'structure.{}.vasp'.format(i)))
        st = map2perflat(perf_st, st_in)
        #st.remove_species(["N", "Li"])
        stbak = st.copy()
        energy, st_rel = energy_calculator.submit(st, os.path.join(rootdir, "activelearn", str(myreplica), "output"))
        st_rel.to("POSCAR", "structure_rel.{}.vasp".format(i))
        energy_corrlist.append([energy_ref[i], energy - e_gs])
        np.savetxt("energy_corr.dat", energy_corrlist)
        dmax = 0
        dmax_id = 0
        for j, site in enumerate(stbak):
            d = site.distance(st_rel[j])
            if d > dmax:
                dmax = d
                dmax_id = j
        relax_max.append([dmax_id, dmax])
        with open("relax_max.dat", "w") as fi:
            for row in relax_max:
                fi.write("{} \t {}\n".format(row[0], row[1]))


if __name__ == '__main__':
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    main_impl(tomlfile)
