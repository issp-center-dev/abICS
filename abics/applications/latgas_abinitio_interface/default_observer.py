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

from __future__ import annotations

import os
import sys
import copy

import numpy as np

from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator

from abics import __version__
from abics.exception import InputError
from abics.util import expand_path
from abics.mc import ObserverBase, MCAlgorithm
from abics.applications.latgas_abinitio_interface.base_solver import create_solver
from abics.applications.latgas_abinitio_interface.run_base_mpi import (
    Runner,
)
from abics.applications.latgas_abinitio_interface.params import DFTParams


def similarity(
    structure: Structure,
    reference: Structure,
    sp: str,
    matcher: StructureMatcher,
):
    """
    Arguments
    =========

    structure: Structure
    references: dict[str, Structure]
        Structure for each element (key) in reference structure
    matcher: StructureMatcher

    """
    matched = 0
    nsites = 0
    sites = matcher.get_mapping(structure, reference)
    if sites is None:
        raise RuntimeError("structure does not match with reference structure")
    species = [str(sp) for sp in structure.species]
    for i in sites:
        if species[i] == sp:
            matched += 1
    nsites += len(sites)
    return matched / nsites


class DefaultObserver(ObserverBase):
    """
    Default observer.

    Attributes
    ----------
    minE : float
        Minimum of energy
    """

    references: dict[str, Structure] | None
    reference_species: list[str]
    calculators: list

    def __init__(self, comm, Lreload=False, params={}):
        """

        Parameters
        ----------
        comm: mpi4py.MPI.Intracomm
            MPI communicator
        Lreload: bool
            Reload or not
        """
        super().__init__(comm, Lreload, params)
        self.minE = 100000.0
        myrank = comm.Get_rank()
        if Lreload:
            with open(os.path.join(str(myrank), "minEfi.dat"), "r") as f:
                self.minE = float(f.readlines()[-1])
            with open(os.path.join(str(myrank), "obs.dat"), "r") as f:
                self.lprintcount = int(f.readlines()[-1].split()[0]) + 1

        self.names = ["energy"]

        params_solvers = params.get("solver", [])
        self.calculators = []
        for params_solver in params_solvers:
            dft_params = DFTParams.from_dict(params_solver)
            solver = create_solver(params_solver["type"], dft_params)
            obsname = params_solver["name"]
            if obsname in self.names:
                raise InputError(f"Duplicated observer name: {obsname}")
            nrunners = len(dft_params.base_input_dir)
            perturbs = [0.0] * nrunners
            perturbs[0] = dft_params.perturb
            self.names.append(obsname)
            runners = [
                Runner(
                    base_input_dir=bid,
                    Solver=solver,
                    nprocs_per_solver=1,
                    comm=comm,
                    perturb=perturb,
                    nthreads_per_proc=1,
                    solver_run_scheme=dft_params.solver_run_scheme,
                    use_tmpdir=dft_params.use_tmpdir,
                )
                for perturb, bid in zip(perturbs, dft_params.base_input_dir)
            ]
            self.calculators.append(
                {"name": obsname, "runners": runners}
            )

        if "similarity" in params:
            params_similarity = params["similarity"]
            reference = Structure.from_file(params_similarity["reference_structure"])
            ignored_species = params_similarity.get("ignored_species", [])
            if isinstance(ignored_species, str):
                ignored_species = [ignored_species]
            reference.remove_species(ignored_species)
            sp_set = {str(sp) for sp in reference.species}
            self.references = {}
            self.reference_species = []
            for sp in sp_set:
                bs = reference.copy()
                bs.remove_species(sp_set - {sp})
                self.reference_species.append(sp)
                self.references[sp] = bs
                self.names.append(f"similarity_{sp}")
            self.matcher = StructureMatcher(
                ltol=0.1,
                primitive_cell=False,
                allow_subset=True,
                comparator=FrameworkComparator(),
                ignored_species=ignored_species,
            )
        else:
            self.references = None
            self.reference_species = []

    def logfunc(self, calc_state: MCAlgorithm) -> tuple[float, ...]:
        assert calc_state.config is not None
        structure: Structure = calc_state.config.structure_norel
        energy = calc_state.energy
        result = [energy]
        if energy < self.minE:
            self.minE = energy
            with open("minEfi.dat", "a") as f:
                f.write(str(self.minE) + "\n")
            structure.to(fmt="POSCAR", filename="minE.vasp")
        for calculator in self.calculators:
            name = calculator["name"]
            runners: list[Runner] = calculator["runners"]
            new_st = structure
            value = 0.0
            for i, runner in enumerate(runners):
                output_dir = os.path.join(os.getcwd(), f"output{i}")
                value, new_st = runner.submit(new_st, output_dir)
            result.append(value)

        if self.references is not None:
            assert(self.reference_species is not None)
            for sp in self.reference_species:
                ref = self.references[sp]
                sim = similarity(structure, ref, sp, self.matcher)
                result.append(sim)
        return tuple(result)

    def writefile(self, calc_state: MCAlgorithm) -> None:
        assert calc_state.config is not None
        calc_state.config.structure.to(
            fmt="POSCAR", filename="structure." + str(self.lprintcount) + ".vasp"
        )
        calc_state.config.structure_norel.to(
            fmt="POSCAR", filename="structure_norel." + str(self.lprintcount) + ".vasp"
        )


class EnsembleErrorObserver(DefaultObserver):
    def __init__(self, comm, energy_calculators, Lreload=False):
        """

        Parameters
        ----------
        comm: mpi4py.MPI.Intracomm
            MPI communicator
        energy_calculators: abics.applications.latgas_abinitio_interface.run_base_mpi.runner
            setup of the energy calculator
        Lreload: bool
            Reload or not
        """
        super().__init__(comm, Lreload)
        self.calculators = energy_calculators
        self.comm = comm

    def logfunc(self, calc_state: MCAlgorithm):
        if calc_state.energy < self.minE:
            self.minE = calc_state.energy
            with open("minEfi.dat", "a") as f:
                f.write(str(self.minE) + "\n")
            calc_state.config.structure.to(fmt="POSCAR", filename="minE.vasp")
        energies = [calc_state.energy]
        npar = self.comm.Get_size()
        if npar > 1:
            assert npar == len(self.calculators)
            myrank = self.comm.Get_rank()
            energy, _ = self.calculators[myrank].submit(
                calc_state.config.structure,
                os.path.join(os.getcwd(), "ensemble{}".format(myrank)),
            )
            energies_tmp = self.comm.allgather(energy)
            std = np.std(energies_tmp, ddof=1)

        else:
            energies_tmp = []
            for i, calculator in enumerate(self.calculators):
                energy, _ = calculator.submit(
                    calc_state.config.structure,
                    os.path.join(os.getcwd(), "ensemble{}".format(i)),
                )
                energies_tmp.append(energy)
            std = np.std(energies_tmp, ddof=1)
        energies.extend(energies_tmp)
        energies.append(std)
        return np.asarray(energies)


class EnsembleParams:
    @classmethod
    def from_dict(cls, d):
        """
        Read information from dictionary

        Parameters
        ----------
        d: dict
            Dictionary

        Returns
        -------
        params: EnsembleParams object
            self
        """

        params = cls()
        base_input_dirs = d.get("base_input_dirs", ["./baseinput"])
        if isinstance(base_input_dirs, str):
            base_input_dirs = [base_input_dirs]
        params.base_input_dirs = base_input_dirs = list(
            map(lambda x: expand_path(x, os.getcwd()), base_input_dirs)
        )
        params.solver = d["type"]
        params.path = expand_path(d["path"], os.getcwd())
        params.perturb = d.get("perturb", 0.1)
        params.solver_run_scheme = d.get("run_scheme", "mpi_spawn_ready")
        params.ignore_species = d.get("ignore_species", None)
        params.constraint_module = d.get("constraint_module", None)
        params.properties = d

        return params

    @classmethod
    def from_toml(cls, f):
        """
        Read information from toml file

        Parameters
        ----------
        f: str
            Name of input toml File

        Returns
        -------
        oDFTParams: DFTParams object
            self
        """
        import toml

        d = toml.load(f)
        return cls.from_dict(d["ensemble"])


# For backward compatibility
if __version__ < "3":
    default_observer = DefaultObserver
    ensemble_error_observer = EnsembleErrorObserver
