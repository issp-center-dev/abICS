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

from typing import Tuple

import os
import numpy as np

from abics import __version__
from abics.util import expand_path
from abics.mc import ObserverBase, MCAlgorithm


class DefaultObserver(ObserverBase):
    """
    Default observer.

    Attributes
    ----------
    minE : float
        Minimum of energy
    """

    def __init__(self, comm, Lreload=False):
        """

        Parameters
        ----------
        comm: mpi4py.MPI.Intracomm
            MPI communicator
        Lreload: bool
            Reload or not
        """
        super(DefaultObserver, self).__init__()
        self.minE = 100000.0
        myrank = comm.Get_rank()
        if Lreload:
            with open(os.path.join(str(myrank), "minEfi.dat"), "r") as f:
                self.minE = float(f.readlines()[-1])
            with open(os.path.join(str(myrank), "obs.dat"), "r") as f:
                self.lprintcount = int(f.readlines()[-1].split()[0]) + 1

    def logfunc(self, calc_state: MCAlgorithm) -> Tuple[float]:
        if calc_state.energy < self.minE:
            self.minE = calc_state.energy
            with open("minEfi.dat", "a") as f:
                f.write(str(self.minE) + "\n")
            calc_state.config.structure_norel.to(fmt="POSCAR", filename="minE.vasp")
        return (calc_state.energy,)

    def writefile(self, calc_state: MCAlgorithm) -> None:
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
        super(EnsembleErrorObserver, self).__init__(comm, Lreload)
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
