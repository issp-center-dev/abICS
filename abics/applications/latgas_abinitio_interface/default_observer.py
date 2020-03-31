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

from abics.mc import observer_base
import os


class default_observer(observer_base):
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
        super(default_observer, self).__init__()
        self.minE = 100000.0
        myrank = comm.Get_rank()
        if Lreload:
            with open(os.path.join(str(myrank), "minEfi.dat"), "r") as f:
                self.minE = float(f.readlines()[-1])
            with open(os.path.join(str(myrank), "obs.dat"), "r") as f:
                self.lprintcount = int(f.readlines()[-1].split()[0]) + 1

    def logfunc(self, calc_state):
        """

        Parameters
        ----------
        calc_state: MCalgo
        Object of Monte Carlo algorithm
        Returns
        -------
        calc_state.energy : float
        Minimum energy
        """
        if calc_state.energy < self.minE:
            self.minE = calc_state.energy
            with open("minEfi.dat", "a") as f:
                f.write(str(self.minE) + "\n")
            calc_state.config.structure.to(fmt="POSCAR", filename="minE.vasp")
        return calc_state.energy

    def writefile(self, calc_state):
        """

        Parameters
        ----------
        calc_state: MCalgo
        Object of Monte Carlo algorithm

        """
        calc_state.config.structure.to(
            fmt="POSCAR", filename="structure." + str(self.lprintcount) + ".vasp"
        )
