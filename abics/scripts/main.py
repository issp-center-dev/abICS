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

from typing import MutableMapping

import sys
import datetime

from mpi4py import MPI
import toml

from abics import __version__


def main_impl(params: MutableMapping):
    solver_type = params["sampling"]["solver"]["type"]
    if solver_type == "potts":
        import abics.scripts.main_potts
        abics.scripts.main_potts.main_potts(params)
    else:
        import abics.scripts.main_dft_latgas
        abics.scripts.main_dft_latgas.main_dft_latgas(params)


def main():
    now = datetime.datetime.now()
    
    if MPI.COMM_WORLD.Get_rank() == 0:
        print(f"Running abics_sampling (abICS v{__version__}) on {now}")
        
    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    if MPI.COMM_WORLD.Get_rank() == 0:
        print("-Reading input from: {}".format(tomlfile))
    params = toml.load(tomlfile)
    main_impl(params)


if __name__ == "__main__":
    main()
