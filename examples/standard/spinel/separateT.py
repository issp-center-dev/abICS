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

from mpi4py import MPI
import shutil
import os, sys
import datetime
from pymatgen import Structure
from abics.mc_mpi import RXParams
import numpy as np

tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
rxparams = RXParams.from_toml(tomlfile)
nreplicas = rxparams.nreplicas

comm = MPI.COMM_WORLD
if comm.Get_size() < nreplicas:
    if comm.Get_rank() == 0:
        print("please run with at least as many processes as the number of replicas")
    sys.exit()

myreplica = comm.Get_rank()

if myreplica == 0:
    if os.path.exists("Tseparate"):
        shutil.move("Tseparate", "Tseparate.bak.{}".format(datetime.datetime.now()))
    os.mkdir("Tseparate")
comm.Barrier()

if myreplica < nreplicas:
    os.mkdir("Tseparate/{}".format(myreplica))
    Trank_hist = np.load(open(str(myreplica) + "/Trank_hist.npy", "rb"))
    os.chdir(str(myreplica))
    for j in range(len(Trank_hist)):
        shutil.copy(
            "structure." + str(j) + ".vasp", "../Tseparate/" + str(Trank_hist[j])
        )

