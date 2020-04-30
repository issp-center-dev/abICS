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
import argparse
from abics.mc_mpi import RXParams, RX_MPI_init
import numpy as np
import scipy.constants as constants


def main():
    # Input parser
    parser = argparse.ArgumentParser(
        description="Reorganize abICS RXMC results by temperature"
    )

    parser.add_argument("inputfi", help="toml input file used for abICS run")

    parser.add_argument(
        "skipsteps",
        nargs="?",
        type=int,
        default=0,
        help="number of thermalization steps to skip in energy averaging." + " Default: 0",
    )

    args = parser.parse_args()
    inputfi = args.inputfi
    nskip = args.skipsteps
    rxparams = RXParams.from_toml(inputfi)
    nreplicas = rxparams.nreplicas
    comm = RX_MPI_init(rxparams)

    myreplica = comm.Get_rank()

    if myreplica == 0:
        if os.path.exists("Tseparate"):
            shutil.move("Tseparate", "Tseparate.bak.{}".format(datetime.datetime.now()))
        os.mkdir("Tseparate")
    comm.Barrier()

    # Separate structure files
    os.mkdir(os.path.join("Tseparate", str(myreplica)))
    Trank_hist = np.load(os.path.join(str(myreplica), "Trank_hist.npy"))
    os.chdir(str(myreplica))
    for j in range(len(Trank_hist)):
        shutil.copy(
            "structure.{}.vasp".format(j),
            os.path.join(os.pardir, "Tseparate", str(Trank_hist[j])),
        )

    # Separate energies
    myreplica_energies = np.load("obs_save.npy")[:, 0]
    for Tid in range(nreplicas):
        T_energies = np.where(Trank_hist == Tid, myreplica_energies, 0)
        T_energies_rcvbuf = np.zeros(T_energies.shape[0], "d")
        comm.Reduce(
            [T_energies, MPI.DOUBLE],
            [T_energies_rcvbuf, MPI.DOUBLE],
            op=MPI.SUM,
            root=Tid,
        )
        if myreplica == Tid:
            np.savetxt(
                os.path.join(os.pardir, "Tseparate", str(Tid), "energies.dat"),
                T_energies_rcvbuf,
            )

    comm.Barrier()

    if myreplica == 0:
        os.chdir(os.path.join(os.pardir, "Tseparate"))
        with open("energies_T.dat", "w") as fi:
            kTs = np.load(os.path.join(os.pardir, "kTs.npy"))
            Ts = kTs / constants.value(u"Boltzmann constant in eV/K")
            for Tid in range(nreplicas):
                energy_mean = np.mean(
                    np.loadtxt(os.path.join(str(Tid), "energies.dat"))[nskip:]
                )
                fi.write("{}\t{}\n".format(Ts[Tid], energy_mean))


if __name__ == "__main__":
    main()
