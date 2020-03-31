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

import numpy as np
import pickle


if __name__ == "__main__":
    energy_lst = pickle.load(open("energy_reps.pickle", "rb"))
    reps = pickle.load(open("latgas_reps.pickle", "rb"))
    energy_lst = energy_lst - np.sum(energy_lst) / len(energy_lst)  # energy_lst[0]
    kb = 8.617e-5
    kTs = np.array([400.0 + 10.0 * i for i in range(150)]) * kb
    tot_pol_expect = np.zeros(len(kTs))
    num_V_O_expect = np.zeros(len(kTs))
    num_O_HO_expect = np.zeros(len(kTs))
    num_O_OH_expect = np.zeros(len(kTs))
    Zs = np.zeros(len(kTs))

    # energy_minfind = energy_lst.copy()
    # for i in range(20):
    #    id = np.argmin(energy_minfind)
    #    print(reps[id], energy_minfind[id])
    #    energy_minfind[id] = 100

    # Polarization
    tot_pol2_rep = (
        np.array([np.abs(rep.count((1, 0)) - rep.count((1, 1))) for rep in reps]) / 4
    )
    # print(tot_pol2_rep)

    # Nearest neighbor occupation
    num_V_O = []
    num_O_HO = []
    num_O_OH = []
    for rep in reps:
        V_id = rep.index(0)
        O_id = rep.index(2)

        if O_id == len(rep) - 1:
            O_id_p1 = 0
        else:
            O_id_p1 = O_id + 1
        if rep[O_id_p1] == 0 or rep[O_id - 1] == 0:
            num_V_O.append(1)
        else:
            num_V_O.append(0)

        num_O_HO_tmp = 0
        if rep[O_id_p1] == (1, 1):
            num_O_HO_tmp += 1
        if rep[O_id - 1] == (1, 0):
            num_O_HO_tmp += 1
        num_O_HO.append(num_O_HO_tmp)

        num_O_OH_tmp = 0
        if rep[O_id_p1] == (1, 0):
            num_O_OH_tmp += 1
        if rep[O_id - 1] == (1, 1):
            num_O_OH_tmp += 1
        num_O_OH.append(num_O_OH_tmp)

    num_V_O = np.array(num_V_O)
    num_O_HO = np.array(num_O_HO)
    num_O_OH = np.array(num_O_OH)
    # print(num_V_O+num_O_HO+num_O_OH)
    # print(reps[1])
    # print(num_V_O[1])
    # print(kTs)
    # print(tot_pol2_rep)
    # print((-energy_lst/kTs[0]))
    # calculate partition function
    for i in range(len(kTs)):
        Zs[i] = np.sum(np.exp(-energy_lst / kTs[i]))
        tot_pol_expect[i] = np.sum(tot_pol2_rep * np.exp(-energy_lst / kTs[i])) / Zs[i]
        num_V_O_expect[i] = np.sum(num_V_O * np.exp(-energy_lst / kTs[i])) / Zs[i]
        num_O_HO_expect[i] = np.sum(num_O_HO * np.exp(-energy_lst / kTs[i])) / Zs[i]
        num_O_OH_expect[i] = np.sum(num_O_OH * np.exp(-energy_lst / kTs[i])) / Zs[i]

        print(
            kTs[i] / kb,
            tot_pol_expect[i],
            num_V_O_expect[i],
            num_O_HO_expect[i],
            num_O_OH_expect[i],
            num_V_O_expect[i] + num_O_HO_expect[i] + num_O_OH_expect[i],
        )
