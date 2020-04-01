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
import random as rand
import sys, os
import copy
from pymatgen.io.vasp import Poscar
from pymatgen.apps.borg.hive import SimpleVaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from abics.mc import model


def gauss(x, x0, sigma):
    return (
        1.0
        / (np.sqrt(2.0 * np.pi) * sigma)
        * np.exp(-np.power((x - x0) / sigma, 2.0) / 2.0)
    )


def match_id(lst, obj):
    mapping = []
    for i in range(len(lst)):
        if lst[i] == obj:
            mapping.append(i)
    return mapping


def nomatch_id(lst, obj):
    for i in range(len(lst)):
        if lst[i] != obj:
            mapping.append(i)
    return mapping


class dft_HAp(model):
    """This class defines the DFT HAp space charge model"""

    model_name = "dft_HAp"

    def __init__(
        self,
        calcode,
        vasp_run,
        base_vaspinput,
        matcher_base,  # matcher, matcher_site,
        queen,
        selective_dynamics=None,
        matcher=None,
    ):
        self.calcode = calcode
        self.matcher_base = matcher_base
        self.matcher = matcher
        # self.matcher_site = matcher_site
        self.drone = SimpleVaspToComputedEntryDrone(inc_structure=True)
        self.queen = queen
        self.base_vaspinput = base_vaspinput
        self.vasp_run = vasp_run
        self.selective_dynamics = selective_dynamics

    def update_basestruct(self, HAp_config):
        basesites = self.matcher_base.get_mapping(
            HAp_config.structure, HAp_config.base_structure
        )
        idx = 0
        for i in basesites:
            HAp_config.base_structure[idx] = HAp_config.structure[i]
            idx += 1

    def energy(self, HAp_config, save_history=False):
        """ Calculate total energy of the space charge model"""

        structure = HAp_config.structure.get_sorted_structure()
        if save_history and self.matcher != None:
            calc_history = HAp_config.calc_history
            # Try to avoid doing dft calculation for the same structure.
            # Assuming that calc_history is a list of ComputedStructureEntries
            for i in range(len(calc_history)):
                if self.matcher.fit(structure, calc_history[i].structure0):
                    print("match found in history")
                    return calc_history[i].energy

        if self.selective_dynamics:
            seldyn_arr = [[True, True, True] for i in range(len(structure))]
            for specie in self.selective_dynamics:
                indices = structure.indices_from_symbol(specie)
                for i in indices:
                    seldyn_arr[i] = [False, False, False]
        else:
            seldyn_arr = None

        poscar = Poscar(structure=structure, selective_dynamics=seldyn_arr)
        vaspinput = self.base_vaspinput
        vaspinput.update({"POSCAR": poscar})
        exitcode = self.vasp_run.submit(vaspinput, os.getcwd() + "/output")
        if exitcode != 0:
            print("something went wrong")
            sys.exit(1)
        queen = BorgQueen(self.drone)
        queen.serial_assimilate("./output")
        results = queen.get_data()[-1]

        if save_history:
            results.structure0 = HAp_config.structure
            HAp_config.calc_history.append(results)

        HAp_config.structure = results.structure
        return np.float64(results.energy)

    def trialstep(self, HAp_config, energy_now):

        e0 = energy_now

        # Back up structure and latgas_rep
        structure0 = HAp_config.structure
        latgas_rep0 = HAp_config.latgas_rep.copy()

        if latgas_rep0.count(0) == 0 or rand.random() < 0.5:
            # Flip an OH
            OH_ids = match_id(latgas_rep0, (1, 0)) + match_id(latgas_rep0, (1, 1))
            flip_id = rand.choice(OH_ids)
            if HAp_config.latgas_rep[flip_id] == (1, 0):
                HAp_config.latgas_rep[flip_id] = (1, 1)
            else:
                HAp_config.latgas_rep[flip_id] = (1, 0)

        else:
            # Exchange V with OH or O
            V_ids = match_id(latgas_rep0, 0)
            ex_id = rand.choice(V_ids)
            other_ids = nomatch_id(lagas_rep0, 0)
            ex_other_id = rand.choice(other_ids)
            HAp_config.latgas_rep[ex_id], HAp_config.latgas_rep[ex_other_id] = (
                HAp_config.latgas_rep[ex_other_id],
                HAp_config.latgas_rep[ex_id],
            )

        HAp_config.set_latgas()

        # do vasp calculation on structure
        e1 = self.energy(HAp_config)

        # return old structure
        structure = HAp_config.structure
        HAp_config.structure = structure0
        latgas_rep = HAp_config.latgas_rep
        HAp_config.latgas_rep = latgas_rep0

        # Simply pass new structure and latgas_rep  in dconfig to be used by newconfig():
        dconfig = structure, latgas_rep

        dE = e1 - e0

        return dconfig, dE

    def newconfig(self, HAp_config, dconfig):
        """Construct the new configuration after the trial step is accepted"""
        HAp_config.structure, HAp_config.latgas_rep = dconfig
        self.update_basestruct(HAp_config)
        return HAp_config


class energy_lst_HAp(dft_HAp):
    def __init__(
        self,
        calcode,
        vasp_run,
        base_vaspinput,
        matcher_base,  # matcher, matcher_site,
        queen,
        reps,
        energy_lst,
        selective_dynamics=None,
        matcher=None,
    ):
        super().__init__(
            calcode,
            vasp_run,
            base_vaspinput,
            matcher_base,  # matcher, matcher_site,
            queen,
            selective_dynamics=None,
            matcher=None,
        )
        self.reps = reps
        self.energy_list = energy_lst

    def energy(self, HAp_config, save_history=False):
        rep_id = self.reps.index(tuple(HAp_config.latgas_rep))
        return np.float64(self.energy_list[rep_id])


class HAp_config:
    """This class defines the HAp config with lattice gas mapping"""

    def __init__(self, base_structure, cellsize=[1, 1, 1]):
        self.calc_history = []
        self.cellsize = cellsize
        self.base_structure = base_structure
        site_centers = [[0.0, 0.0, 0.25], [0.0, 0.0, 0.75]]
        self.latgas_rep = [(1, 0)] * np.prod(cellsize) * len(site_centers)
        site_centers = np.array(site_centers)
        site_centers_sc = np.zeros(
            (np.prod(cellsize) * site_centers.shape[0], 3), dtype=float
        )
        idx = 0
        for i in range(cellsize[0]):
            for j in range(cellsize[1]):
                for k in range(cellsize[2]):
                    for l in range(site_centers.shape[0]):
                        site_centers_sc[idx] = site_centers[l] + np.array([i, j, k])
                        idx += 1
        site_centers_sc /= np.array(cellsize)
        self.site_centers = site_centers_sc
        self.base_structure.make_supercell([cellsize[0], cellsize[1], cellsize[2]])
        self.supercell = self.base_structure.lattice.matrix
        group_dict = {
            0: [],
            (1, 0): [
                ["O", np.array([0.0, 0.0, -0.915 / 2])],
                ["H", np.array([0.0, 0.0, 0.915 / 2])],
            ],
            (1, 1): [
                ["O", np.array([0.0, 0.0, 0.915 / 2])],
                ["H", np.array([0.0, 0.0, -0.915 / 2])],
            ],
            2: [["O", np.array([0.0, 0.0, 0.0])]],
        }
        invSuper = np.linalg.inv(self.supercell)
        for key in group_dict.keys():
            for atom in group_dict[key]:
                atom[1] = np.dot(atom[1], invSuper)
        self.group_dict = group_dict

    def set_latgas(self, latgas_rep=False):
        if latgas_rep:
            self.latgas_rep = latgas_rep
        numsites = len(self.latgas_rep)
        assert numsites == self.site_centers.shape[0]
        self.structure = copy.deepcopy(self.base_structure)
        for isite in range(numsites):
            gid = self.latgas_rep[isite]
            for atom in self.group_dict[gid]:
                self.structure.append(
                    atom[0],
                    atom[1] + self.site_centers[isite],
                    properties={"velocities": [0, 0, 0]},
                )


def observables(MCcalc, outputfi):
    energy = MCcalc.energy
    nup = MCcalc.config.latgas_rep.count((1, 0))
    ndown = MCcalc.config.latgas_rep.count((1, 1))
    tot_pol = nup - ndown
    # energy2 = energy**2.0
    # xparam = MCcalc.model.xparam(MCcalc.config)
    outputfi.write(
        "\t".join(
            [
                str(observable)
                for observable in [MCcalc.kT, energy, nup, ndown, tot_pol, tot_pol ** 2]
            ]
        )
        + "\n"
    )
    outputfi.flush()
    return [MCcalc.kT, energy, nup, ndown, tot_pol, tot_pol ** 2]
