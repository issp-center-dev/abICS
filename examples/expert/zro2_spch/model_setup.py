import numpy as np
import random as rand
import sys, os
import copy
import pickle

# from mpi4py import MPI

from pymatgen import Lattice, Structure, Element, PeriodicSite
from pymatgen.io.vasp import Poscar, VaspInput
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
from pymatgen.apps.borg.hive import SimpleVaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from mc import model


def gauss(x, x0, sigma):
    return (
        1.0
        / (np.sqrt(2.0 * np.pi) * sigma)
        * np.exp(-np.power((x - x0) / sigma, 2.0) / 2.0)
    )


class dft_zro2_spch(model):
    """This class defines the DFT ZrO2 space charge model"""

    model_name = "dft_zro2_spch"

    def __init__(
        self,
        calcode,
        vasp_run,
        base_vaspinput,
        matcher,
        matcher_site,
        queen,
        selective_dynamics=None,
    ):
        self.calcode = calcode
        self.matcher = matcher
        self.matcher_site = matcher_site
        self.drone = SimpleVaspToComputedEntryDrone(inc_structure=True)
        self.queen = queen
        self.base_vaspinput = base_vaspinput
        self.vasp_run = vasp_run
        self.selective_dynamics = selective_dynamics

    def energy(self, spch_config):
        """ Calculate total energy of the space charge model"""

        structure = spch_config.structure.get_sorted_structure()

        # if len(spinel_config.calc_history) >= 20:
        #    print("truncate calc_history")
        #    del spinel_config.calc_history[0:10]
        # calc_history = spinel_config.calc_history
        # if calc_history:
        #    # Try to avoid doing dft calculation for the same structure.
        #    # Assuming that calc_history is a list of ComputedStructureEntries
        #    for i in range(len(calc_history)):
        #        if self.matcher.fit(structure, calc_history[i].structure):
        #            print("match found in history")
        #            return calc_history[i].energy
        # print("before poscar")
        if self.selective_dynamics:
            seldyn_arr = [[True, True, True] for i in range(len(structure))]
            for specie in self.selective_dynamics:
                indices = structure.indices_from_symbol(specie)
                for i in indices:
                    seldyn_arr[i] = [False, False, False]
        else:
            seldyn_arr = None

        poscar = Poscar(structure=structure, selective_dynamics=seldyn_arr)
        # print("before vaspinput")
        vaspinput = self.base_vaspinput
        vaspinput.update({"POSCAR": poscar})
        exitcode = self.vasp_run.submit(vaspinput, os.getcwd() + "/output")
        # print("vasp exited with exit code", exitcode)
        if exitcode != 0:
            print("something went wrong")
            sys.exit(1)
        queen = BorgQueen(self.drone)
        # print(os.getcwd())
        queen.serial_assimilate("./output")
        # print(queen.get_data())
        # results = self.queen.get_data()[-1]
        results = queen.get_data()[-1]
        # calc_history.append(results)
        spch_config.structure = results.structure
        # print(results.energy)
        # sys.stdout.flush()

        return np.float64(results.energy)

    # def xparam(self,spinel_config):
    #     '''Calculate number of B atoms in A sites'''

    #     asites = self.matcher_site.get_mapping(spinel_config.structure,
    #                                            spinel_config.Asite_struct)
    #     #print asites
    #     #print spinel_config.structure
    #     #print spinel_config.Asite_struct
    #     x = 0
    #     for i in asites:
    #         if spinel_config.structure.species[i] == Element(spinel_config.Bspecie):
    #             x += 1
    #     x /= float(len(asites))
    #     return x

    def trialstep(self, spch_config, energy_now):

        e0 = energy_now

        structure = spch_config.structure.copy()
        structure_Osite = spch_config.structure.copy()
        structure_Osite.remove_species(["Pt", "Zr"])
        base_Osite = spch_config.Osite_struct.copy()

        # base_Osite.append("He", [0,0,0])
        # base_Osite.append("He", [0.4,0.5,0])
        # base_Osite.append("He", [0.3, 0.2, 0.4])
        # Osite.append("He", [0,0,0])
        # Osite.append("He", [0.4,0.5,0])
        # Osite.append("He", [0.3, 0.2, 0.4])
        # Figure out where the vacancies are by comparing with base_structure
        # filled_sites = self.matcher_site.get_mapping(base_Osite,
        #                                             Osite)
        # if not filled_sites.any():
        #    print("something went wrong at matching")
        #    sys.exit(1)
        # vac_struct = base_Osite
        # vac_struct.remove_sites(filled_sites)
        # Make sure the correct number of vacancies have been detected
        # assert len(vac_struct)==int(spch_config.N_Ovac*len(spch_config.Osite_struct))
        # choose one vacancy and one O atom randomly and flip
        # Osites = structure.indices_from_symbol("O")
        # O_flip = rand.choice(Osites)
        # del structure[O_flip]

        # choose one O atom and move it
        Osites = structure.indices_from_symbol("O")
        O_flip = rand.choice(Osites)
        del structure[O_flip]

        ntrial = 0
        while True:
            ntrial += 1
            Omove_candidate = PeriodicSite(
                "O",
                [
                    rand.random(),
                    rand.random(),
                    rand.uniform(spch_config.ZrO2boxZ[0], spch_config.ZrO2boxZ[1]),
                ],
                structure[0].lattice,
            )
            okflag = True
            for site in structure:
                if site.distance_and_image(Omove_candidate)[0] < 1.5:
                    okflag = False
                    break

            if okflag:
                break

        print("It took ", ntrial, " trials to find space for O")
        print(Omove_candidate.frac_coords)

        structure.append(
            "O", Omove_candidate.frac_coords, properties={"velocities": [0, 0, 0]}
        )
        # print(vac_flip.frac_coords)
        # structure.append("O", vac_candidate.frac_coords,properties={"velocities":[0]})
        # structure.to(fmt="poscar", filename="POSCAR.1")

        # Backup structure of previous step
        structure0 = spch_config.structure
        spch_config.structure = structure

        # do vasp calculation on structure
        e1 = self.energy(spch_config)

        # return old spch structure
        spch_config.structure = structure0

        # Simply pass new structure in dconfig to be used by newconfig():
        dconfig = structure

        dE = e1 - e0

        return dconfig, dE

    def newconfig(self, spch_config, dconfig):
        """Construct the new configuration after the trial step is accepted"""
        spch_config.structure = dconfig
        return spch_config


class spch_config:
    """This class defines the metal/ZrO2 space charge config"""

    def __init__(self, base_structure, N_Ovac, cellsize=1, Vholder="He"):
        self.base_structure = base_structure
        self.base_structure.make_supercell([cellsize, cellsize, 1])
        self.Osite_struct = self.base_structure.copy()
        self.Osite_struct.remove_species(["Pt", "Zr"])
        self.Zrsite_struct = self.base_structure.copy()
        self.Zrsite_struct.remove_species(["O", "Pt"])
        self.Vholder = Vholder
        self.structure = None
        self.calc_history = []
        self.N_Ovac = N_Ovac
        minZ_O = 2.0
        maxZ_O = -1.0
        for site in self.Osite_struct:
            if site.frac_coords[2] < minZ_O:
                minZ_O = site.frac_coords[2]
            if site.frac_coords[2] > maxZ_O:
                maxZ_O = site.frac_coords[2]
        self.ZrO2boxZ = np.array([minZ_O - 0.05, maxZ_O + 0.05])

        sigma = 1.0
        X, Y = np.mgrid[0:23:0.05, 0:1]
        self.base_O_distribution = np.vstack((X.flatten(), Y.flatten())).T
        for site in self.Osite_struct:
            x0 = site.coords[2]
            for grid in self.base_O_distribution:
                grid[1] += gauss(grid[0], x0, sigma)

    def prepare_Ovac(self):
        # Prepare a structure with N_ovac*N_Osites oxygen vacancies
        Osites = self.base_structure.indices_from_symbol("O")
        N_Vsites = int(self.N_Ovac * len(Osites))
        Vacsites = rand.sample(Osites, N_Vsites)
        Vacsites.sort()
        Vacsites.reverse()
        self.structure = self.base_structure.copy()
        for site in Vacsites:
            self.structure.pop(site)

    def prepare_ordered(self):
        self.structure = self.base_structure.copy()

    # def __str__(self):
    #    s = ""
    #    for i in range(self.lenX):
    #        for j in range(self.lenY):
    #            if self.config[i,j] < 0:
    #                s += "-"
    #            else:
    #                s += "+"
    #        s += "\n"
    #    return s


# def writeEandX(calc):
#     with open("energy.out", "a") as f:
#         f.write(str(calc.energy)+"\n")
#         f.flush()
#     with open("xparam.out", "a") as f:
#         xparam = calc.model.xparam(calc.config)
#         f.write(str(xparam)+"\n")
#         f.flush()


def observables(MCcalc, outputfi):
    energy = MCcalc.energy
    # energy2 = energy**2.0
    # xparam = MCcalc.model.xparam(MCcalc.config)
    outputfi.write(
        "\t".join([str(observable) for observable in [MCcalc.kT, energy]]) + "\n"
    )
    outputfi.flush()
    X, Y = np.mgrid[0:23:0.05, 0:1]
    sigma = 1.0
    data = np.vstack((X.flatten(), Y.flatten())).T
    data_base = copy.deepcopy(data)
    #    T_to_rank = pickle.load
    config = copy.deepcopy(MCcalc.config)
    structure = config.structure
    structure.remove_species(["Pt", "Zr"])
    for site in structure:
        x0 = site.coords[2]
        for grid in data:
            grid[1] += gauss(grid[0], x0, sigma)
    if hasattr(config, "base_O_distribution"):
        vac_dist = config.base_O_distribution[:, 1] - data[:, 1]
    else:
        sigma = 1.0
        X, Y = np.mgrid[0:23:0.05, 0:1]
        MCcalc.config.base_O_distribution = np.vstack((X.flatten(), Y.flatten())).T
        for site in MCcalc.config.Osite_struct:
            x0 = site.coords[2]
            for grid in MCcalc.config.base_O_distribution:
                grid[1] += gauss(grid[0], x0, sigma)
        vac_dist = MCcalc.config.base_O_distribution[:, 1] - data[:, 1]
    return vac_dist
