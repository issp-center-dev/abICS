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

from itertools import product
import numpy as np
import random as rand
import sys
import os
import copy

# from mpi4py import MPI

from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
import pymatgen.analysis.structure_analyzer as analy

from abics.exception import InputError
from abics.mc import model
from abics.util import read_vector, read_matrix, read_tensor


def gauss(x, x0, sigma):
    """
    Gaussian function


    Parameters
    ----------
    x: float
        The position
    x0: float
        The position of the center of the peak
    sigma: float
        The standard deviation
    Returns
    -------
        value: float
    """
    return (
        1.0
        / (np.sqrt(2.0 * np.pi) * sigma)
        * np.exp(-np.power((x - x0) / sigma, 2.0) / 2.0)
    )


def match_id(lst, obj):
    """
    Return the index list of lst which matches obj.

    Parameters
    ----------
    lst: list
    obj: object

    Returns
    -------
        The index list of lst which matches obj.
    """
    mapping = []
    for i in range(len(lst)):
        if lst[i] == obj:
            mapping.append(i)
    return mapping


def nomatch_id(lst, obj):
    """
    Return the index list of lst which does not match obj.

    Parameters
    ----------
    lst: list
    obj: object

    Returns
    -------
        The index list of lst which does not match obj.
    """
    mapping = []
    for i in range(len(lst)):
        if lst[i] != obj:
            mapping.append(i)
    return mapping


def match_latgas_group(latgas_rep, group):
    """
    Return the index list of latgas_rep which matches group.name.

    Parameters
    ----------
    latgas_rep: list
    group: object
        Notes: the group must have name
    Returns
    -------
        The index list of latgas_rep which matches group.name.
    """
    mapping = []
    for i in range(len(latgas_rep)):
        if latgas_rep[i][0] == group.name:
            mapping.append(i)
    return mapping


class dft_latgas(model):
    """
    This class defines the DFT lattice gas mapping  model
    """

    model_name = "dft_latgas"

    def __init__(
        self,
        abinitio_run,
        save_history=True,
        l_update_basestruct=False,
        check_ion_move=False,
        ion_move_tol=0.7,
    ):
        """

        Parameters
        ----------
        abinitio_run: runner object
            Runner (manager) of external solver program
        save_history: boolean
        l_update_basestruct: boolean
        check_ion_move: boolean
        ion_move_tol: float
        """
        self.matcher = StructureMatcher(primitive_cell=False, allow_subset=False)
        self.abinitio_run = abinitio_run
        self.save_history = save_history
        self.l_update_basestruct = l_update_basestruct
        self.check_ion_move = check_ion_move
        self.ion_move_tol = ion_move_tol

    def energy(self, config):
        """
        Calculate total energy

        Parameters
        ----------
        config: config object
            Configurations

        Returns
        -------
        energy: float
        """

        config.structure.sort(key=lambda site: site.species_string)
        structure = config.structure
        if self.save_history:
            calc_history = config.calc_history
            # Try to avoid doing dft calculation for the same structure.
            # Assuming that calc_history is a list of ComputedStructureEntries
            for i in range(len(calc_history)):
                if self.matcher.fit(structure, calc_history[i][1]):
                    print("match found in history")
                    sys.stdout.flush()
                    config.structure = calc_history[i][2]
                    return calc_history[i][0]

        structure0 = structure
        energy, structure = self.abinitio_run.submit(
            structure, os.path.join(os.getcwd(), "output")
        )
        if self.check_ion_move:
            relax_analy = analy.RelaxationAnalyzer(structure0, structure)
            data = relax_analy.get_percentage_bond_dist_changes()
            breakflag = False
            for ion in self.check_ion_move:
                if breakflag:
                    break
                for index in structure0.indices_from_symbol(ion):
                    if breakflag:
                        break
                    for i in data[index].keys():
                        if data[index][i] > self.ion_move_tol:
                            energy = float("inf")
                            print("ion relaxed out of initial site")
                            sys.stdout.flush()
                            breakflag = True
                            break
        if self.save_history:
            config.calc_history.append((energy, structure0.copy(), structure.copy()))
            if len(config.calc_history) == 25:
                del config.calc_history[0:5]

        config.structure = structure
        return np.float64(energy)

    def trialstep(self, config, energy_now):
        """

        Parameters
        ----------
        config: config object
            Configurations
        energy_now: float
            Present energy

        Returns
        -------
        dconfig: float
            Difference of configurations
        dE : float
            Difference of energies
        """

        e0 = energy_now

        # Back up structure and defect_sublattices
        structure0 = copy.deepcopy(config.structure)
        defect_sublattices0 = copy.deepcopy(config.defect_sublattices)

        # for defect_sublattice in [rand.choice(config.defect_sublattices)]: #config.defect_sublattices:

        defect_sublattice = rand.choice(config.defect_sublattices)
        latgas_rep = defect_sublattice.latgas_rep
        # print(latgas_rep)

        # If there is more than one group on the defect_sublattice,
        # we either change orientation of one group, or  exchange groups between sites
        if len(defect_sublattice.groups) > 1:
            random_divide = 0.5
        else:
            random_divide = 2.0
        if defect_sublattice.groups_orr and rand.random() < random_divide:
            # Change orientation of one group with orientation attributes
            # Let's first locate sites that have such groups
            rot_ids = []
            for group in defect_sublattice.groups_orr:
                rot_ids += match_latgas_group(latgas_rep, group)
            # Choose the site to rotate and rotate
            rot_id = rand.choice(rot_ids)
            rot_group = defect_sublattice.group_dict[latgas_rep[rot_id][0]]
            rot_group_orr = latgas_rep[rot_id][1]
            # Remove the current orientation from the list of possible orientations
            or_choices = [
                orx for orx in range(rot_group.orientations) if orx != rot_group_orr
            ]
            latgas_rep[rot_id] = [rot_group.name, rand.choice(or_choices)]

        else:
            # Exchange different groups between sites
            ex1_group, ex2_group = rand.sample(defect_sublattice.groups, 2)
            ex1_id = rand.choice(match_latgas_group(latgas_rep, ex1_group))
            ex2_id = rand.choice(match_latgas_group(latgas_rep, ex2_group))
            latgas_rep[ex1_id], latgas_rep[ex2_id] = (
                latgas_rep[ex2_id],
                latgas_rep[ex1_id],
            )
        # print(latgas_rep)
        config.set_latgas()

        # do vasp calculation on structure
        e1 = self.energy(config)

        # return old structure
        structure = config.structure
        config.structure = structure0
        defect_sublattices = config.defect_sublattices
        config.defect_sublattices = defect_sublattices0

        # Simply pass new structure and latgas_rep  in dconfig to be used by newconfig():
        dconfig = structure, defect_sublattices
        if e0 == float("inf"):
            dE = e1
        else:
            dE = e1 - e0

        return dconfig, dE

    def newconfig(self, config, dconfig):
        """
        Update config by the trial step, dconfig

        Parameters
        ----------
        config: config object
            Configuration
        dconfig: config object
            Difference of configuration

        Returns
        -------
        config: config object
            New configuration
        """
        config.structure, config.defect_sublattices = dconfig
        if self.l_update_basestruct:
            self.update_basestruct(config)
        return config

class energy_lst(dft_latgas):
    def __init__(
        self,
        calcode,
        vasp_run,
        base_vaspinput,
        matcher_base,  # matcher, matcher_site,
        queen,
        reps,
        energy_lst,
        matcher=None
    ):
        """

        Parameters
        ----------
        calcode:
        vasp_run: runner object
            Runner (manager) of external solver program
        base_vaspinput:
        matcher_base:
        queen:
        reps:
        energy_lst: list
            Energy list
        matcher:
        """
        super().__init__(
            calcode,
            vasp_run,
            base_vaspinput,
            matcher_base,  # matcher, matcher_site,
            queen,
            matcher=None,
        )
        self.reps = reps
        self.energy_list = energy_lst

    def energy(self, config, save_history=False):
        """

        Parameters
        ----------
        config: config object
            Configuration
        save_history: boolean

        Returns
        -------
        energy: float

        """
        rep_id = self.reps.index(tuple(config.latgas_rep))
        return np.float64(self.energy_list[rep_id])


class group(object):
    def __init__(self, name, species, *,
            coords=None,
            relaxations=None,
            magnetizations=None):
        """

        Parameters
        ----------
        name: str
            The name of atomic group
        species: str
            The atomic species belonging to the atom group
        coords: numpy array
            The coordinates of each atom in the atom group.
        relaxations: numpy array
            Whether to perform structure optimization or not
        magnetization: numpy array
            Magnetizations (inbalance of up/down spins)
        """
        self.name = name
        self.species = species
        self.coords = np.array(coords) if coords is not None else np.zeros((1, 3))
        self.relaxations = np.array(relaxations) if relaxations is not None else np.ones((1, 3), dtype=bool)
        self.magnetizations = np.array(magnetizations) if magnetizations is not None else np.zeros(1)
        self.orientations = len(coords)
        if self.orientations == 0:
            self.orientations = 1
        self.natoms = len(species)


class defect_sublattice(object):
    def __init__(self, site_centers, groups):
        """

        Parameters
        ----------
        site_centers: list
            Center coordinates at each groups
        groups: list
            List of groups
        """
        self.site_centers = np.array(site_centers)
        self.groups = groups
        self.groups_orr = []
        self.group_dict = {}
        for group in groups:
            if group.orientations > 1:
                self.groups_orr.append(group)
            self.group_dict[group.name] = group

    @classmethod
    def from_dict(cls, d):
        """

        Parameters
        ----------
        d: dict

        Returns
        -------
        site_centers: list
            Center coordinates at each groups
        groups: list
            List of groups
        """
        site_centers = read_matrix(d["coords"])
        groups = []
        for g in d["groups"]:
            name = g["name"]
            species = g.get("species", [name])
            n = len(species)

            coords = read_tensor(g.get("coords", [[[0,0,0]]]), rank=3)
            m = coords.shape[1]
            if m != n:
                raise InputError(
                    'number of atoms mismatch in group [{}]: "species"={}, "coords"={}'.format(
                        name, n, m
                    )
                )

            relaxation = read_matrix(g.get("relaxation", np.ones((n, 3))), dtype=bool)
            m = relaxation.shape[0]
            if m != n:
                raise InputError(
                    'number of atoms mismatch in group [{}]: "species"={}, "relaxation"={}'.format(
                        name, n, m
                    )
                )

            mag = read_vector(g.get("magnetization", np.zeros(n)))
            m = len(mag)
            if m != n:
                raise InputError(
                    'number of atoms mismatch in group [{}]: "species"={}, "magnetization"={}'.format(
                        name, n, m
                    )
                )
            groups.append(group(name, species, coords=coords, relaxations=relaxation, magnetizations=mag))
        return cls(site_centers, groups)


def base_structure(lat, dict_str):
    """

    Parameters
    ----------
    lat: pymatgen.Lattice
    dict_str: dict
        Dictionary of base structure

    Returns
    -------
    st: pymatgen.Structure
    """
    if len(dict_str) == 1 and not dict_str[0]:
        return Structure(
            lattice=lat,
            species=[],
            coords=[],
            site_properties={"seldyn": np.zeros((0, 3), dtype=bool),
                             "magnetization": np.zeros(0)}
        )
    elems = []
    coords = []
    relaxations = []
    magnetizations = []
    for tc in dict_str:
        if "type" not in tc:
            raise InputError('"type" is not found in "base_structure"')
        sp = tc["type"]
        crds = read_matrix(tc["coords"])
        n = crds.shape[0]

        if "relaxation" in tc:
            relax = read_matrix(tc["relaxation"], dtype=bool)
            m = relax.shape[0]
            if m != n:
                raise InputError(
                    'number of base atoms mismatch: "coords"={}, "relaxation"={}'.format(
                        n, m
                    )
                )
        else:
            relax = np.ones((n, 3), dtype=bool)  # all True

        if "magnetization" in tc:
            mag = tc["magnetization"]
            if not isinstance(mag, list):
                raise InputError('"magnetization" should be a list of floats')
            try:
                mag = [float(x) for x in mag]
            except ValueError:
                raise InputError('"magnetization" should be a list of floats')
            m = len(mag)
            if m != n:
                raise InputError(
                    'number of base atoms mismatch: "coords"={}, "magnetization"={}'.format(
                        n, m
                    )
                )
        else:
            mag = np.zeros(n)

        elems.append([sp] * n)
        coords.append(crds)
        relaxations.append(relax)
        magnetizations.append(mag)
    elems = sum(elems, [])
    coords = np.concatenate(coords, axis=0)
    relaxations = np.concatenate(relaxations, axis=0)
    magnetizations = np.concatenate(magnetizations, axis=0)

    return Structure(
        lattice=lat,
        species=elems,
        coords=coords,
        site_properties={"seldyn": relaxations, "magnetization": magnetizations},
    )


class config:
    """This class defines the config with lattice gas mapping"""

    def __init__(
        self,
        base_structure,
        defect_sublattices,
        num_defects,
        cellsize=[1, 1, 1],
        perfect_structure=None,
    ):
        """

        Parameters
        ----------
        base_structure : pymatgen.Structure
            Structure of base sites (unsampled sites)
        defect_sublattices : defect_sublattice
            Structure of defects (sampled sites)
        num_defects : dict
            {group name: number of defects}
        cellsize : list, optional
            Cell size, by default [1, 1, 1]
        perfect_structure : pymatgen.Structure, optional
            Strucure of all sites (union of base and defect), by default None
        """
        try:
            num_defect_sublat = len(defect_sublattices)
        except TypeError:
            num_defect_sublat = 1
            defect_sublattices = [defect_sublattices]
            num_defects = [num_defects]
        self.matcher_base = StructureMatcher(primitive_cell=False, allow_subset=True)
        self.matcher_frame = StructureMatcher(
            stol=0.4,
            primitive_cell=False,
            allow_subset=True,
            comparator=FrameworkComparator(),
        )

        self.calc_history = []
        self.cellsize = cellsize
        self.base_structure = base_structure
        if self.base_structure.num_sites == 0:
            # we need at least one site for make_supercell
            self.base_structure.append("H", np.array([0, 0, 0]))
            self.base_structure.make_supercell([cellsize[0], cellsize[1], cellsize[2]])
            self.base_structure.remove_sites(range(self.base_structure.num_sites))
        else:
            self.base_structure.make_supercell([cellsize[0], cellsize[1], cellsize[2]])
        if perfect_structure:
            self.perfect_structure = perfect_structure
            self.perfect_structure.make_supercell(
                [cellsize[0], cellsize[1], cellsize[2]]
            )
        self.supercell = self.base_structure.lattice.matrix
        self.n_sublat = num_defect_sublat
        invSuper = np.linalg.inv(self.supercell)

        # add supercell information to defect_sites
        num_sites = 0
        ntot_defects = 0
        sublat_id = 0
        for defect_sublattice in defect_sublattices:
            site_centers = defect_sublattice.site_centers
            defect_sublattice.site_centers_sc = np.zeros(
                (np.prod(cellsize) * site_centers.shape[0], 3), dtype=float
            )
            idx = 0
            for (idx, (i, j, k, l)) in enumerate(product(range(cellsize[0]),
                                                         range(cellsize[1]),
                                                         range(cellsize[2]),
                                                         range(site_centers.shape[0]))):
                defect_sublattice.site_centers_sc[idx] = site_centers[l] + np.array([i, j, k])
                idx += 1
            defect_sublattice.site_centers_sc /= np.array(cellsize)
            num_sites += len(defect_sublattice.site_centers_sc)

            # change cartesian coordinates to fractional
            for group in defect_sublattice.groups:
                for i in range(group.orientations):
                    for j in range(group.natoms):
                        group.coords[i][j] = np.dot(group.coords[i][j], invSuper)

            # fill the lattice gas representation list
            latgas_rep = []
            for group in num_defects[sublat_id].keys():
                latgas_rep.extend(
                    [[group, 0] for i in range(num_defects[sublat_id][group])]
                )
                ntot_defects += num_defects[sublat_id][group]
            defect_sublattice.latgas_rep = latgas_rep

            sublat_id += 1

        self.defect_sublattices = defect_sublattices

        assert num_sites == ntot_defects
        self.set_latgas()

    def set_latgas(self, defect_sublattices=False):
        """

        Parameters
        ----------
        defect_sublattices: pymatgen.Structure
        """
        if defect_sublattices:
            self.defect_sublattices = defect_sublattices
        assert len(self.defect_sublattices) == self.n_sublat

        self.structure = copy.deepcopy(self.base_structure)
        for defect_sublattice in self.defect_sublattices:
            latgas_rep = defect_sublattice.latgas_rep
            assert len(latgas_rep) == len(defect_sublattice.site_centers_sc)
            for isite in range(len(latgas_rep)):
                grp_name = latgas_rep[isite][0]
                orr = latgas_rep[isite][1]
                group = defect_sublattice.group_dict[grp_name]
                for j in range(group.natoms):
                    self.structure.append(
                        group.species[j],
                        group.coords[orr][j] + defect_sublattice.site_centers_sc[isite],
                        properties={"seldyn": group.relaxations[j, :],
                                    "magnetization": group.magnetizations[j]},
                    )

    def shuffle(self):
        for defect_sublattice in self.defect_sublattices:
            latgas_rep = defect_sublattice.latgas_rep
            rand.shuffle(latgas_rep)
            for site in latgas_rep:
                group = defect_sublattice.group_dict[site[0]]
                norr = group.orientations
                site[1] = rand.randrange(norr)
        self.set_latgas()

    def count(self, group_name, orientation):
        """

        Parameters
        ----------
        group_name: str
            The name of the group
        orientation:

        Returns
        -------

        """
        num_grp = []

        for defect_sublattice in self.defect_sublattices:
            num_grp.append(
                defect_sublattice.latgas_rep.count([group_name, orientation])
            )
        return num_grp

    def update_basestruct(self):
        basesites = self.matcher_base.get_mapping(self.structure, self.base_structure)
        idx = 0
        for i in basesites:
            self.base_structure[idx] = self.structure[i]
            idx += 1

    def defect_sublattice_structure(self, sublat_id):
        """

        Parameters
        ----------
        sublat_id: int
            index of sublattice

        Returns
        -------
        sublattice_structure: pymatgen.Structure
            sublattice structure object

        """
        assert sublat_id < self.n_sublat
        sublattice_structure = self.structure.copy()
        base_sites = self.matcher_base.get_mapping(self.structure, self.base_structure)
        sublattice_structure.remove_sites(base_sites)
        return sublattice_structure

    @property
    def vacancy_structure(self):
        filledsites = self.matcher_frame.get_mapping(
            self.perfect_structure, self.structure
        )
        vac_structure = self.perfect_structure.copy()
        vac_structure.remove_sites(filledsites)
        return vac_structure


class ObserverParams:
    def __init__(self):
        pass

    @classmethod
    def from_dict(cls, d):
        """

        Parameters
        ----------
        d: dict
            Dictionary

        Returns
        -------
        oparams: ObserverParams
            self
        """
        if "observer" in d:
            d = d["observer"]
        params = cls()
        # params.ignored_species = d["ignored_species"]
        return params

    @classmethod
    def from_toml(cls, f):
        """

        Parameters
        ----------
        f: str
            Name of input toml File

        Returns
        -------
        oparams : ObserverParams
            self
        """
        import toml

        return cls.from_dict(toml.load(f))
