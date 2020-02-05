import numpy as np
import random as rand
import sys
import os
import copy

# from mpi4py import MPI

from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
import pymatgen.analysis.structure_analyzer as analy

from abics.mc import model
from abics.util import read_coords


def gauss(x, x0, sigma):
    """

    Parameters
    ----------
    x
    x0
    sigma

    Returns
    -------

    """
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
    mapping = []
    for i in range(len(lst)):
        if lst[i] != obj:
            mapping.append(i)
    return mapping


def match_latgas_group(latgas_rep, group):
    mapping = []
    for i in range(len(latgas_rep)):
        if latgas_rep[i][0] == group.name:
            mapping.append(i)
    return mapping


def g_r(structure, specie1, specie2, grid_1D):
    X = grid_1D.x
    dr = grid_1D.dx

    lattice = structure.lattice
    types_of_specie = [element.symbol for element in structure.types_of_specie]
    assert specie1 in types_of_specie
    assert specie2 in types_of_specie

    structure1 = structure.copy()
    structure2 = structure.copy()

    # print(specie1)
    not_specie1 = copy.copy(types_of_specie)
    not_specie1.remove(specie1)
    not_specie2 = types_of_specie
    not_specie2.remove(specie2)
    # print(not_specie1,not_specie2)

    structure1.remove_species(not_specie1)
    structure2.remove_species(not_specie2)

    num_specie1 = structure1.num_sites
    num_specie2 = structure2.num_sites

    dist = lattice.get_all_distances(structure1.frac_coords, structure2.frac_coords)
    dist_bin = np.around(dist / dr)
    # print(num_specie1,num_specie2)
    g = np.zeros(len(X))
    for i in range(len(X)):
        g[i] = np.count_nonzero(dist_bin == i + 1)

    if specie1 == specie2:
        pref = lattice.volume / (
            4.0 * np.pi * X * X * dr * num_specie1 * (num_specie1 - 1)
        )
    else:
        pref = lattice.volume / (4.0 * np.pi * X * X * dr * num_specie1 * num_specie2)

    g *= pref
    return g


class dft_latgas(model):
    """This class defines the DFT lattice gas mapping  model"""

    model_name = "dft_latgas"

    def __init__(
        self,
        abinitio_run,
        selective_dynamics=None,
        save_history=True,
        l_update_basestruct=False,
        check_ion_move=False,
        ion_move_tol=0.7,
    ):
        self.matcher = StructureMatcher(primitive_cell=False, allow_subset=False)
        self.abinitio_run = abinitio_run
        self.selective_dynamics = selective_dynamics
        self.save_history = save_history
        self.l_update_basestruct = l_update_basestruct
        self.check_ion_move = check_ion_move
        self.ion_move_tol = ion_move_tol

    def energy(self, config):
        """ Calculate total energy"""

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

        if self.selective_dynamics:
            seldyn_arr = [[True, True, True] for i in range(len(structure))]
            for specie in self.selective_dynamics:
                indices = structure.indices_from_symbol(specie)
                for i in indices:
                    seldyn_arr[i] = [False, False, False]
        else:
            seldyn_arr = None

        structure0 = structure
        energy, structure = self.abinitio_run.submit(
            structure, os.path.join(os.getcwd(), "output"), seldyn_arr
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
        """Construct the new configuration after the trial step is accepted"""
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

    def energy(self, config, save_history=False):
        rep_id = self.reps.index(tuple(config.latgas_rep))
        return np.float64(self.energy_list[rep_id])


class group(object):
    def __init__(self, name, species, coords=np.array([[[0.0, 0.0, 0.0]]])):
        self.name = name
        self.species = species
        self.coords = np.array(coords)
        self.orientations = len(coords)
        if self.orientations == 0:
            self.orientations = 1
        self.natoms = len(species)


class defect_sublattice(object):
    def __init__(self, site_centers, groups):
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
        site_centers = read_coords(d["coords"])
        groups = []
        for g in d["groups"]:
            name = g["name"]
            species = g.get("species", [name])
            coords = np.array(g.get("coords", [[[0.0, 0.0, 0.0]]]))
            groups.append(group(name, species, coords))
        return cls(site_centers, groups)


def base_structure(lat, dict_str):
    st = Structure(lat, [], [])
    if dict_str[0] == {}:
        return st
    for tc in dict_str:
        sp = tc["type"]
        coords = read_coords(tc["coords"])
        if len(coords.shape) == 1:
            coords = np.reshape(coords, (1, 3))
        n = coords.shape[0]
        for i in range(n):
            st.append(sp, coords[i, :])
    return st


class config:
    """This class defines the config with lattice gas mapping"""

    def __init__(
        self,
        base_structure,
        defect_sublattices,
        num_defects,
        cellsize=[1, 1, 1],
        perf_structure=None,
    ):
        """

        Parameters
        ----------
        base_structure
        defect_sublattices
        num_defects
        cellsize
        perf_structure
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
        if perf_structure:
            self.perf_structure = perf_structure
            self.perf_structure.make_supercell([cellsize[0], cellsize[1], cellsize[2]])
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
            for i in range(cellsize[0]):
                for j in range(cellsize[1]):
                    for k in range(cellsize[2]):
                        for l in range(site_centers.shape[0]):
                            defect_sublattice.site_centers_sc[idx] = site_centers[
                                l
                            ] + np.array([i, j, k])
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
        defect_sublattices

        Returns
        -------

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
                        group.coords[orr][j] + defect_sublattice.site_centers_sc[isite]
                        # properties={"velocities":[0,0,0]}
                    )

    def shuffle(self):
        for defect_sublattice in self.defect_sublattices:
            latgas_rep = defect_sublattice.latgas_rep
            rand.shuffle(latgas_rep)
            for site in latgas_rep:
                group = defect_sublattice.group_dict[site[0]]
                norr = group.orientations
                site[1] = rand.randrange(norr)
            # print(latgas_rep)
        self.set_latgas()

    def count(self, group_name, orientation):
        """

        Parameters
        ----------
        group_name: str

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
            self.perf_structure, self.structure
        )
        # print(self.perf_structure)
        # print(self.structure)
        # print(filledsites)
        vac_structure = self.perf_structure.copy()
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
