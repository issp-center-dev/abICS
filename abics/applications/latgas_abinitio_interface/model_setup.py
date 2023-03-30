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

from __future__ import annotations

import sys
import os
import copy
import itertools

import numpy as np
import numpy.random as rand
from scipy.special import gammaln

# from mpi4py import MPI

from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
import pymatgen.analysis.structure_analyzer as analy

from abics import __version__
from abics.exception import InputError
from abics.mc import Model
from abics.util import read_vector, read_matrix, read_tensor
from abics.applications.latgas_abinitio_interface import map2perflat

import logging
logger = logging.getLogger("main")

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
    return [ idx for idx,v in enumerate(lst) if v == obj ]

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
    return [ idx for idx,v in enumerate(lst) if v != obj ]

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
    return [ idx for idx,v in enumerate(latgas_rep) if v[0] == group.name ]

def perturb_structure(st: Structure, distance: float) -> None:
    """
    Perform random perturbation of the atomic coordinates.
    All atoms will be moved at the same distance.

    Which components will be perturbed is specified by a boolean array stored as st.site_properties["seldyn"].
    If not stored, all the components will be perturbed.
    Argument st will be mutated.

    Parameters
    ----------
    st: Structure
    distance: float
        strength of perturb
    """
    N = st.num_sites
    seldyn = np.array(st.site_properties.get("seldyn", np.ones((N, 3))), dtype=np.float64)
    assert seldyn.shape == (N, 3)
    seldyn *= rand.randn(N, 3)
    for i in range(N):
        r = seldyn[i, :]
        norm = np.linalg.norm(r)
        if norm != 0.0:
            st.sites[i].coords += seldyn[i, :] / norm * distance


class DFTLatticeGas(Model):
    """
    This class defines the DFT lattice gas mapping model
    """

    model_name = "DFTLatticeGas"

    def __init__(
        self,
        abinitio_run,
        save_history=True,
        l_update_basestruct=False,
        check_ion_move=False,
        ion_move_tol=0.7,
        enable_grandcanonical=False,
        gc_ratio=0.3,
        debug=False,
    ):
        """

        Parameters
        ----------
        abinitio_run:          runner object
            Runner (manager) of external solver program
        save_history:          boolean
        l_update_basestruct:   boolean
        check_ion_move:        boolean
        ion_move_tol:          float
        enable_grandcanonical: boolean
        """
        self.matcher = StructureMatcher(primitive_cell=False, allow_subset=False)
        self.abinitio_run = abinitio_run
        self.save_history = save_history
        self.l_update_basestruct = l_update_basestruct
        self.check_ion_move = check_ion_move
        self.ion_move_tol = ion_move_tol
        self.enable_grandcanonical = enable_grandcanonical
        self.gc_ratio = gc_ratio
        self.debug = debug

        if self.debug:
            logger.debug("debug mode enabled")
        logger.debug("enable_grandcanonical = {}".format(self.enable_grandcanonical))
        logger.debug("gc_ratio = {}".format(self.gc_ratio))

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
            # Try to avoid doing dft calculation for the same structure.
            # Assuming that calc_history is a list of ComputedStructureEntries
            for hist in config.calc_history:
                if self.matcher.fit(structure, hist[1]):
                    logger.info("match found in history")
                    config.structure = hist[2]
                    return hist[0]

        structure0 = structure
        energy, structure = self.abinitio_run.submit(
            structure, os.path.join(os.getcwd(), "output")
        )
        if self.check_ion_move:
            relax_analy = analy.RelaxationAnalyzer(structure0, structure)
            data = relax_analy.get_percentage_bond_dist_changes()
            for ion in self.check_ion_move:
                for index in structure0.indices_from_symbol(ion):
                    if any([v > self.ion_move_tol for v in data[index].values()]):
                        energy = float("inf")
                        logger.info("ion relaxed out of initial site")
                        break
                else:
                    continue
                break

        if self.save_history:
            config.calc_history.append((energy, structure0.copy(), structure.copy()))
            if len(config.calc_history) == 25:
                del config.calc_history[0:5]
        config.structure_norel = structure0
        config.structure = structure

        # grand potential
        if self.enable_grandcanonical:
            energy -= self._calc_muN_term(config)

        return np.float64(energy)

    def _calc_num_species(self, config, species):
        n = 0
        for sublat in config.defect_sublattices:
            n += len([ i for i,v in enumerate(sublat.latgas_rep) if v[0] == species ])
        return n

    def _calc_muN_term(self, config):
        muN = 0.0
        for chem in config.chemical_potential:
            sps = chem['species']
            mu = chem['mu']
            for sp in sps:
                muN += mu * self._calc_num_species(config, sp)
        return muN

    def _calc_weight(self, config):
        w = 0.0
        for sublat in config.defect_sublattices:
            nsite = len(sublat.latgas_rep)
            w += gammaln(nsite+1)
            for grp in sublat.groups:
                ngrp = len(match_latgas_group(sublat.latgas_rep, grp))
                w -= gammaln(ngrp+1)
        return w

    def _find_num_group_list(self, config):
        num = []
        for sublat in config.defect_sublattices:
            for grp in sublat.groups:
                ngrp = len(match_latgas_group(sublat.latgas_rep, grp))
                num.append(ngrp)
        return num

    def _calc_dweight(self, initial, final):
        w = 0.0
        for i,f in zip(initial, final):
            while i<f:
                w -= np.log(f)
                f -= 1
            while i>f:
                w += np.log(i)
                i -= 1
        return w

    def _try_rotate(self, defect_sublattice):
        latgas_rep = defect_sublattice.latgas_rep

        rot_ids = []
        for group in defect_sublattice.groups_orr:
            rot_ids += match_latgas_group(latgas_rep, group)

        if len(rot_ids) == 0:
            logger.debug("try_rotate: list empty. retry")
            if self.debug:
                return False, { 'mode': 'rotate', 'result': 'empty' }
            else:
                return False

        # Choose the site to rotate and rotate
        rot_id = rand.choice(rot_ids)
        rot_group = defect_sublattice.group_dict[latgas_rep[rot_id][0]]
        rot_group_orr = latgas_rep[rot_id][1]

        # Remove the current orientation from the list of possible orientations
        or_choices = [
            orx for orx in range(rot_group.orientations) if orx != rot_group_orr
        ]
        new_orr = rand.choice(or_choices)
        latgas_rep[rot_id] = [rot_group.name, new_orr]

        logger.debug("try_rotate: rotate id={}, group={}, from={}, to={}".format(
            rot_id, rot_group.name, rot_group_orr, new_orr))

        if self.debug:
            return True, { 'mode': 'rotate', 'result': 'ok', 'rot_id': rot_id, 'rot_group': rot_group.name, 'from': rot_group_orr, 'to': new_orr }
        else:
            return True

    def _try_exchange(self, defect_sublattice):
        latgas_rep = defect_sublattice.latgas_rep

        # Exchange different groups between sites
        ex1_group, ex2_group = rand.choice(
            defect_sublattice.groups, 2, replace=False
        )
        ex1_idlist = match_latgas_group(latgas_rep, ex1_group)
        ex2_idlist = match_latgas_group(latgas_rep, ex2_group)

        if len(ex1_idlist) == 0 or len(ex2_idlist) == 0:
            logger.debug("try_exchange: list empty {}({}), {}({}). retry".format(
                ex1_group.name, len(ex1_idlist), ex2_group.name, len(ex2_idlist)))
            if self.debug:
                return False, { 'mode': 'exchange', 'result': 'empty' }
            else:
                return False

        ex1_id = rand.choice(ex1_idlist)
        ex2_id = rand.choice(ex2_idlist)
        latgas_rep[ex1_id], latgas_rep[ex2_id] = (
            latgas_rep[ex2_id],
            latgas_rep[ex1_id],
        )

        logger.debug("try_exchange: exchange {}({}) and {}({})".format(
            ex1_id, ex1_group.name, ex2_id, ex2_group.name))

        if self.debug:
            return True, { 'mode': 'exchange', 'result': 'ok', 'ex1_id': ex1_id, 'ex1_group': ex1_group.name, 'ex2_id': ex2_id, 'ex2_group': ex2_group.name }
        else:
            return True

    def _try_add_remove(self, config, spec, do_add):
        if do_add:
            return self._try_add(config, spec)
        else:
            return self._try_remove(config, spec)

    def _try_remove(self, config, spec):
        # remove (a set of) particles
        logger.debug("src:")
        for sublat in config.defect_sublattices:
            logger.debug("state = {}".format([s for s,t in sublat.latgas_rep]))

        groups = spec[0]
        logger.debug("try_remove: groups={}".format(groups))

        # find sublattices which each species belongs to
        # sublats := [ (grp, [ sublat_id ]) ]
        sublats = [ (grp, [ sublat_id
                            for sublat_id, sublat in enumerate(config.defect_sublattices)
                            if grp in sublat.group_dict.keys() ])
                    for grp in groups ]

        # sublattice, and indexes of group and vacancy in latgas_rep
        # for each sublattice which each group belongs to
        # rep_table := [ (grp, [ (sublat_id, rep_index) ]) ]
        rep_table = []
        for grp, sublat_ids in sublats:
            tbl = []
            for sublat_id in sublat_ids:
                lat = config.defect_sublattices[sublat_id]
                idxs = match_latgas_group(lat.latgas_rep, lat.group_dict[grp])
                tbl += [ (sublat_id, idx) for idx in idxs ]
            rep_table += [(grp, tbl)]
        logger.debug("try_remove: rep_table={}".format(rep_table))

        # check if there are any choice
        is_unremovable = any([len(tbl) == 0 for grp, tbl in rep_table])
        if is_unremovable:
            logger.debug("try_remove: there is no choice")
            if self.debug:
                return False, { 'mode': 'remove', 'result': 'nochoice', 'spec': spec }
            else:
                return False

        # choose sublattice and latgas_rep index
        # for each species and vacancy associated with the species
        choice = []
        for grp, tbl in rep_table:
            sublat_id, idx = tbl[rand.choice(len(tbl))]
            choice += [(grp, sublat_id, idx)]
        logger.debug("try_remove: choice={}".format(choice))

        # update config
        for grp, sublat_id, idx in choice:
            lat = config.defect_sublattices[sublat_id]
            lat.latgas_rep[idx] = [lat.vac_group, 0]

        logger.debug("dst:")
        for sublat in config.defect_sublattices:
            logger.debug("state = {}".format([s for s,t in sublat.latgas_rep]))

        if self.debug:
            return True, { 'mode': 'remove', 'result': 'ok', 'spec': spec, 'choice': choice }
        else:
            return True

    def _try_add(self, config, spec):
        # add (a set of) particles
        logger.debug("src:")
        for sublat in config.defect_sublattices:
            logger.debug("  {}".format([s for s,t in sublat.latgas_rep]))

        groups = spec[1]
        logger.debug("try_add: groups={}".format(groups))

        # find sublattices which each species belongs to
        # sublats := [ (grp, [ sublat_id ]) ]
        sublats = [ (grp, [ sublat_id
                            for sublat_id, sublat in enumerate(config.defect_sublattices)
                            if grp in sublat.group_dict.keys() ])
                    for grp in groups ]

        # sublattice, and indexes of group and vacancy in latgas_rep
        # for each sublattice which each group belongs to
        # rep_table := [ (grp, [ (sublat_id, rep_index) ]) ]
        vac_table = []
        for grp, sublat_ids in sublats:
            tbl = []
            for sublat_id in sublat_ids:
                lat = config.defect_sublattices[sublat_id]
                idxs = match_latgas_group(lat.latgas_rep, lat.group_dict[lat.vac_group])
                tbl += [ (sublat_id, idx) for idx in idxs ]
            vac_table += [(grp, tbl)]
        logger.debug("try_add: vac_table={}".format(vac_table))

        # check if there are any choice
        is_unaddable = any([len(tbl) == 0 for grp, tbl in vac_table])
        if is_unaddable:
            logger.debug("try_add: there is no choice")
            if self.debug:
                return False, { 'mode': 'add', 'result': 'nochoice', 'spec': spec }
            else:
                return False

        # number of species
        nspecies = len(sublats)

        # number of sublattices to which any of species belong to
        nsublat = len(set(sum([ sublat_ids for grp,sublat_ids in sublats ], [])))

        # sublattices of each species are disjoint or not
        is_disjoint = all([ set(a[1]) & set(b[1]) == set() for a, b in itertools.combinations(sublats, 2) ])

        logger.debug("try_add: nspecies={}, nsublat={}, is_disjoint={}".format(nspecies, nsublat, is_disjoint))

        if nspecies == 1:
            # there is only one species: no conflict
            choice = []
            for grp, tbl in vac_table:
                sublat_id, idx = tbl[rand.choice(len(tbl))]
                choice += [(grp, sublat_id, idx)]
            logger.debug("try_add: nspecies=1 branch. choice={}".format(choice))

        elif nsublat == 1:
            # if there is one sublattice accommodating all species of concern,
            # set of vacancies for each species is identical.

            # assert identical
            assert(all([ set(a[1]) == set(b[1]) for a,b in itertools.combinations(vac_table, 2) ]))

            vac_list = vac_table[0][1]

            if len(vac_list) < nspecies:
                logger.debug("try_add: nsublat=1, there is no choice")
                if self.debug:
                    return False, { 'mode': 'add', 'result': 'nochoice', 'spec': spec }
                else:
                    return False

            ch = rand.choice(len(vac_list), nspecies, replace=False)
            logger.debug(">>> {}".format(ch))
            choice = [ (grp, *(vac_list[chidx])) for (grp,tbl),chidx in zip(vac_table, ch) ]
            logger.debug("try_add: nsublat=1 branch. choice={}".format(choice))

        elif is_disjoint:
            # if sublattices of each species are disjoint, may choose independently
            choice = []
            for grp, tbl in vac_table:
                sublat_id, idx = tbl[rand.choice(len(tbl))]
                choice += [(grp, sublat_id, idx)]
            logger.debug("try_add: is_disjoint branch. choice={}".format(choice))

        else:
            # general case. retry until consistent choice is made
            retry = 0
            max_retry = 5
            while True:
                # choose sublattice and latgas_rep index
                # for each species and vacancy associated with the species
                choice = []
                for grp, tbl in vac_table:
                    sublat_id, idx = tbl[rand.choice(len(tbl))]
                    choice += [(grp, sublat_id, idx)]

                # check consistency
                is_dup = False
                for a,b in itertools.combinations(choice, 2):
                    if a[1] == b[1] and a[2] == b[2]:
                        is_dup = True
                        break
                if is_dup == False:
                    break

                logger.debug("try_add: duplicated. retry. choice={}".format(choice))

                retry += 1
                if retry >= max_retry:
                    logger.debug("try_add: no choice. quit")
                    if self.debug:
                        return False, { 'mode': 'add', 'result': 'dup', 'spec': spec }
                    else:
                        return False
            logger.debug("try_add: general branch. choice={}".format(choice))

        # update config
        for grp, sublat_id, idx in choice:
            lat = config.defect_sublattices[sublat_id]
            lat.latgas_rep[idx] = [grp, rand.choice(lat.group_dict[grp].orientations)]

        logger.debug("dst:")
        for sublat in config.defect_sublattices:
            logger.debug("  {}".format([s for s,t in sublat.latgas_rep]))

        if self.debug:
            return True, { 'mode': 'add', 'result': 'ok', 'spec': spec, 'choice': choice }
        else:
            return True

    def _try_replace(self, config, spec):
        # replace set of particles to another set
        # single species version, i.e. A -> B
        logger.debug("src:")
        for sublat in config.defect_sublattices:
            logger.debug("  {}".format([s for s,t in sublat.latgas_rep]))

        sp_from = spec["from"]
        sp_to   = spec["to"]

        logger.debug("try_replace: groups from={}, to={}".format(sp_from, sp_to))

        if len(sp_from) != len(sp_to):
            logger.error("try_replace: supports only from single species to single species case")
            if self.debug:
                return False, { 'mode': 'replace', 'result': 'unsupported', 'spec': spec }
            else:
                return False

        # find sublattices which each species belongs to
        # sublats := [ (grp, [ sublat_id ]) ]
        sublats = [ (grp, [ sublat_id
                            for sublat_id, sublat in enumerate(config.defect_sublattices)
                            if grp in sublat.group_dict.keys() ])
                    for grp in sp_from ]

        # sublattice, and indexes of group and vacancy in latgas_rep
        # for each sublattice which each group belongs to
        # rep_table := [ (grp, [ (sublat_id, rep_index) ]) ]
        rep_table = []
        for grp, sublat_ids in sublats:
            tbl = []
            for sublat_id in sublat_ids:
                lat = config.defect_sublattices[sublat_id]
                idxs = match_latgas_group(lat.latgas_rep, lat.group_dict[grp])
                tbl += [ (sublat_id, idx) for idx in idxs ]
            rep_table += [(grp, tbl)]
        logger.debug("try_replace: rep_table={}".format(rep_table))

        # check if there are any choice
        is_unremovable = any([len(tbl) == 0 for grp, tbl in rep_table])
        if is_unremovable:
            logger.debug("try_replace: there is no choice")
            if self.debug:
                return False, { 'mode': 'replace', 'result': 'nochoice', 'spec': spec }
            else:
                return False

        # choose sublattice and latgas_rep index
        # for each species and vacancy associated with the species
        choice = []
        for grp, tbl in rep_table:
            sublat_id, idx = tbl[rand.choice(len(tbl))]
            choice += [(grp, sublat_id, idx)]
        logger.debug("try_replace: choice={}".format(choice))

        # assign new groups to vacant places
        group_assign = rand.permutation(sp_to)

        logger.debug("try_replace: assign={}".format(group_assign))

        # check if choice is valid
        noassign = False
        for k, newgrp in enumerate(group_assign):
            grp, sublat_id, idx = choice[k]

            if newgrp in config.defect_sublattices[sublat_id].group_dict.keys():
                pass
            else:
                noassign = True
                break
        if noassign:
            logger.debug("try_replace: assignment not satisfied")
            if self.debug:
                return False, { 'mode': 'replace', 'result': 'noassign', 'spec': spec }
            else:
                return False

        # update config
        for k, newgrp in enumerate(group_assign):
            grp, sublat_id, idx = choice[k]
            lat = config.defect_sublattices[sublat_id]
            lat.latgas_rep[idx] = [newgrp, rand.choice(lat.group_dict[newgrp].orientations)]

        logger.debug("dst:")
        for sublat in config.defect_sublattices:
            logger.debug("  {}".format([s for s,t in sublat.latgas_rep]))

        if self.debug:
            return True, { 'mode': 'replace', 'result': 'ok', 'spec': spec, 'choice': choice }
        else:
            return True

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

        logger.debug(">>> start trialstep")

        if self.enable_grandcanonical is True and config.chemical_potential is None:
            raise InputError("chemical potential missing")

        for subl in config.defect_sublattices:
            logger.debug("state = {}".format([ v[0] for v in subl.latgas_rep ]))
        logger.debug("energy = {:16.12e}".format(energy_now))

        e0 = energy_now

        if self.enable_grandcanonical:
            # w0 = self._calc_weight(config)
            numvec0 = self._find_num_group_list(config)

        # Back up structure and defect_sublattices
        structure0 = copy.deepcopy(config.structure)
        structure_norel0 = copy.deepcopy(config.structure_norel)
        defect_sublattices0 = copy.deepcopy(config.defect_sublattices)

        nretry = 0
        while True:
            # for defect_sublattice in [rand.choice(config.defect_sublattices)]: #config.defect_sublattices:
            logger.debug("=== trial {}".format(nretry))
            nretry += 1
            # print(latgas_rep)

            trial = True
            trial_info = None

            if not self.enable_grandcanonical or rand.rand() >= self.gc_ratio:
                # canonical move

                # If there is more than one group on the defect_sublattice,
                # we either change orientation of one group, or  exchange groups between sites

                defect_sublattice = rand.choice(config.defect_sublattices)
                #latgas_rep = defect_sublattice.latgas_rep

                if self.debug:
                    for i, x in enumerate(config.defect_sublattices):
                        if x == defect_sublattice:
                            sublat_id = i
                            break
                    else:
                        sublat_id = -1
                        logger.error("sublattice id not found")

                if len(defect_sublattice.groups) == 1:
                    logger.debug("trialstep: try rotate only")
                    if self.debug:
                        trial, trial_info = self._try_rotate(defect_sublattice)
                    else:
                        trial = self._try_rotate(defect_sublattice)
                elif not defect_sublattice.groups_orr:
                    logger.debug("trialstep: try exchange only")
                    if self.debug:
                        trial, trial_info = self._try_exchange(defect_sublattice)
                    else:
                        trial = self._try_exchange(defect_sublattice)
                else:
                    if rand.rand() < 0.5:
                        logger.debug("trialstep: try rotate")
                        if self.debug:
                            trial, trial_info = self._try_rotate(defect_sublattice)
                        else:
                            trial = self._try_rotate(defect_sublattice)
                    else:
                        logger.debug("trialstep: try exchange")
                        if self.debug:
                            trial, trial_info = self._try_exchange(defect_sublattice)
                        else:
                            trial = self._try_exchange(defect_sublattice)

                if self.debug:
                    trial_info['sublat_id'] = sublat_id

            else:
                # grandcanonical move

                # # choose a species set
                # #chem = rand.choice(config.chemical_potential)
                # chem_id = rand.randint(len(config.chemical_potential))
                # chem = config.chemical_potential[chem_id]
                # do_add = rand.rand() < 0.5

                # if do_add:
                #     logger.debug("trialstep: try grandcanonical move, type=add, chem={}".format(chem))
                #     if self.debug:
                #         trial, trial_info = self._try_add(config, chem)
                #     else:
                #         trial = self._try_add(config, chem)
                # else:
                #     logger.debug("trialstep: try grandcanonical move, type=remove, chem={}".format(chem))
                #     if self.debug:
                #         trial, trial_info = self._try_remove(config, chem)
                #     else:
                #         trial = self._try_remove(config, chem)


                # choose a move
                move_id = rand.randint(len(config.gc_move))
                move = config.gc_move[move_id]
                do_forward = (move['reverse'] == False) or (rand.rand() < 0.5)

                spec = (move['from'],move['to']) if do_forward else (move['to'],move['from'])

                if len(spec[0]) == 0:
                    # [] - >[A,..]
                    trial = self._try_add(config, spec)

                elif len(spec[1]) == 0:
                    # [A,..] -> []
                    trial = self._try_remove(config, spec)

                elif len(spec[0]) == len(spec[1]):
                    # [A] -> [B]
                    trial = self._try_replace(config, spec)

                else:
                    # most general case: to be implemented
                    raise RuntimeError("unsupported move")

                if self.debug:
                    trial, trial_info = trial

            if trial is False:
                #- retry
                # continue
                #- do not retry. let it be rejected
                break

            # print(latgas_rep)
            constraint_fulfilled = config.set_latgas()
            if constraint_fulfilled:
                break
            else:
                # restore from backup and retry
                config.structure = copy.deepcopy(structure0)
                config.defect_sublattices = copy.deepcopy(defect_sublattices0)
                continue

        if trial is True:
            # do vasp calculation on structure
            e1 = self.energy(config)

            if self.enable_grandcanonical:
                # w1 = self._calc_weight(config)
                # dW = w1 - w0
                numvec1 = self._find_num_group_list(config)
                dW = self._calc_dweight(numvec0, numvec1)
            else:
                dW = 0.0
        else:
            e1 = e0  # nominal value
            dW = float("nan")

        # return old structure
        structure = config.structure
        structure_norel = config.structure_norel
        config.structure = structure0
        config.structure_norel = structure_norel0
        defect_sublattices = config.defect_sublattices
        config.defect_sublattices = defect_sublattices0

        # Simply pass new structure and latgas_rep  in dconfig to be used by newconfig():
        dconfig = structure, structure_norel, defect_sublattices
        if e0 == float("inf"):
            dE = e1
        else:
            dE = e1 - e0

        for subl in defect_sublattices:
            logger.debug("state = {}".format([ v[0] for v in subl.latgas_rep ]))
        logger.debug("energy = {:16.12f}".format(e1))

        logger.debug("trialstep: dE={:16.12e}, dW={:.5e}".format(dE, dW))
        logger.debug("<<< end trialstep")

        if self.debug:
            return dconfig, dE, dW, trial_info
        else:
            return dconfig, dE, dW


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
        config.structure, config.structure_norel, config.defect_sublattices = dconfig
        if self.l_update_basestruct:
            self.update_basestruct(config)
        return config


class EnergyList(DFTLatticeGas):
    def __init__(
        self,
        calcode,
        vasp_run,
        base_vaspinput,
        matcher_base,  # matcher, matcher_site,
        queen,
        reps,
        energy_lst,
        matcher=None,
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


class Group:
    def __init__(
        self, name, species, *, coords=None, relaxations=None, magnetizations=None
    ):
        """

        Parameters
        ----------
        name: str
            The name of atomic group
        species: list of str
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
        self.relaxations = (
            np.array(relaxations)
            if relaxations is not None
            else np.ones((1, 3), dtype=bool)
        )
        self.magnetizations = (
            np.array(magnetizations) if magnetizations is not None else np.zeros(1)
        )
        self.orientations = len(self.coords)
        if self.orientations == 0:
            self.orientations = 1
        self.natoms = len(species)


class DefectSublattice:
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
        self.groups_orr = [ group for group in groups if group.orientations > 1 ]
        self.group_dict = { group.name: group for group in groups }

    @classmethod
    def from_dict(cls, d) -> "DefectSublattice":
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

            coords = read_tensor(g.get("coords", [[[0, 0, 0]]]), rank=3)
            m = coords.shape[1]
            if n != 0 and m != n:
                raise InputError(
                    'number of atoms mismatch in group [{}]: "species"={}, "coords"={}'.format(
                        name, n, m
                    )
                )

            relaxation = read_matrix(g.get("relaxation", np.ones((n, 3))), dtype=bool)
            m = relaxation.shape[0]
            if n != 0 and m != n:
                raise InputError(
                    'number of atoms mismatch in group [{}]: "species"={}, "relaxation"={}'.format(
                        name, n, m
                    )
                )

            mag = read_vector(g.get("magnetization", np.zeros(n)))
            m = len(mag)
            if n != 0 and m != n:
                raise InputError(
                    'number of atoms mismatch in group [{}]: "species"={}, "magnetization"={}'.format(
                        name, n, m
                    )
                )
            groups.append(
                Group(
                    name,
                    species,
                    coords=coords,
                    relaxations=relaxation,
                    magnetizations=mag,
                )
            )
        return cls(site_centers, groups)


def base_structure(lat, dict_str) -> Structure:
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
            site_properties={
                "seldyn": np.zeros((0, 3), dtype=bool),
                "magnetization": np.zeros(0),
            },
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


class Config:
    """This class defines the config with lattice gas mapping"""

    def __init__(
        self,
        base_structure: Structure,
        defect_sublattices,
        num_defects,
        cellsize=[1, 1, 1],
        constraint_func=bool,
        constraint_energy=None,
        perfect_structure=None,
        chemical_potential=None,
        gc_move=None,
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
        constraint_func : function, optional
            function that takes pymatgen.Structure as argument and returns True
            when constraints are satisfied, by default bool
        perfect_structure : pymatgen.Structure, optional
            Strucure of all sites (union of base and defect), by default None
        chemical_potential : dict
            { tuple of species names: chemical potential }
        gc_move : list
            grandcanonical moves: { species from, species to, reversible }
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

        self.calc_history: list = []
        self.cellsize = cellsize
        self.base_structure = base_structure
        self.constraint_func = constraint_func
        self.constraint_energy = constraint_energy
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
            for (idx, (i, j, k, l)) in enumerate(
                itertools.product(
                    range(cellsize[0]),
                    range(cellsize[1]),
                    range(cellsize[2]),
                    range(site_centers.shape[0]),
                )
            ):
                defect_sublattice.site_centers_sc[idx] = site_centers[l] + np.array(
                    [i, j, k]
                )
                idx += 1
            defect_sublattice.site_centers_sc /= np.array(cellsize)
            num_sites += len(defect_sublattice.site_centers_sc)

            # add vacancy
            num_sites_sublat = len(defect_sublattice.site_centers_sc)
            num_defects_sublat = 0
            for group, num in num_defects[sublat_id].items():
                num_defects_sublat += num
            # check if vacancy group has already been defined
            for grp in defect_sublattice.groups:
                if len(grp.species) == 0:
                    vac_name = grp.name
                    break
            else:
                # otherwise, create a new group for vacancy
                vac_name = '@'
                vac_group = Group(name=vac_name, species=[])
                defect_sublattice.groups.append(vac_group)
                defect_sublattice.group_dict[vac_name] = vac_group
            # append vacancy even if sites are full
            nvac = num_sites_sublat - num_defects_sublat
            num_defects[sublat_id][vac_name] = num_defects[sublat_id].get(vac_name, 0) + nvac
            defect_sublattice.vac_group = vac_name

            # change cartesian coordinates to fractional
            for group in defect_sublattice.groups:
                for i in range(group.orientations):
                    for j in range(group.natoms):
                        group.coords[i][j] = np.dot(group.coords[i][j], invSuper)

            # fill the lattice gas representation list
            latgas_rep = []
            for group, num in num_defects[sublat_id].items():
                latgas_rep.extend(
                    [[group, 0] for i in range(num)]
                )
                ntot_defects += num
            defect_sublattice.latgas_rep = latgas_rep

            sublat_id += 1

        self.defect_sublattices = defect_sublattices

        # # table group.name to sublattice which the group belongs to
        # self.sublattice_dict = {}
        # for defect_sublattice in self.defect_sublattices:
        #     for grp in defect_sublattice.groups:
        #         self.sublattice_dict[grp.name] = self.sublattice_dict.get(grp.name, []) + [defect_sublattice]

        assert num_sites == ntot_defects
        constraint_fullfilled = self.set_latgas()
        if not constraint_fullfilled:
            self.shuffle()

        # store chemical potential parameters and grandcanonical move
        self.chemical_potential = chemical_potential
        self.gc_move = gc_move

    def set_latgas(self, defect_sublattices=False):
        """

        Parameters
        ----------
        defect_sublattices: pymatgen.Structure

        Returns
        -------
        True if resulting structure satisifies constraints, False if not
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
                        properties={
                            "seldyn": group.relaxations[j, :],
                            "magnetization": group.magnetizations[j],
                        },
                    )
        return self.constraint_func(self.structure)

    def reset_from_structure(self, st_in):
        st = self.dummy_structure()
        st = map2perflat(st, st_in)
        st.remove_species(["X"])

        for defect_sublattice in self.defect_sublattices:    
            # Extract group information for this sublattice
            d_sp2grp = {}
            sp = set()
            for group in defect_sublattice.groups:
                if len(group.species) > 1:
                    raise InputError(
                        'Cannot set initial structure when multi-atom groups are used'
                        )
                if len(group.species) == 0:
                    d_sp2grp["X0+"] = group.name
                    sp.add("X")
                else:
                    d_sp2grp[group.species[0]] = group.name
                    sp.add(group.species[0])

            # Map initial structure to sublattice
            sublattice_dummy_st = self.dummy_structure_from_sublattice(defect_sublattice)
            d_matrix = st.lattice.get_all_distances(st.frac_coords, sublattice_dummy_st.frac_coords)
            mapping = np.where(d_matrix < 1e-4)
            num_match = len(mapping[0])
            for i in range(num_match): 
                sublattice_dummy_st.replace(mapping[1][i], st[mapping[0][i]].species_string)
            if not sp.issuperset(set(sublattice_dummy_st.symbol_set)):
                raise InputError(
                    "Initial structure contains species {} not specified in config".format(
                        sp ^ set(sublattice_dummy_st.symbol_set)
                        ) + \
                    "\n(if species 'X' is listed here, it means the number of atoms don't match)"
                )
                
            # Set lattice gas representation list
            latgas_rep = []
            for isite, site in enumerate(sublattice_dummy_st):
                latgas_rep.append([d_sp2grp[site.species_string], 0])
            defect_sublattice.latgas_rep = latgas_rep

        constraint_fullfilled = self.set_latgas()
        if not constraint_fullfilled:
            raise InputError(
                "initial structure violates constraints specified by constraint_func"
                )

    def dummy_structure(self):
        """

        Parameters
        ----------

        Returns
        -------
        Structure where all atoms and vacancies are replaced by dummy atoms
        """
        dummy_structure = copy.deepcopy(self.base_structure)
        for defect_sublattice in self.defect_sublattices:
            latgas_rep = defect_sublattice.latgas_rep
            assert len(latgas_rep) == len(defect_sublattice.site_centers_sc)
            for isite in range(len(latgas_rep)):
                grp_name = latgas_rep[isite][0]
                orr = latgas_rep[isite][1]
                group = defect_sublattice.group_dict[grp_name]
                if group.natoms == 0:
                    dummy_structure.append(
                        "X",
                        defect_sublattice.site_centers_sc[isite],
                        properties={
                            "seldyn": (True, True, True),
                            "magnetization": (0, 0, 0),
                        },
                    )
                for j in range(group.natoms):
                    dummy_structure.append(
                        "X",
                        group.coords[orr][j] + defect_sublattice.site_centers_sc[isite],
                        properties={
                            "seldyn": group.relaxations[j, :],
                            "magnetization": group.magnetizations[j],
                        },
                    )
        return dummy_structure

    def dummy_structure_sp(self, species_in):
        """

        Parameters
        ----------
        species (str): name of species for constructing dummy lattice

        Returns
        -------
        Structure where all atoms and vacancies are replaced by dummy atoms
        """
        dummy_structure = copy.deepcopy(self.base_structure)
        sp_remove = filter(lambda sp: sp != species_in, dummy_structure.symbol_set)
        dummy_structure.remove_species(sp_remove)
        for defect_sublattice in self.defect_sublattices:
            latgas_rep = defect_sublattice.latgas_rep
            assert len(latgas_rep) == len(defect_sublattice.site_centers_sc)
            # find all species that can reside on this lattice
            species_sublattice = set()
            groups = defect_sublattice.group_dict.keys()
            for grp_name in groups:
                group = defect_sublattice.group_dict[grp_name]
                if group.natoms > 1:
                    logger.info("dummy_structure_sp does not support multi-atom groups")
                if group.natoms == 1:
                    species_sublattice.add(group.species[0])
            if species_in not in species_sublattice:
                continue
            for isite in range(len(latgas_rep)):
                dummy_structure.append(
                    "X",
                    defect_sublattice.site_centers_sc[isite],
                    properties={
                        "seldyn": (True, True, True),
                        "magnetization": (0, 0, 0),
                    },
                )

        return dummy_structure

    def dummy_structure_from_sublattice(self, defect_sublattice):
        """

        Parameters
        ----------
        defect_sublattice (str): defect_sublattice for constructing dummy lattice

        Returns
        -------
        Structure where all atoms and vacancies are replaced by dummy atoms
        """
        dummy_structure = Structure(
            lattice=self.base_structure.lattice,
            species=[],
            coords=[],
            )
        

        latgas_rep = defect_sublattice.latgas_rep
        assert len(latgas_rep) == len(defect_sublattice.site_centers_sc)
        # find all species that can reside on this lattice
        species_sublattice = set()
        groups = defect_sublattice.group_dict.keys()
        for grp_name in groups:
            group = defect_sublattice.group_dict[grp_name]
            if group.natoms > 1:
                print("dummy_structure_sp does not support multi-atom groups")
                sys.exit(1)
            if group.natoms == 1:
                species_sublattice.add(group.species[0])
        for isite in range(len(latgas_rep)):
            dummy_structure.append(
                "X",
                defect_sublattice.site_centers_sc[isite],
                properties={
                    "seldyn": (True, True, True),
                    "magnetization": (0, 0, 0),
                },
            )
                
        return dummy_structure

    def shuffle(self):
        max_trial = 1000
        num_trial = 0
        while num_trial < max_trial:
            for defect_sublattice in self.defect_sublattices:
                latgas_rep = defect_sublattice.latgas_rep
                rand.shuffle(latgas_rep)
                for site in latgas_rep:
                    group = defect_sublattice.group_dict[site[0]]
                    norr = group.orientations
                    site[1] = rand.randint(norr)
            if self.set_latgas(): 
                return 0, 'Configuration initialized randomly'
            else:
                if self.constraint_energy: 
                    e = self.constraint_energy(self.structure)
                    num_trial = 0
                    while num_trial < max_trial: 
                        num_trial += 1
                        defect_sublattice = rand.choice(self.defect_sublattices)
                        latgas_rep = defect_sublattice.latgas_rep
                        i0, i1 = rand.choice(len(latgas_rep), 2)
                        latgas_rep[i0], latgas_rep[i1] = latgas_rep[i1], latgas_rep[i0]
                        if self.set_latgas():
                            msg = "---constraint_module.constraint_energy was used" + \
                                " to find configuration that follows constraints with {} trials".format(num_trial)
                            return 0, msg
                        e1 = self.constraint_energy(self.structure)
                        if e1 == 0:
                            msg = "There's something wrong: constraint_func and constraint_energy are inconsistent"
                            return 1, msg
                        if e1 > e:
                            latgas_rep[i0], latgas_rep[i1] = latgas_rep[i1], latgas_rep[i0]
                        else:
                            e = e1
                    msg = "Failed to find configuration with constraint_module.constraint_energy after {} trials".format(num_trial)
                    return 1, msg
                    
                else:
                    num_trial += 1
        msg = "Failed to find configuration that follows constraints by random shuffling. Try setting constraint_energy."       
        return 1, msg



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

        d = toml.load(f)
        return cls.from_dict(d["observer"])


# For backward compatibility
if __version__ < "3":
    dft_latgas = DFTLatticeGas
    energy_lst = EnergyList
    group = Group
    defect_sublattice = DefectSublattice
    config = Config
