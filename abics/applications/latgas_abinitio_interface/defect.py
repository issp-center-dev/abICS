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

from typing import Any, MutableMapping
import os, sys
import numpy as np

from pymatgen.core import Lattice, Structure

from .model_setup import Config, DefectSublattice, base_structure

from abics.exception import InputError
from abics.util import read_matrix, expand_path

import logging
logger = logging.getLogger("main")


class DFTConfigParams:

    lat: Lattice
    supercell: list[int]
    base_structure: Structure
    defect_sublattices: list[DefectSublattice]
    num_defects: list[dict[str, int]]

    def __init__(self, dconfig: MutableMapping[str, Any]) -> None:
        """
        Get information from dictionary

        Parameters
        ----------
        dconfig: dict
            Dictionary
        """

        if "unitcell" not in dconfig:
            raise InputError('"config.unitcell" is not found in the "config" section.')
        self.lat = Lattice(read_matrix(dconfig["unitcell"]))

        self.supercell = dconfig.get("supercell", [1, 1, 1])
        if not isinstance(self.supercell, list):
            raise InputError('"config.supercell" should be a list of integers')

        constraint_module = dconfig.get("constraint_module", False)
        if not isinstance(constraint_module, bool):
            raise InputError('"config.constraint_module" should be true or false')
        self.constraint_energy = None
        if constraint_module:
            sys.path.append(os.getcwd())
            import constraint_module
            self.constraint_func = constraint_module.constraint_func
            if 'constraint_energy' in dir(constraint_module):
                self.constraint_energy = constraint_module.constraint_energy
        else:
            self.constraint_func = bool

        # alternative: provide module and function name as a parameter value
        constraint_package = dconfig.get("constraint", None)
        if constraint_package:
            import importlib
            _cpkg = constraint_package.split('.')
            _cpkg_mod_name = '.'.join(_cpkg[0:-1])
            _cpkg_func_name = _cpkg[-1]
            try:
                _cpkg_mod = importlib.import_module(_cpkg_mod_name)
                _cpkg_func = getattr(_cpkg_mod, _cpkg_func_name)
            except BaseException:
                raise InputError("constraint module could not be loaded.")
            self.constraint_func = _cpkg_func #overwrite

        init_structure_path = dconfig.get("init_structure", None)
        if init_structure_path is not None:
            init_structure_path = expand_path(init_structure_path, os.getcwd())
            if not os.path.isfile(init_structure_path):
                raise InputError(
                    'cannot find {}; "config.init_structure" should be path to a structure file'.format(init_structure_path)
                    )
            try:
                self.init_structure = Structure.from_file(init_structure_path)
            except:
                raise InputError('failed reading {}; check that the file is pymatgen-compatible.'.format(init_structure_path))
        else:
            self.init_structure = None

        if "base_structure" not in dconfig:
            raise InputError('"base_structure" is not found in the "config" section.')
        self.base_structure = base_structure(self.lat, dconfig["base_structure"])

        if "defect_structure" not in dconfig:
            raise InputError('"defect_structure" is not found in the "config" section.')
        self.defect_sublattices = [
            DefectSublattice.from_dict(ds) for ds in dconfig["defect_structure"]
        ]
        self.num_defects = [
            {g["name"]: g["num"] for g in ds["groups"]}
            for ds in dconfig["defect_structure"]
        ]

        if "chemical_potential" in dconfig:
            self.chemical_potential = []
            for entry in dconfig["chemical_potential"]:
                species = entry.get("species", None)
                mu = entry.get("mu", None)
                if type is None or mu is None:
                    raise InputError("invalid input for chemical potentials")
                species = [species] if type(species) is str else species
                self.chemical_potential.append({'species': species, 'mu': mu})
        else:
            self.chemical_potential = None
        logger.debug("chemical_potential = {}".format(self.chemical_potential))

        def _uniq_count(xs):
            return [(i,j) for i,j in zip(*np.unique(xs,return_counts=True))]

        if "grandcanonical_move" in dconfig:
            self.gc_move = []
            for entry in dconfig["grandcanonical_move"]:
                if "species" in entry.keys():
                    species = entry["species"]
                    species = [ species ] if type(species) is not list else species
                    self.gc_move.append({ 'from': _uniq_count(species), 'to': [], 'reverse': True })
                else:
                    sp_from = entry.get("from", [])
                    sp_to   = entry.get("to", [])
                    #sp_rev  = entry.get("reverse", True)
                    sp_rev  = True  # consider only reversible case
                    if sp_from is [] and sp_to is []:
                        raise("invalid input for grandcanonical_move")
                    sp_from = [ sp_from ] if type(sp_from) is not list else sp_from
                    sp_to   = [ sp_to ]   if type(sp_to)   is not list else sp_to
                    self.gc_move.append({ 'from': _uniq_count(sp_from), 'to': _uniq_count(sp_to), 'reverse': sp_rev })
        else:
            if self.chemical_potential is None:
                self.gc_move = None
            else:
                # taken from chemical_potential table: assume add/remove of species
                self.gc_move = []
                for entry in self.chemical_potential:
                    species = entry["species"]
                    self.gc_move.append({ 'from': _uniq_count(species), 'to': [], 'reverse': True })
        logger.debug("grandcanonical_move = {}".format(self.gc_move))

    @classmethod
    def from_dict(cls, d: MutableMapping) -> "DFTConfigParams":
        return cls(d)

    @classmethod
    def from_toml(cls, f: str) -> "DFTConfigParams":
        """

        Get information from toml file.

        Parameters
        ----------
        f: str
            Name of input toml File

        Returns
        -------
        cls.from_dict: DFTConfigParams object
            Parameters for DFT solver
        """
        import toml

        d = toml.load(f)
        return cls(d["config"])


def defect_config(cparam: DFTConfigParams) -> Config:
    """
    Get configuration information

    Parameters
    ----------
    cparam: DFTConfigParams object
        Parameters for DFT solver

    Returns
    -------
    spinel_config: config object
        spinel configure object
    """
    spinel_config = Config(
        base_structure = cparam.base_structure,
        defect_sublattices = cparam.defect_sublattices,
        num_defects = cparam.num_defects,
        cellsize = cparam.supercell,
        constraint_func = cparam.constraint_func,
        constraint_energy = cparam.constraint_energy,
        chemical_potential = cparam.chemical_potential,
        gc_move = cparam.gc_move,
    )
    spinel_config.shuffle()
    return spinel_config
