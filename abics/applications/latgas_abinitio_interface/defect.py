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

from pymatgen.core import Lattice, Structure

from .model_setup import Config, DefectSublattice, base_structure

from abics.exception import InputError
from abics.util import read_matrix


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
        if constraint_module:
            import os, sys

            sys.path.append(os.getcwd())
            from constraint_module import constraint_func

            self.constraint_func = constraint_func
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
        cparam.base_structure,
        cparam.defect_sublattices,
        cparam.num_defects,
        cparam.supercell,
        cparam.constraint_func,
    )
    spinel_config.shuffle()
    return spinel_config
