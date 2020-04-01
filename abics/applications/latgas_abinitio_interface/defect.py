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

from pymatgen import Lattice

from .model_setup import config, defect_sublattice, base_structure

from abics.exception import InputError
from abics.util import read_matrix


class DFTConfigParams:
    def __init__(self, dconfig):
        """
        Get information from dictionary

        Parameters
        ----------
        dconfig: dict
            Dictionary
        """

        if "config" in dconfig:
            dconfig = dconfig["config"]

        if "unitcell" not in dconfig:
            raise InputError('"unitcell" is not found in the "config" section.')
        self.lat = Lattice(read_matrix(dconfig["unitcell"]))

        self.supercell = dconfig.get("supercell", [1, 1, 1])
        if "base_structure" not in dconfig:
            raise InputError('"base_structure" is not found in the "config" section.')
        self.base_structure = base_structure(self.lat, dconfig["base_structure"])

        if "defect_structure" not in dconfig:
            raise InputError('"defect_structure" is not found in the "config" section.')
        self.defect_sublattices = [
            defect_sublattice.from_dict(ds) for ds in dconfig["defect_structure"]
        ]
        self.num_defects = [
            {g["name"]: g["num"] for g in ds["groups"]} for ds in dconfig["defect_structure"]
        ]

    @classmethod
    def from_dict(cls, d):
        return cls(d)

    @classmethod
    def from_toml(cls, f):
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

        return cls(toml.load(f))


def defect_config(cparam: DFTConfigParams):
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
    spinel_config = config(
        cparam.base_structure,
        cparam.defect_sublattices,
        cparam.num_defects,
        cparam.supercell,
    )
    spinel_config.shuffle()
    return spinel_config
