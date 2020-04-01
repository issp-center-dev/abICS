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

import os

from ...util import expand_path

class DFTParams:
    def __init__(self):
        pass

    @classmethod
    def from_dict(cls, d):
        """
        Read information from dictionary

        Parameters
        ----------
        d: dict
            Dictionary

        Returns
        -------
        params: DFTParams object
            self
        """
        if 'solver' in d:
            d = d['solver']
        params = cls()
        params.base_input_dir = expand_path(d.get('base_input_dir', './baseinput'), os.getcwd())
        params.solver = d['type']
        params.path = expand_path(d['path'], os.getcwd())
        params.perturb = d.get('perturb', 0.1)
        params.solver_run_scheme = d.get('run_scheme',
                                         'mpi_spawn_ready')
        return params

    @classmethod
    def from_toml(cls, f):
        """
        Read information from toml file

        Parameters
        ----------
        f: str
            Name of input toml File

        Returns
        -------
        oDFTParams: DFTParams object
            self
        """
        import toml
        return cls.from_dict(toml.load(f))
