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
        self.base_input_dir = []
        self.solver = ""
        self.path = ""
        self.perturb = 0.0
        self.solver_run_scheme = ""
        self.properties = {}

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
        base_input_dir = d.get('base_input_dir', ['./baseinput'])
        if isinstance(base_input_dir, str):
            base_input_dir = [base_input_dir]
        params.base_input_dir = base_input_dir = list(map(lambda x: expand_path(x, os.getcwd()), base_input_dir))
        params.solver = d['type']
        params.path = expand_path(d['path'], os.getcwd())
        params.perturb = d.get('perturb', 0.1)
        params.solver_run_scheme = d.get('run_scheme',
                                         'mpi_spawn_ready')
        params.properties = d
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
