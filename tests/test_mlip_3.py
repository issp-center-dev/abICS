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
import tempfile
import unittest

from abics.applications.latgas_abinitio_interface.base_solver import create_solver
from abics.applications.latgas_abinitio_interface.mlip_3 import MLIP3Solver
from abics.applications.latgas_abinitio_interface.params import DFTParams


class TestMLIP3(unittest.TestCase):
    def setUp(self):
        params = DFTParams.from_dict({"type": "mlip_3", "path": ".", "run_scheme": "subprocess"})
        self.solver = create_solver(params.solver, params)

    def test_create_solver(self):
        self.assertIsInstance(self.solver, MLIP3Solver)

    def test_rerun_is_not_implemented(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(NotImplementedError) as cm:
                self.solver.input.update_info_from_files(tmpdir, 1)
            self.assertIn(os.path.basename(tmpdir), str(cm.exception))


if __name__ == "__main__":
    unittest.main()
