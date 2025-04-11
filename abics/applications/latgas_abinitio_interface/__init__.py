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

from .map2perflat import *

from .base_solver import register_solver
from .base_trainer import register_trainer

register_solver("vasp", "VASPSolver", "abics.applications.latgas_abinitio_interface.vasp")
register_solver("qe", "QESolver", "abics.applications.latgas_abinitio_interface.qe")
register_solver("openmx", "OpenMXSolver", "abics.applications.latgas_abinitio_interface.openmx")
register_solver("aenet", "AenetSolver", "abics.applications.latgas_abinitio_interface.aenet")
register_solver("aenetPyLammps", "AenetPyLammpsSolver", "abics.applications.latgas_abinitio_interface.aenet_pylammps")
register_solver("nequip", "NequipSolver", "abics.applications.latgas_abinitio_interface.nequip")
register_solver("mlip_3", "MLIP3Solver", "abics.applications.latgas_abinitio_interface.mlip_3")
register_solver("User", "UserFunctionSolver", "abics.applications.latgas_abinitio_interface.user_function_solver")

register_trainer("aenet", "AenetTrainer", "abics.applications.latgas_abinitio_interface.aenet_trainer")
register_trainer("nequip", "NequipTrainer", "abics.applications.latgas_abinitio_interface.nequip_trainer")
register_trainer("mlip_3", "MLIP3Trainer", "abics.applications.latgas_abinitio_interface.mlip_3_trainer")