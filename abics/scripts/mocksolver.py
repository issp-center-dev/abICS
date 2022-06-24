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

"""
Mock for energy calculator
"""

from abics.applications.latgas_abinitio_interface.mocksolver import MockSolver

def main_impl(output_dir):
    solver = MockSolver()
    solver.calc_energy("", output_dir)

def main():
    import sys
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    main_impl(output_dir)

if __name__ == "__main__":
    main()
