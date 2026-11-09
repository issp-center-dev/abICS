# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2020- The University of Tokyo
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

"""``abics_al``: run the active-learning loop as one program.

Usage::

    abics_al [input.toml] [--max-iterations N]

Exit status 0: the target number of Monte Carlo runs exists; 2: stopped
because of the wall time (run again to continue); 1: a stage failed.
"""

import argparse
import datetime
import logging
import sys

import toml

from abics import __version__
import abics.loggers as loggers
from abics.applications.latgas_abinitio_interface.al_driver import DriverError, run_driver

logger = logging.getLogger("main")


def main():
    parser = argparse.ArgumentParser(
        description="Run the active-learning loop (abics_mlref -> abics_train -> abics_sampling -> ...) "
        "as a single program; see the [al] section of the input file."
    )
    parser.add_argument("tomlfile", nargs="?", default="input.toml", help="abICS input file (default: input.toml)")
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=None,
        help="number of Monte Carlo runs to complete (overrides [al] max_iterations)",
    )
    args = parser.parse_args()

    params_root = toml.load(args.tomlfile)
    log_params = dict(params_root.get("log", {}))
    # a single non-MPI process: log without MPI unless configured otherwise
    log_params.setdefault("console", "serial")
    log_params.setdefault("logfile_mode", "serial")
    loggers.set_log_handles(app_name="al_driver", level=logging.INFO, params=log_params)

    logger.info(f"Running abics_al (abICS v{__version__}) on {datetime.datetime.now()}")
    logger.info(f"-Reading input from: {args.tomlfile}")
    try:
        state = run_driver(params_root, args.tomlfile, args.max_iterations)
    except (DriverError, RuntimeError, ValueError) as e:
        logger.error(str(e))
        sys.exit(1)
    logger.info(f"Exiting on {datetime.datetime.now()} ({state}).")
    sys.exit(0 if state == "done" else 2)


if __name__ == "__main__":
    main()
