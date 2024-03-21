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

from typing import MutableMapping

import sys
import datetime

import toml

from abics import __version__

import logging
import abics.loggers as loggers
# from pathlib import Path

logger = logging.getLogger("main")

def main_impl(params: MutableMapping):
    solver_type = params["sampling"]["solver"]["type"]
    if solver_type == "potts":
        raise NotImplementedError("Potts solver is not implemented yet.")
    else:
        import abics.scripts.postproc_dft_latgas
        abics.scripts.postproc_dft_latgas.postproc_dft_latgas(params)


def main():
    now = datetime.datetime.now()

    tomlfile = sys.argv[1] if len(sys.argv) > 1 else "input.toml"
    params = toml.load(tomlfile)

    loggers.set_log_handles(
        app_name = "postproc",
        level = logging.INFO,
        logfile_path = None,
        logfile_mode = "master",
        params=params.get("log", {}))

    logger.info(f"Running abics_postproc (abICS v{__version__}) on {now}")
    logger.info("-Reading input from: {}".format(tomlfile))

    main_impl(params)


if __name__ == "__main__":
    main()
