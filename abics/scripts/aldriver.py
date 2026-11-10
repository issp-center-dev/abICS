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

    abics_al [input.toml] [--max-iterations N] [--evaluate-final-model] [--forget-durations]

Exit status 0: the target number of Monte Carlo runs exists; 2: stopped
because of the wall time (run again to continue); 1: a stage failed or
the input is wrong.
"""

import argparse
import datetime
import logging
import os
import sys

import toml

from abics import __version__
import abics.loggers as loggers
from abics.applications.latgas_abinitio_interface.al_driver import DRIVER_DIR, DriverError, run_driver

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
        help="number of Monte Carlo runs to complete, >= 1 (overrides [al] max_iterations)",
    )
    parser.add_argument(
        "--evaluate-final-model",
        action="store_true",
        help="after the last Monte Carlo run, run one more reference-data step to evaluate the final model "
        "(same as [al] evaluate_final_model = true)",
    )
    parser.add_argument(
        "--forget-durations",
        action="store_true",
        help="forget the recorded stage durations used for the wall-time decisions",
    )
    args = parser.parse_args()
    if args.max_iterations is not None and args.max_iterations < 1:
        parser.error("--max-iterations must be >= 1")

    try:
        params_root = toml.load(args.tomlfile)
    except (OSError, toml.TomlDecodeError) as e:
        print(f"abics_al: cannot read {args.tomlfile}: {e}", file=sys.stderr)
        sys.exit(1)

    # The driver logs into its own file: the stage programs started from the same
    # input file would truncate a shared [log] logfile_path.
    os.makedirs(DRIVER_DIR, exist_ok=True)
    log_params = dict(params_root.get("log", {}))
    shared_logfile = log_params.get("logfile_path")
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    log_params["logfile_path"] = os.path.join(DRIVER_DIR, f"abics_al_{stamp}.log")
    log_params.setdefault("console", "serial")
    log_params["logfile_mode"] = "serial"
    loggers.set_log_handles(app_name="al_driver", level=logging.INFO, params=log_params)

    logger.info(f"Running abics_al (abICS v{__version__}) on {datetime.datetime.now()}")
    logger.info(f"-Reading input from: {args.tomlfile}")
    logger.info(f"-Log file: {log_params['logfile_path']}")
    if shared_logfile:
        logger.warning(f"-[log] logfile_path = '{shared_logfile}' is used by abics_train and abics_sampling, "
                       "each of which overwrites it; their full output is kept in al_driver/*.log")
    try:
        state = run_driver(
            params_root,
            args.tomlfile,
            args.max_iterations,
            args.forget_durations,
            True if args.evaluate_final_model else None,
        )
    except DriverError as e:
        logger.error(str(e))
        sys.exit(1)
    except Exception as e:  # noqa: BLE001 - message first, traceback kept in the log
        logger.error(f"abics_al stopped: {type(e).__name__}: {e}")
        logger.error("traceback:", exc_info=True)
        sys.exit(1)
    logger.info(f"Exiting on {datetime.datetime.now()} ({state}).")
    sys.exit(0 if state == "done" else 2)


if __name__ == "__main__":
    main()
