"""
Fault-injecting variant of ``mock.py`` for the fault-injection test
(``run_faults.py``).

It computes the same energy as ``mock.py`` but first reads ``faults.json``
next to this file and misbehaves as told there.  Keys (sample indices are
the numbers of the ``input<i>`` directories; ``"all"`` means every sample)::

    {
      "fail_first":   [0, 3],   # exit with status 1 on the first attempt
      "always_fail":  [5],      # exit with status 1 on every attempt
      "hang_first":   [2],      # sleep hang_seconds on the first attempt (to be killed by the timeout)
      "hang_seconds": 60,
      "sleep":        3,        # seconds every run takes (to make a stage long enough to be interrupted)
      "fail_delay":   3         # seconds before an injected failure (default 3, see below)
    }

The attempt number is kept in ``attempt.count`` in the working directory,
which abics_mlref runs in place (no scratch_dir), so it survives a retry.

An injected failure happens ``fail_delay`` seconds after the start: a task
that ends within ``fail_time`` (2 s) of its launch is taken by libEnsemble's
MPIExecutor for a launch failure and relaunched by the executor itself (up
to 5 times, with the same output files), which abics_mlref does not count
as an attempt.
"""

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mock import calc_energy  # noqa: E402

from pymatgen.core import Structure  # noqa: E402


def sample_index(workdir):
    """``AL<n>/<r>/input<i>/baseinput<k>`` -> i"""
    name = os.path.basename(os.path.dirname(os.path.abspath(workdir)))
    return int(name[len("input"):])


def selected(spec, key, index):
    value = spec.get(key, [])
    return value == "all" or index in value


def main():
    workdir = sys.argv[1] if len(sys.argv) > 1 else "."
    try:
        with open(os.path.join(HERE, "faults.json")) as f:
            spec = json.load(f)
    except OSError:
        spec = {}

    count_path = os.path.join(workdir, "attempt.count")
    attempt = 1
    if os.path.exists(count_path):
        with open(count_path) as f:
            attempt = int(f.read().strip() or 0) + 1
    with open(count_path, "w") as f:
        f.write(f"{attempt}\n")

    index = sample_index(workdir)
    if selected(spec, "always_fail", index) or (attempt == 1 and selected(spec, "fail_first", index)):
        time.sleep(float(spec.get("fail_delay", 3.0)))
        print(f"mock_faulty: injected failure for sample {index}, attempt {attempt}", file=sys.stderr)
        sys.exit(1)
    if attempt == 1 and selected(spec, "hang_first", index):
        time.sleep(float(spec.get("hang_seconds", 60)))
    if spec.get("sleep", 0):
        time.sleep(float(spec["sleep"]))

    st = Structure.from_file(os.path.join(workdir, "structure.vasp"))
    with open(os.path.join(workdir, "energy.dat"), "w") as f:
        f.write("{:.15f}\n".format(calc_energy(st)))


if __name__ == "__main__":
    main()
