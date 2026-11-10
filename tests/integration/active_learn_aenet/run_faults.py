#!/usr/bin/env python3
"""
Fault-injection test of the integrated ``abics_mlref`` and the ``abics_al``
driver (requires libEnsemble; run through ``run_faults.sh`` in this
directory).

The reference solver is ``mock_faulty.py``, which fails, hangs or sleeps
as told by ``faults.json``; the driver and its stages are real processes.
Each scenario starts from a clean directory and checks the files that the
programs leave behind (``abics_mlref.done`` markers with their attempt
counts, ``ALloop.progress``, the lock file, moved-aside ``MC<n>``) and the
exit status:

1. solver runs that fail once are retried and the result equals a
   fault-free run byte for byte;
2. solver runs that exceed ``[mlref.integrated] timeout`` are killed and retried;
3. solver runs killed from outside (SIGKILL) are retried;
4. a run that always fails stops the driver with status 1, a missing input
   directory is detected on the rerun, and after the fault is removed the
   rerun does not repeat the finished runs;
5. the driver killed (SIGTERM) during the reference-data step leaves no
   lock and no processes, and the rerun keeps the finished runs;
6. the driver killed during sampling: the rerun refuses while ``MC<n>`` is
   recent, then moves it aside and repeats the sampling;
7. the wall time: no progress is an error without resubmission; a stop
   after a finished stage exits with status 2 and runs ``resubmit_command``;
8. ``ALloop.progress`` and the directories disagree: the driver refuses;
9. the driver killed with SIGKILL during sampling: the sampling stage
   survives, the rerun refuses while it runs (``al_driver/stage.pid``) and
   continues once it has ended.
"""

import contextlib
import filecmp
import glob
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time

import toml

HERE = os.path.dirname(os.path.abspath(__file__))
BASE_TOML = "input_driver.toml"
FAULTS_JSON = "faults.json"
LOG_DIR = "faults_logs"
RESUBMITTED = "faults_resubmitted"
NDATA = 10
NREPLICAS = 2

GENERATED_FILES = ["ALloop.progress", "rundirs.txt", "runtask.log", "al_summary.dat", RESUBMITTED]
GENERATED_DIRS = ["AL0", "AL1", "AL2", "MC0", "MC1", "MC2", "aenetXSF", "train0", "generate0",
                  "baseinput", "al_driver", "faults_aside"]
REFERENCE_DIR = "faults_ref_AL0"
"""The fault-free AL0 of scenario 0, kept until the end for the comparisons."""


class Failure(Exception):
    pass


def check(cond, msg):
    if not cond:
        raise Failure(msg)


def banner(text):
    print(f"\n##### {text}", flush=True)


# ---------------------------------------------------------------------------
# set-up helpers
# ---------------------------------------------------------------------------


def cleanup(reference=False):
    for name in GENERATED_FILES:
        if os.path.exists(name):
            os.remove(name)
    for name in GENERATED_DIRS + ([REFERENCE_DIR] if reference else []):
        shutil.rmtree(name, ignore_errors=True)


def set_faults(**spec):
    with open(FAULTS_JSON, "w") as f:
        json.dump(spec, f)


def write_input(name, timeout=None, nsteps=None, **al):
    """A variant of input_driver.toml with mock_faulty.py as the solver."""
    params = toml.load(BASE_TOML)
    params["mlref"]["solver"]["path"] = "./mock_faulty.py"
    if timeout is not None:
        params["mlref"]["integrated"]["timeout"] = timeout
    if nsteps is not None:
        params["sampling"]["nsteps"] = nsteps
    params["al"].update(al)
    path = f"faults_{name}.toml"
    with open(path, "w") as f:
        toml.dump(params, f)
    return path


# ---------------------------------------------------------------------------
# processes
# ---------------------------------------------------------------------------


def run(args, expect, timeout=900, logname=None):
    """Run a command to completion; its output goes to faults_logs/<logname>.log and is returned."""
    logname = logname or "_".join(args)
    logfile = os.path.join(LOG_DIR, re.sub(r"[^\w.-]+", "_", logname) + ".log")
    print(f"$ {' '.join(args)}  (expect status {expect}; log {logfile})", flush=True)
    with open(logfile, "w") as fo:
        proc = subprocess.run(args, stdout=fo, stderr=subprocess.STDOUT, timeout=timeout)
    with open(logfile) as f:
        output = f.read()
    if proc.returncode != expect:
        print(output[-4000:])
        raise Failure(f"{' '.join(args)}: exit status {proc.returncode}, expected {expect}")
    return output


@contextlib.contextmanager
def background(args, logname):
    """Run a command in the background, in a session of its own.

    The command and its direct children that do not change their group
    (the libEnsemble workers of abics_al, mpirun) share its process group,
    so :func:`group_alive` sees them whatever their command lines look
    like, and the whole group is killed if the scenario fails before the
    command has ended.  Open MPI gives every rank a group of its own and
    abics_al starts its training and sampling stages in sessions of their
    own, so those are checked by name (:func:`processes`).
    """
    logfile = os.path.join(LOG_DIR, logname + ".log")
    print(f"$ {' '.join(args)} &  (log {logfile})", flush=True)
    fo = open(logfile, "w")
    proc = subprocess.Popen(args, stdout=fo, stderr=subprocess.STDOUT, start_new_session=True)
    proc.logfile = logfile
    try:
        yield proc
    finally:
        if proc.poll() is None or group_alive(proc):
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait()
        fo.close()


def group_alive(proc):
    """True while any process of the group started by :func:`background` is alive."""
    return bool(subprocess.run(["pgrep", "-g", str(proc.pid)], capture_output=True, text=True).stdout.split())


def finish(proc, expect, timeout=900):
    try:
        rc = proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        raise Failure(f"background {' '.join(proc.args)} did not finish within {timeout} s")
    with open(proc.logfile) as f:
        output = f.read()
    if rc != expect:
        print(output[-4000:])
        raise Failure(f"background {' '.join(proc.args)}: exit status {rc}, expected {expect}")
    return output


def wait_for_group_exit(proc, timeout, what):
    t = wait_until(lambda: not group_alive(proc), timeout, f"the process group of {what} has ended")
    print(f"  the process group of {what} (libEnsemble workers, launchers) ended within {t:.1f} s")


def wait_until(cond, timeout, what):
    t0 = time.time()
    while not cond():
        if time.time() - t0 > timeout:
            raise Failure(f"timed out after {timeout} s waiting until {what}")
        time.sleep(0.2)
    return time.time() - t0


def processes(pattern):
    """``(pid, command line)`` of the live processes whose command line contains ``pattern``."""
    pids = subprocess.run(["pgrep", "-f", pattern], capture_output=True, text=True).stdout.split()
    result = []
    for pid in pids:
        pid = int(pid)
        if pid == os.getpid():
            continue
        args = subprocess.run(["ps", "-o", "args=", "-p", str(pid)], capture_output=True, text=True).stdout.strip()
        if args:
            result.append((pid, args))
    return result


def solver_processes():
    """The python processes running mock_faulty.py (not the mpirun that started them)."""
    return [(pid, args) for pid, args in processes("mock_faulty.py")
            if "python" in os.path.basename(args.split()[0]).lower()]


# ---------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------


def done_markers(al="AL0"):
    """``{baseinput0 directory: (mtime, attempts)}`` of the finished solver runs of ``al``."""
    markers = {}
    for path in glob.glob(os.path.join(al, "*", "input*", "baseinput0", "abics_mlref.done")):
        with open(path) as f:
            m = re.search(r"after (\d+) attempt", f.read())
        if m is None:
            continue  # being written right now: not finished yet
        markers[os.path.dirname(path)] = (os.stat(path).st_mtime_ns, int(m.group(1)))
    return markers


def attempt_count(stage_dir):
    with open(os.path.join(stage_dir, "attempt.count")) as f:
        return int(f.read())


def sample_of(stage_dir):
    return int(os.path.basename(os.path.dirname(stage_dir))[len("input"):])


def stage_dirs(indices, al="AL0"):
    """The ``baseinput0`` directories of the given sample indices in every replica."""
    return {os.path.join(al, str(r), f"input{i:02d}", "baseinput0") for r in range(NREPLICAS) for i in indices}


def check_attempts(markers, retried):
    """Every run in ``retried`` (a set of baseinput0 directories) took 2 attempts, all others 1."""
    check(len(markers) == NDATA * NREPLICAS, f"{len(markers)} finished runs, expected {NDATA * NREPLICAS}")
    check(retried <= set(markers), f"retried runs without a marker: {sorted(retried - set(markers))}")
    for stage_dir, (_, attempts) in sorted(markers.items()):
        expected = 2 if stage_dir in retried else 1
        check(attempts == expected, f"{stage_dir}: {attempts} attempt(s) recorded, expected {expected}")
        check(attempt_count(stage_dir) == expected,
              f"{stage_dir}: mock_faulty.py was started {attempt_count(stage_dir)} times, expected {expected}")
    print(f"  {len(markers)} runs finished; retried once: {sorted(retried)}")


def check_unchanged(before, after):
    """The runs finished before an interruption were not run again."""
    for stage_dir, (mtime, attempts) in before.items():
        check(stage_dir in after, f"{stage_dir}: marker disappeared")
        check(after[stage_dir] == (mtime, attempts), f"{stage_dir}: the finished run was repeated")
    print(f"  {len(before)} runs finished before the interruption were not repeated")


def progress():
    if not os.path.exists("ALloop.progress"):
        return []
    with open("ALloop.progress") as f:
        return [line.strip() for line in f if line.strip()]


def check_outputs(tomlfile, nloops):
    run(["python3", "../check_active_learn.py", tomlfile, str(nloops)], expect=0,
        logname=f"check_{tomlfile}_{nloops}")
    print(f"  check_active_learn.py: {nloops} loop(s) OK")


def compare_with_reference(ref=REFERENCE_DIR, new="AL0"):
    files = sorted(glob.glob(os.path.join(ref, "*", "structure.*.vasp")) + glob.glob(os.path.join(ref, "*", "energy_corr.dat")))
    check(files, f"no reference files under {ref}")
    for f in files:
        g = f.replace(ref, new, 1)
        check(os.path.exists(g), f"{g} missing")
        check(filecmp.cmp(f, g, shallow=False), f"{g} differs from the fault-free run {f}")
    print(f"  {len(files)} files identical to the fault-free run")


def age(path, seconds=300.0):
    """Make ``path`` and everything below it look ``seconds`` old."""
    t = time.time() - seconds
    os.utime(path, (t, t))
    for root, dirs, files in os.walk(path):
        for name in dirs + files:
            os.utime(os.path.join(root, name), (t, t))


# ---------------------------------------------------------------------------
# scenarios
# ---------------------------------------------------------------------------


def scenario_reference(tomlfile):
    banner("0. fault-free AL0 as the reference")
    cleanup(reference=True)
    set_faults()
    run(["abics_mlref", tomlfile], expect=0, logname="s0_mlref")
    check_attempts(done_markers(), retried=set())
    shutil.move("AL0", REFERENCE_DIR)
    print("  reference finished in 1 attempt each")


def scenario_fail_once(tomlfile):
    banner("1. solver runs failing once are retried")
    cleanup()
    retried = [0, 3, 7]
    set_faults(fail_first=retried)
    run(["abics_mlref", tomlfile], expect=0, logname="s1_mlref")
    check(progress() == ["AL0"], f"ALloop.progress = {progress()}")
    check_attempts(done_markers(), stage_dirs(retried))
    compare_with_reference()


def scenario_timeout():
    banner("2. solver runs exceeding the timeout are killed and retried")
    cleanup()
    # the timeout must leave a healthy run (python start, pymatgen import) enough
    # time on a slow CI runner, and the hang must be long enough to tell apart
    timeout, hang = 20, 600
    tomlfile = write_input("timeout", timeout=timeout)
    retried = [2, 5]
    set_faults(hang_first=retried, hang_seconds=hang)
    t0 = time.time()
    run(["abics_mlref", tomlfile], expect=0, logname="s2_mlref")
    elapsed = time.time() - t0
    check(elapsed < hang, f"abics_mlref took {elapsed:.0f} s: the hanging runs were not killed by the timeout")
    check_attempts(done_markers(), stage_dirs(retried))
    compare_with_reference()
    print(f"  finished in {elapsed:.0f} s although {NREPLICAS * len(retried)} runs would have hung for {hang} s")


def scenario_external_kill(tomlfile):
    banner("3. solver runs killed from outside are retried")
    cleanup()
    set_faults(sleep=8)
    with background(["abics_al", tomlfile], "s3_driver") as proc:
        wait_until(lambda: len(solver_processes()) >= 2 or proc.poll() is not None, 180, "2 solver runs are active")
        check(proc.poll() is None, "abics_al ended before any solver run could be killed")
        # A task that dies within 2 s of its launch is relaunched by libEnsemble's executor
        # itself (see mock_faulty.py): kill only runs that have been alive for longer.
        victims = solver_processes()
        time.sleep(3)
        alive = dict(solver_processes())
        killed = set()
        for pid, args in victims:
            if alive.get(pid) != args:
                continue
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                continue
            killed.add(os.path.relpath(args.split()[-1], HERE))
        check(killed, "no solver run could be killed")
        print(f"  killed the solver runs in {sorted(killed)}")
        finish(proc, expect=0)
    check(progress() == ["AL0", "train", "MC0"], f"ALloop.progress = {progress()}")
    check_attempts(done_markers(), killed)
    check_outputs(tomlfile, 1)


def scenario_permanent_failure(tomlfile):
    banner("4. a run that always fails, a missing input directory, and the rerun")
    cleanup()
    set_faults(always_fail=[5])
    output = run(["abics_al", tomlfile], expect=1, logname="s4_driver_fails")
    check(not os.path.exists("ALloop.progress"), "ALloop.progress was written although solver runs failed")
    check(not os.path.exists(os.path.join("al_driver", "abics_al.lock")), "the lock file was left behind")
    check("2 of 20 solver runs failed" in output, "the failed runs were not counted in the log")
    check(re.search(r"AL0/0/input05: stage 0: solver exited with status 1 \(after 2 attempt", output)
          and re.search(r"AL0/1/input05: stage 0: solver exited", output), "the failed runs were not listed")
    before = done_markers()
    check_attempts_partial(before, failed={5})

    banner("4b. an input directory removed by hand is detected")
    os.makedirs("faults_aside")
    shutil.move(os.path.join("AL0", "0", "input05"), "faults_aside")
    output = run(["abics_al", tomlfile], expect=1, logname="s4_driver_missing_input")
    check("has 9 input directories, expected ndata = 10" in output, "the missing input directory was not reported")
    check("cannot be resumed" in output, "the resume was not refused")
    shutil.move(os.path.join("faults_aside", "input05"), os.path.join("AL0", "0"))

    banner("4c. the fault removed: the rerun finishes without repeating the finished runs")
    set_faults()
    run(["abics_al", tomlfile], expect=0, logname="s4_driver_rerun")
    after = done_markers()
    check_unchanged(before, after)
    check(len(after) == NDATA * NREPLICAS, f"{len(after)} finished runs, expected {NDATA * NREPLICAS}")
    for r in range(NREPLICAS):
        stage_dir = os.path.join("AL0", str(r), "input05", "baseinput0")
        check(attempt_count(stage_dir) == 3, f"{stage_dir}: started {attempt_count(stage_dir)} times, expected 3")
    check_outputs(tomlfile, 1)
    compare_with_reference()


def check_attempts_partial(markers, failed):
    expected = NDATA * NREPLICAS - NREPLICAS * len(failed)
    check(len(markers) == expected, f"{len(markers)} finished runs, expected {expected}")
    for stage_dir, (_, attempts) in markers.items():
        check(sample_of(stage_dir) not in failed, f"{stage_dir}: a failing run has a done marker")
        check(attempts == 1, f"{stage_dir}: {attempts} attempt(s), expected 1")
    for r in range(NREPLICAS):
        for i in failed:
            stage_dir = os.path.join("AL0", str(r), f"input{i:02d}", "baseinput0")
            check(attempt_count(stage_dir) == 2, f"{stage_dir}: started {attempt_count(stage_dir)} times, expected 2")
    print(f"  {len(markers)} runs finished, samples {sorted(failed)} failed after 2 attempts")


def scenario_driver_killed_in_mlref(tomlfile):
    banner("5. the driver killed during the reference-data step")
    cleanup()
    set_faults(sleep=3)
    with background(["abics_al", tomlfile], "s5_driver_killed") as proc:
        wait_until(lambda: len(done_markers()) >= 4 or proc.poll() is not None, 300, "4 solver runs have finished")
        check(proc.poll() is None, "abics_al ended before it could be killed")
        before = done_markers()
        proc.send_signal(signal.SIGTERM)
        output = finish(proc, expect=1, timeout=120)
        check("interrupted by SIGTERM" in output, "the interruption was not reported")
        check(not os.path.exists(os.path.join("al_driver", "abics_al.lock")), "the lock file was left behind")
        check(not os.path.exists("ALloop.progress"), "ALloop.progress was written by the killed run")
        wait_for_group_exit(proc, 30, "the killed abics_al")
        check(not solver_processes(), "solver runs are still alive")

    banner("5b. the rerun keeps the finished runs and completes the loop")
    set_faults()
    run(["abics_al", tomlfile], expect=0, logname="s5_driver_rerun")
    check_unchanged(before, done_markers())
    check(progress() == ["AL0", "train", "MC0", "AL1", "train", "MC1"], f"ALloop.progress = {progress()}")
    check_outputs(tomlfile, 2)


def scenario_driver_killed_in_sampling():
    banner("6. the driver killed during sampling")
    cleanup()
    set_faults()
    tomlfile = write_input("sampling", nsteps=640)
    with background(["abics_al", tomlfile], "s6_driver_killed") as proc:
        wait_until(lambda: os.path.isdir("MC0") or proc.poll() is not None, 300, "MC0 is created")
        check(proc.poll() is None, "abics_al ended before it could be killed")
        check(progress() == ["AL0", "train"], f"ALloop.progress = {progress()} when MC0 appeared")
        proc.send_signal(signal.SIGTERM)
        output = finish(proc, expect=1, timeout=120)
        check("interrupted by SIGTERM" in output and "terminating the stage" in output,
              "the interruption of the sampling stage was not reported")
        wait_for_group_exit(proc, 60, "the killed abics_al")
        check(not processes("abics_sampling"), "abics_sampling is still alive")
    check(progress() == ["AL0", "train"], f"ALloop.progress = {progress()} after the kill")
    check(os.path.isdir("MC0"), "MC0 disappeared")

    banner("6b. the rerun refuses while MC0 looks alive")
    output = run(["abics_al", tomlfile], expect=1, logname="s6_driver_refuses")
    check("a sampling run may still be active" in output, "the recent MC0 was not refused")
    check(os.path.isdir("MC0") and not glob.glob("al_driver/MC0.interrupted-*"), "MC0 was moved although recent")

    banner("6c. an old unfinished MC0 is moved aside and the sampling repeated")
    age("MC0")
    output = run(["abics_al", tomlfile], expect=0, logname="s6_driver_rerun")
    check("the sampling run did not finish" in output, "the interrupted sampling was not reported")
    moved = glob.glob("al_driver/MC0.interrupted-*")
    check(len(moved) == 1, f"moved-aside directories: {moved}")
    check(progress() == ["AL0", "train", "MC0", "AL1", "train", "MC1"], f"ALloop.progress = {progress()}")
    check_outputs(tomlfile, 2)
    print(f"  interrupted sampling kept in {moved[0]}")


def scenario_walltime():
    banner("7. wall time: no progress is an error, no resubmission")
    cleanup()
    set_faults()
    tomlfile = write_input("walltime_noprogress", max_iterations=1, walltime=30, walltime_margin=5,
                           stage_estimates={"mlref": 1000}, resubmit_command=f"touch {RESUBMITTED}")
    output = run(["abics_al", tomlfile], expect=1, logname="s7_noprogress")
    check("Not resubmitting" in output, "the no-progress case was not explained")
    check(not os.path.exists(RESUBMITTED), "the resubmission command ran although nothing was done")
    check(not os.path.isdir("AL0"), "AL0 was started although it did not fit")

    banner("7b. stop after a finished stage: status 2 and the resubmission command")
    tomlfile = write_input("walltime", max_iterations=1, walltime=30, walltime_margin=5,
                           stage_estimates={"mlref": 1, "train": 1000}, resubmit_command=f"touch {RESUBMITTED}")
    output = run(["abics_al", tomlfile], expect=2, logname="s7_walltime")
    check("Stopping before the train stage of step 0" in output, "the stop was not reported")
    check(os.path.exists(RESUBMITTED), "the resubmission command did not run")
    check(progress() == ["AL0"], f"ALloop.progress = {progress()}")
    with open(os.path.join("al_driver", "durations.json")) as f:
        durations = json.load(f)
    check("mlref" in durations and durations["mlref"] > 0, f"durations.json = {durations}")
    check(not os.path.exists(os.path.join("al_driver", "abics_al.lock")), "the lock file was left behind")

    banner("7c. the resubmitted run (no wall time limit) finishes the loop")
    tomlfile = write_input("walltime_rerun", max_iterations=1)
    run(["abics_al", tomlfile], expect=0, logname="s7_rerun")
    check(progress() == ["AL0", "train", "MC0"], f"ALloop.progress = {progress()}")
    check(os.path.exists(os.path.join("al_driver", "model_metadata.json")), "model_metadata.json missing")
    check_outputs(tomlfile, 1)
    return tomlfile


def scenario_inconsistent_state(tomlfile):
    banner("8. ALloop.progress and the directories disagree")
    shutil.rmtree("MC0")
    output = run(["abics_al", tomlfile], expect=1, logname="s8_inconsistent")
    check("ALloop.progress ends with MC0 but 0 MC run(s) exist" in output, "the inconsistency was not reported")
    check(not os.path.exists(os.path.join("al_driver", "abics_al.lock")), "the lock file was left behind")


def scenario_driver_sigkilled_in_sampling():
    banner("9. the driver killed with SIGKILL during sampling: the orphaned stage is detected")
    cleanup()
    set_faults()
    tomlfile = write_input("sigkill", nsteps=1280)
    # this scenario's sampling stage (the driver passes the input file as an absolute path), not another run's
    sampling = f"abics_sampling .*{re.escape(tomlfile)}"
    stage_record = os.path.join("al_driver", "stage.pid")
    try:
        with background(["abics_al", tomlfile], "s9_driver_sigkilled") as proc:
            wait_until(lambda: os.path.isdir("MC0") or proc.poll() is not None, 300, "MC0 is created")
            check(proc.poll() is None, "abics_al ended before it could be killed")
            proc.send_signal(signal.SIGKILL)
            finish(proc, expect=-signal.SIGKILL, timeout=60)
            # the driver could not stop the stage: the lock and the stage record stay, the sampling goes on
            check(os.path.exists(os.path.join("al_driver", "abics_al.lock")), "the lock is missing after SIGKILL")
            check(os.path.exists(stage_record), "stage.pid is missing after SIGKILL")
            check(processes(sampling), "the sampling stage did not survive the SIGKILL of the driver")

            banner("9b. the rerun refuses while the orphaned sampling runs")
            output = run(["abics_al", tomlfile], expect=1, logname="s9_driver_refuses")
            check("Removing the stale lock file" in output, "the dead driver's lock was not recognized as stale")
            check("is still running (process group" in output, "the orphaned sampling stage was not detected")
            check(not os.path.exists(os.path.join("al_driver", "abics_al.lock")), "the lock was left after the refusal")
            check(os.path.exists(stage_record), "stage.pid was removed although the stage runs")

            banner("9c. the orphaned sampling ends; the rerun removes the stale record and continues")
            wait_until(lambda: not processes(sampling), 300, "the orphaned sampling has ended")
    finally:
        # on a failure the orphaned stage (in a session of its own) must not outlive the scenario
        if os.path.exists(stage_record):
            with open(stage_record) as f:
                fields = f.read().split()
            if len(fields) > 1 and fields[1].isdigit():
                try:
                    os.killpg(int(fields[1]), signal.SIGKILL)
                except (ProcessLookupError, PermissionError):
                    pass
    check(progress() == ["AL0", "train", "MC0"], f"ALloop.progress = {progress()} after the orphan ended")
    output = run(["abics_al", tomlfile], expect=0, logname="s9_driver_rerun")
    check("Removing the stale" in output and "stage.pid" in output, "the stale stage record was not reported")
    check(not os.path.exists(os.path.join("al_driver", "stage.pid")), "stage.pid was left behind")
    check(progress() == ["AL0", "train", "MC0", "AL1", "train", "MC1"], f"ALloop.progress = {progress()}")
    check_outputs(tomlfile, 2)


def main():
    """``run_faults.py [scenario ...]``: all scenarios, or only the numbered ones (1-9)."""
    os.chdir(HERE)
    shutil.rmtree(LOG_DIR, ignore_errors=True)
    os.makedirs(LOG_DIR)
    selected = [int(a) for a in sys.argv[1:]] or list(range(1, 10))
    t0 = time.time()
    try:
        tomlfile = write_input("driver1", max_iterations=1)
        tomlfile2 = write_input("driver2", max_iterations=2)
        if {1, 2, 4} & set(selected):
            scenario_reference(tomlfile)
        for number in selected:
            if number == 1:
                scenario_fail_once(tomlfile)
            elif number == 2:
                scenario_timeout()
            elif number == 3:
                scenario_external_kill(tomlfile)
            elif number == 4:
                scenario_permanent_failure(tomlfile)
            elif number == 5:
                scenario_driver_killed_in_mlref(tomlfile2)
            elif number == 6:
                scenario_driver_killed_in_sampling()
            elif number == 7:
                last = scenario_walltime()
            elif number == 8:
                if 7 not in selected:
                    last = scenario_walltime()
                scenario_inconsistent_state(last)
            elif number == 9:
                scenario_driver_sigkilled_in_sampling()
            else:
                raise Failure(f"unknown scenario {number}")
    except Failure as e:
        print(f"\nFAULT TEST FAILED: {e}  (outputs kept in {HERE}, logs in {LOG_DIR}/)")
        sys.exit(1)
    cleanup(reference=True)
    for path in glob.glob("faults_*.toml") + [FAULTS_JSON]:
        os.remove(path)
    shutil.rmtree(LOG_DIR, ignore_errors=True)
    print(f"\nFAULT TEST OK: scenarios {selected} in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
