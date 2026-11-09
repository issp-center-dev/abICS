#!/bin/sh
# Same loop as run.sh, with abics_mlref in integrated mode (requires libEnsemble).
# The first active-learning step is also run in the staged mode with the same
# [mlref] seed, and its training data must be identical to the integrated one.

export OMP_NUM_THREADS=1

set -e

cleanup() {
  rm -f ./ALloop.progress ./rundirs.txt ./runtask.log
  rm -rf ./AL0 ./AL1 ./MC0 ./MC1 ./aenetXSF ./train0 ./generate0 ./baseinput
}

cleanup
rm -rf ./ref_staged_AL0

# --- staged AL0 as the reference (abics_mlref -> mock solver -> abics_mlref) ---
echo start staged AL0 reference
mpiexec -np 2 --oversubscribe abics_mlref input.toml
sh parallel_run.sh
mpiexec -np 2 --oversubscribe abics_mlref input.toml
mv AL0 ref_staged_AL0
cleanup

# --- the integrated loop ---
sh ./AL_integrated.sh
echo compare integrated AL0 with the staged reference
python3 - <<'PYEOF'
import filecmp, glob, os, sys
bad = []
for ref in sorted(glob.glob("ref_staged_AL0/*/structure.*.vasp") + glob.glob("ref_staged_AL0/*/energy_corr.dat")):
    new = ref.replace("ref_staged_AL0", "AL0", 1)
    if not filecmp.cmp(ref, new, shallow=False):
        bad.append(new)
if bad:
    print("DIFFERENT from the staged mode:", *bad, sep="\n  ")
    sys.exit(1)
print("integrated AL0 identical to the staged AL0 (structures and energy_corr.dat)")
PYEOF
sh ./MC.sh
sh ./AL_integrated.sh
sh ./MC.sh

# Verify the outputs of AL0 -> train -> MC0 -> AL1 -> train -> MC1
python3 ../check_active_learn.py input_integrated.toml 2
echo OK
