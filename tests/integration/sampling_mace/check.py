import os
import numpy as np

# The run must have produced the RXMC temperature file and per-replica observables.
assert os.path.exists("kTs.npy"), "kTs.npy was not produced"

obs = np.loadtxt(os.path.join("0", "obs.dat"))
obs = np.atleast_2d(obs)
assert obs.shape[0] >= 1, "no samples were recorded in 0/obs.dat"

# Column 3 is the total energy (see obs.dat header). For the pretrained CHGNet
# smoke test we do not hard-code values (they depend on the model version); we
# only check the pipeline produced finite, physically plausible energies.
energy = obs[:, 3]
assert np.all(np.isfinite(energy)), "non-finite energies in obs.dat"
assert np.all(energy < 0.0), "expected negative total energies for MgAl2O4"

print("OK")
