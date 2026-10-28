import os
import numpy as np

obs = np.loadtxt(os.path.join("0", "obs.dat"))
num_Al = int(obs[-1, 4])
num_Mg = int(obs[-1, 5])

np.testing.assert_equal(num_Al, 0)
np.testing.assert_equal(num_Mg, 24)

# array-valued observable: columns 7 and 8 of obs.dat are num_Al_Mg[0], num_Al_Mg[1]
np.testing.assert_equal(obs[:, 6], obs[:, 4])
np.testing.assert_equal(obs[:, 7], obs[:, 5])
stats = np.loadtxt("num_Al_Mg.dat")
np.testing.assert_equal(stats.shape[1], 1 + 6 * 2)

dat = np.loadtxt("num_Al_Mg.dat")
npz = np.load("num_Al_Mg.npz")
np.testing.assert_allclose(npz["temperature"], dat[:, 0])
np.testing.assert_equal(npz["mean"].shape, (len(dat), 2))
for i in range(2):
    np.testing.assert_allclose(npz["mean"][:, i], dat[:, 1 + 6 * i])
    np.testing.assert_allclose(npz["second_moment"][:, i], dat[:, 3 + 6 * i])
    np.testing.assert_allclose(npz["fluctuation"][:, i], dat[:, 5 + 6 * i])
    np.testing.assert_allclose(npz["fluctuation_err"][:, i], dat[:, 6 + 6 * i])
np.testing.assert_equal(np.load("energy.npz")["mean"].ndim, 1)

print("OK")
