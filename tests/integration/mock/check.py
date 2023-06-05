import os
import numpy as np

obs = np.loadtxt(os.path.join("0", "obs.dat"))
num_Al = int(obs[-1, 4])
num_Mg = int(obs[-1, 5])

np.testing.assert_equal(num_Al, 0)
np.testing.assert_equal(num_Mg, 24)

print("OK")
