import numpy as np

result = np.loadtxt("./result.dat")
logZ = np.loadtxt("./logZ.dat")

ref_result = np.loadtxt("./ref_result.dat")
ref_logZ = np.loadtxt("./ref_logZ.dat")

np.testing.assert_array_almost_equal(result, ref_result)
np.testing.assert_array_almost_equal(logZ, ref_logZ)

print("OK")
