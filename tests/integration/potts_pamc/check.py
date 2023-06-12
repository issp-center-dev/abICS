import numpy as np

for name in ("energy", "magnetization", "logZ"):
    res = np.loadtxt(f"./{name}.dat")
    ref = np.loadtxt(f"./ref_{name}.dat")
    np.testing.assert_array_almost_equal(res, ref)
    print(f"{name} OK")
print("All OK")
