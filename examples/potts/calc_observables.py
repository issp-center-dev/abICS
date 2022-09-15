import os
import numpy as np


def calc(ene, ms, T: float, nbootstrap: int = 10):
    N: int = ms.size
    nskip = N//2
    e2s = ene**2
    m2s = ms**2
    m4s = m2s**2

    obs = np.zeros((nbootstrap, 6))

    for ib in range(nbootstrap):
        idx = np.random.randint(nskip, N, size=N-nskip)
        e = np.mean(ene[idx])
        e2 = np.mean(e2s[idx])
        spec = (e2 - e * e) / (T * T)
        m = np.mean(ms[idx])
        m2 = np.mean(m2s[idx])
        m4 = np.mean(m4s[idx])
        u = m4 / (m2**2)
        # chi = (m2) / T
        chi = (m2 - m * m) / T
        obs[ib, :] = [e, spec, m, m2, chi, u]
    obs_mean = np.mean(obs, axis=0)
    obs_std = np.std(obs, axis=0)
    return obs_mean, obs_std


print("# 1: T")
for i, n in enumerate(["e", "spec", "m", "m2", "chi", "u"]):
    j = i+1
    print(f"# {2*j},{2*j+1}: {n}") 

Ts = np.load("kTs.npy")
for irep, T in enumerate(Ts):
    A = np.loadtxt(os.path.join("Tseparate", str(irep), "energies.dat"))
    obs_mean, obs_std = calc(A[:, 0], A[:, 1], T)
    print(f"{T}", end="")
    for i in range(len(obs_mean)):
        print(f" {obs_mean[i]} {obs_std[i]}", end="")
    print()
