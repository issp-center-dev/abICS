from py_mc.mc import binning
import numpy as np
import sys
from mpi4py import MPI

nreplicas = int(sys.argv[1])
comm = MPI.COMM_WORLD
if comm.Get_size() != nreplicas:
    print("please run with the same number of processes as the number of replicas")
    sys.exit()


myreplica = comm.Get_rank()
obs = np.load(open(str(myreplica) + "/obs_save.npy", "rb"))
Trank_hist = np.load(open(str(myreplica) + "/Trank_hist.npy", "rb"))
kT_hist = np.load(open(str(myreplica) + "/kT_hist.npy", "rb"))
kTs = np.load("kTs.npy")
nsamples = obs.shape[0]
throwout = 100  # int(nsamples*0.5)
nobs = obs.shape[1]
# bin_level = int(np.log2(nsamples//30))
if myreplica == 0:
    rcvbuffer = np.empty(obs.shape[0])
    binning_fi = open("binning.dat", "w")
    obs_fi = open("obs_mean.dat", "w")


i_sample = 0
for i in range(nreplicas):  # we work on one T at a time
    # if myreplica == 0: print("working on Temp.  "+str(i)+"\n")
    if myreplica == 0:
        print(kTs[i], end="\t", file=obs_fi)
    for j in range(nobs):
        obs_T = np.where(Trank_hist == i, obs[:, j], 0.0)
        # print(obs_T[0], "test")
        if myreplica == 0:
            comm.Reduce(obs_T, rcvbuffer, op=MPI.SUM, root=0)

            print(rcvbuffer[throwout:].mean(), end="\t", file=obs_fi)
            # print(rcvbuffer)
            # bin_data = binning(rcvbuffer, bin_level)
            # print(np.max(bin_data), end='\t', file=obs_fi)
            # binning_fi.write("# "+str(i) + "," + str(j) + "\n")
            # binning_fi.write("\n".join([str(x) for x in bin_data])+"\n\n\n")

        else:
            comm.Reduce(obs_T, None, op=MPI.SUM, root=0)
    if myreplica == 0:
        print(file=obs_fi)
