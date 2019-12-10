import numpy as np

nreplicas = 8
skip = 20

obsset = set()
datalist = [[] for i in range(nreplicas)]
for i in range(nreplicas):
    data = open(str(i) + "/obs.dat", "r")
    for _ in range(skip):
        next(data)
    for line in data:
        tmp = line.split()
        obsset.add(tmp[0])
        datalist[i].append([tmp[0], np.array([float(val) for val in tmp[1:]])])
    data.close()

obs = dict(zip(obsset, np.zeros(len(obsset))))
obs_nsample = dict(zip(obsset, np.zeros(len(obsset))))


for i in range(nreplicas):
    for j in range(len(datalist[i])):
        obs[datalist[i][j][0]] += datalist[i][j][1]
        obs_nsample[datalist[i][j][0]] += 1

# nsample = len(datalist[0])
# print(nsample)
# print(obs_nsample)


for kT in sorted(obs.keys()):
    print(kT, end="\t")
    print("\t".join([str(obs[kT][i] / obs_nsample[kT]) for i in range(len(obs[kT]))]))
