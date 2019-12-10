import numpy as np
nreplicas = 8
skip = 0

obsset = set()
datalist = [[] for i in range(nreplicas)]
for i in range(nreplicas):
    data = open(str(i)+"/obs.dat", 'r')
    for _ in range(skip):
        next(data)
    for line in data:
        tmp = line.split()
        obsset.add(tmp[0])
        datalist[i].append([tmp[0], np.array([float(val) for val in tmp[1:]])])
    data.close()
print obsset
obs = dict(zip(obsset, " "*(len(obsset))))
print obs.keys()
obs_nsample = dict(zip(obsset, np.zeros(len(obsset))))


for j in range(len(datalist[0])):
    for i in range(nreplicas):
        obs[datalist[i][j][0]] += "\t".join([str(datalist[i][j][1][k])
                                             for k in range(len(datalist[i][j][1]))])+"\n"
        #obs_nsample[datalist[i][j][0]] += 1

#nsample = len(datalist[0])
# print(nsample)
# print(obs_nsample)


for kT in sorted(obs.keys()):
    with open(str(kT)+"_converge.dat", "w") as f:
        f.write(obs[kT])
