# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

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
print (obsset)
obs = dict(zip(obsset, " "*(len(obsset))))
print (obs.keys())
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
