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

import itertools


def list_append(lst, item):
    lst.append(item)
    return lst


def list_extend(lst, item):
    lst.extend(item)
    return lst


def make_ising_list(n, latgas_vals=[-1, 1]):
    if n == 1:
        return [[val] for val in latgas_vals]
    lst = []
    for val in latgas_vals:
        lst0 = make_ising_list(n - 1)
        lst0 = [list_append(rep, val) for rep in lst0]
        lst.extend(lst0)
    return lst


def make_sumlist(n, val):
    # make list of n integers that sum up to val
    # ordering is distinguished here
    if n == 1:
        return [[val]]
    lst = []
    for i in range(val + 1):
        lst0 = [list_extend([i], sl) for sl in make_sumlist(n - 1, val - i)]
        lst.extend(lst0)
    return lst


# print(make_sumlist(3,5))


def make_numberlist(spec_id=[0, 1, 2], nspec=[1, 4, 1], spec_orient=[1, 2, 1]):
    # make list for the number of each specie-orientation in the system
    # from the number of each species and possible number of orientations
    # also returns accompanying list of specie-orientation designators
    lst0 = [[]]
    spec_orient_id = []
    for i in range(len(spec_id)):
        if spec_orient[i] > 1:
            # calculate all possible orientation combinations
            for j in range(spec_orient[i]):
                spec_orient_id.append((spec_id[i], j))
            tmp = []
            for sl in make_sumlist(spec_orient[i], nspec[i]):
                for ls in lst0:
                    tmp.append(list_extend(ls.copy(), sl))
            lst0 = tmp

        else:
            spec_orient_id.append(spec_id[i])
            lst0 = [list_append(ls, nspec[i]) for ls in lst0]
            # print(lst0)
    return lst0, spec_orient_id


# make_list2()
# print(make_numberlist(spec_orient=[1,2,2]))


def make_latgas_reps(spec_id=[0, 1, 2], nspec=[1, 4, 1], spec_orient=[1, 2, 1]):
    numlist, spec_orient_id = make_numberlist(spec_id, nspec, spec_orient)
    reps = []
    tmp = []
    for numls in numlist:
        for i in range(len(spec_orient_id)):
            tmp += [spec_orient_id[i]] * numls[i]
        reps.extend(list(set(tuple(itertools.permutations(tmp)))))
        tmp = []
    return reps


# print(make_latgas_reps(spec_id=[0], nspec=[6], spec_orient=[2]))
# print(len(make_latgas_reps(spec_id=[0], nspec=[6], spec_orient=[2])))
"""
def make_latgas_reps(spec_id=[0,1,2], nspec=[1,4,1], spec_orient=[1,2,1]):
    reps = []
    # first make permutations disregarding orientation
    for i in range(len(spec_id)):
        reps += [spec_id[i]] * nspec[i]
    reps = list(set(tuple(itertools.permutations(reps))))
    tmp = []
    # then consider orientation for each species
    for rep in reps:
        rep_list = list(rep)
        for ispec in range(len(spec_id)):
            if spec_orient[ispec] > 1:
                ising_list = make_ising_list(nspec[ispec], range(spec_orient[ispec]))
                for i in range(nspec[ispec]):
                    id = rep_list.index(spec_id[ispec])
                    rep_list
                

    return rep

print(make_latgas_reps())
print(len(make_latgas_reps()))
"""
# print(list(itertools.permutations([1,1,1,1])))
