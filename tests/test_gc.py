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

# ----------
# Unit test for GrandCanonical features
#   pytest -v -s --log-cli-level=DEBUG tests/test_gc.py
# ----------

import unittest

import os,sys
import numpy as np
import numpy.random as random
import toml
import copy

from abics.applications.latgas_abinitio_interface.defect import DFTConfigParams, defect_config
from abics.applications.latgas_abinitio_interface.model_setup import DFTLatticeGas, Config, match_latgas_group

random_seed = 123456789

import logging
logger = logging.getLogger(__name__)


class StabCalculator:
    def __init__(self):
        pass
    def submit(self, structure, path):
        print("StabCalculator.submit")
        e0 = 0.0
        return e0, structure

class Params:
    def __init__(self):
        pass
    def create_params(self, sublat, chem=None):

        params = toml.loads(
            """
[config]
unitcell = [[8.1135997772, 0.0000000000, 0.0000000000],
            [0.0000000000, 8.1135997772, 0.0000000000],
            [0.0000000000, 0.0000000000, 8.1135997772]]
supercell = [1,1,1]

[[config.base_structure]]
type = "O"
coords = [
     [0.237399980, 0.237399980, 0.237399980],
     [0.762599945, 0.762599945, 0.762599945],
     [0.512599945, 0.012600004, 0.737399936],
     [0.487399966, 0.987399936, 0.262599975],
     [0.012600004, 0.737399936, 0.512599945],
     [0.987399936, 0.262599975, 0.487399966],
     [0.737399936, 0.512599945, 0.012600004],
     [0.262599975, 0.487399966, 0.987399936],
     [0.987399936, 0.487399966, 0.262599975],
     [0.012600004, 0.512599945, 0.737399936],
     [0.487399966, 0.262599975, 0.987399936],
     [0.512599945, 0.737399936, 0.012600004],
     [0.262599975, 0.987399936, 0.487399966],
     [0.737399936, 0.012600004, 0.512599945],
     [0.237399980, 0.737399936, 0.737399936],
     [0.762599945, 0.262599975, 0.262599975],
     [0.512599945, 0.512599945, 0.237399980],
     [0.487399966, 0.487399966, 0.762599945],
     [0.012600004, 0.237399980, 0.012600004],
     [0.987399936, 0.762599945, 0.987399936],
     [0.987399936, 0.987399936, 0.762599945],
     [0.012600004, 0.012600004, 0.237399980],
     [0.487399966, 0.762599945, 0.487399966],
     [0.512599945, 0.237399980, 0.512599945],
     [0.737399936, 0.237399980, 0.737399936],
     [0.262599975, 0.762599945, 0.262599975],
     [0.237399980, 0.512599945, 0.512599945],
     [0.762599945, 0.487399966, 0.487399966],
     [0.762599945, 0.987399936, 0.987399936],
     [0.237399980, 0.012600004, 0.012600004],
     [0.737399936, 0.737399936, 0.237399980],
     [0.262599975, 0.262599975, 0.762599945],
     ]
            """
        )

        _tbl = toml.loads(
            """
# [[config.defect_structure]]
coords = [
     [0.000000000, 0.000000000, 0.000000000],
     [0.749999940, 0.249999985, 0.499999970],
     [0.249999985, 0.749999940, 0.499999970],
     [0.249999985, 0.499999970, 0.749999940],
     [0.749999940, 0.499999970, 0.249999985],
     [0.499999970, 0.749999940, 0.249999985],
     [0.499999970, 0.249999985, 0.749999940],
     [0.000000000, 0.499999970, 0.499999970],
     [0.749999940, 0.749999940, 0.000000000],
     [0.249999985, 0.249999985, 0.000000000],
     [0.249999985, 0.000000000, 0.249999985],
     [0.749999940, 0.000000000, 0.749999940],
     [0.499999970, 0.000000000, 0.499999970],
     [0.000000000, 0.749999940, 0.749999940],
     [0.000000000, 0.249999985, 0.249999985],
     [0.499999970, 0.499999970, 0.000000000],
     [0.374999970, 0.374999970, 0.374999970],
     [0.624999940, 0.624999940, 0.624999940],
     [0.374999970, 0.874999940, 0.874999940],
     [0.624999940, 0.124999993, 0.124999993],
     [0.874999940, 0.874999940, 0.374999970],
     [0.124999993, 0.124999993, 0.624999940],
     [0.874999940, 0.374999970, 0.874999940],
     [0.124999993, 0.624999940, 0.124999993],
     ]
        """
        )

        # create defect_structure
        params['config']['defect_structure'] = []

        nsublat = len(sublat)
        ntbl = len(_tbl['coords']) // nsublat

        print("number of sublattice = {}, number of sites for each = {}".format(nsublat, ntbl))

        for i in range(nsublat):
            t = {}
            t['coords'] = _tbl['coords'][i*ntbl:(i+1)*ntbl]
            t['groups'] = sublat[i]
            params['config']['defect_structure'].append(t)

        if chem:
            params['config']['chemical_potential'] = chem

        return params

class TestGrandCanonical(unittest.TestCase):
    def setUp(self):
        print()
        print("initialize random number sequence")
        random.seed(random_seed)
        pass

    def tearDown(self):
        pass

    def try_add_remove(self, proto, chem_tbl, mode, choice, src, expect_stat, expect_dst):

        # examples of arguments:
        #   proto -- list of species on each sublattice
        #     [ [ 'A' ], [ 'A', 'B' ] ]
        #   chem_tbl -- table of chemical potentials
        #     [ { 'species': ['A'], 'mu': 0.1 }, { 'species': ['A','B'], 'mu': 0.2 } ]
        #   mode -- add or remove
        #     add = True, remove = False
        #   choice -- move choice from chem table, specify id
        #     0 (for ['A'])
        #   src -- input state
        #     [ [ 4 ], [ 2, 2 ] ]
        #   stat -- expected status
        #     True/False
        #   expect -- expected result
        #     [ [[ 5 ], [ 2, 2 ]], [[ 4 ], [ 3, 2 ]] ]

        # create model
        model = DFTLatticeGas(StabCalculator(),
                              save_history = False,
                              enable_grandcanonical = True,
                              gc_ratio = 1.0,
                              debug = True)

        # proto = [ ['A'], ['A','B'] ]
        # -> [ { name: A, num: NA1 } ], [ { name: A, num: NA2 }, { name: B, num: NB } ]
        pp = [ [ { 'name': q, 'num': src[i][k] } for k,q in enumerate(p) ] for i,p in enumerate(proto) ]

        p = Params().create_params(pp, chem_tbl)
        # print(p)

        # create parameter set
        params = DFTConfigParams.from_dict(p['config'])

        # create config
        config = defect_config(params)

        # backup
        config_latgas_bak = [
            copy.deepcopy(sublat.latgas_rep) for sublat in config.defect_sublattices
        ]

        # try
        stat, info = model._try_add_remove(config, chem_tbl[choice], mode)
        print("--- stat =", stat, info)

        # new state
        config_latgas = [
            copy.deepcopy(sublat.latgas_rep) for sublat in config.defect_sublattices
        ]

        # check if move is successful/unsatisfied as expected
        self.assertEqual(stat, expect_stat)

        if stat:
            # check number of particles
            x = [ [ len([ i for i,v in enumerate(subl.latgas_rep) if v[0] == grp ])
                    for grp in proto[k] ]
                  for k,subl in enumerate(config.defect_sublattices) ]

            print("--- dst =", x)
            self.assertTrue(x in expect_dst)

            # check info to find if the changes are consistent with the info
            for lat,(rep0,rep1) in enumerate(zip(config_latgas_bak, config_latgas)):
                for rep_idx,(v0,v1) in enumerate(zip(rep0, rep1)):
                    for grp,sublat_id,idx in info['choice']:
                        if lat == sublat_id and rep_idx == idx:
                            if mode: # try add
                                self.assertTrue(v0[0] == '@' and v1[0] == grp)
                            else:
                                self.assertTrue(v0[0] == grp and v1[0] == '@')
                            break
                    else:
                        self.assertTrue(v0[0] == v1[0])
        else:
            # check info to find if the changes are consistent with the info
            for lat,(rep0,rep1) in enumerate(zip(config_latgas_bak, config_latgas)):
                for rep_idx,(v0,v1) in enumerate(zip(rep0, rep1)):
                    self.assertTrue(v0[0] == v1[0])

        pass

    def try_add(self, proto, chem_tbl, choice, src, expect_stat, expect_dst):
        return self.try_add_remove(proto, chem_tbl, True, choice, src, expect_stat, expect_dst)

    def try_remove(self, proto, chem_tbl, choice, src, expect_stat, expect_dst):
        return self.try_add_remove(proto, chem_tbl, False, choice, src, expect_stat, expect_dst)

    def try_exchange(self, proto, choice, src, nrepeat=10):

        # examples of arguments:
        #   proto -- list of species on each sublattice
        #     [ [ 'A' ], [ 'A', 'B' ] ]
        #   choice -- sublattice id
        #     0
        #   src -- input state
        #     [ [ 4 ], [ 2, 2 ] ]
        #   nrepeat -- number of runs
        #     10 (default)

        # create model
        model = DFTLatticeGas(StabCalculator(),
                              save_history = False,
                              enable_grandcanonical = False,
                              gc_ratio = 0.0,
                              debug = True)

        # proto = [ ['A'], ['A','B'] ]
        # -> [ { name: A, num: NA1 } ], [ { name: A, num: NA2 }, { name: B, num: NB } ]
        pp = [ [ { 'name': q, 'num': src[i][k] } for k,q in enumerate(p) ] for i,p in enumerate(proto) ]

        p = Params().create_params(pp)
        # print(p)

        # create parameter set
        params = DFTConfigParams.from_dict(p['config'])

        # create config
        config = defect_config(params)

        for iter in range(nrepeat):
            # backup
            config_latgas_bak = [
                copy.deepcopy(sublat.latgas_rep) for sublat in config.defect_sublattices
            ]

            # try
            stat, info = model._try_exchange(config.defect_sublattices[choice])
            print("--- stat =", stat, info)

            # new state
            config_latgas = [
                copy.deepcopy(sublat.latgas_rep) for sublat in config.defect_sublattices
            ]

            # check if move is successful/unsatisfied as expected
            #XXX self.assertEqual(stat, expect_stat)

            if stat:
                # check number of particles
                x = [ [ len([ i for i,v in enumerate(subl.latgas_rep) if v[0] == grp ])
                        for grp in proto[k] ]
                      for k,subl in enumerate(config.defect_sublattices) ]

                # print("--- dst =", x)
                self.assertEqual(x, src)

                # check info to find if the changes are consistent with the info
                ex1 = info['ex1_id']
                ex1_group = info['ex1_group']
                ex2 = info['ex2_id']
                ex2_group = info['ex2_group']

                for lat,(rep0,rep1) in enumerate(zip(config_latgas_bak, config_latgas)):
                    if lat == choice:
                        for rep_idx,(v0,v1) in enumerate(zip(rep0, rep1)):
                            if rep_idx == ex1:
                                self.assertEqual(v0[0], ex1_group)
                                self.assertEqual(v1[0], ex2_group)
                            elif rep_idx == ex2:
                                self.assertEqual(v0[0], ex2_group)
                                self.assertEqual(v1[0], ex1_group)
                            else:
                                self.assertEqual(v0[0], v1[0])
                    else:
                        for rep_idx,(v0,v1) in enumerate(zip(rep0, rep1)):
                            self.assertEqual(v0[0], v1[0])
            else:
                # check info to find if the changes are consistent with the info
                for lat,(rep0,rep1) in enumerate(zip(config_latgas_bak, config_latgas)):
                    for rep_idx,(v0,v1) in enumerate(zip(rep0, rep1)):
                        self.assertTrue(v0[0] == v1[0])

        pass

    def try_rotate(self, proto, group, choice, src, nrepeat=1):

        # examples of arguments:
        #   proto -- list of species on each sublattice
        #     [ [ 'A' ], [ 'A', 'B' ] ]
        #   group -- group structure
        #     { 'A': n } where 'A' has n species
        #   choice -- sublattice id
        #     0
        #   src -- input state
        #     [ [ 4 ], [ 2, 2 ] ]
        #   nrepeat -- number of runs
        #     10 (default)

        # create model
        model = DFTLatticeGas(StabCalculator(),
                              save_history = False,
                              enable_grandcanonical = False,
                              gc_ratio = 0.0,
                              debug = True)

        # proto = [ ['A'], ['A','B'] ]
        # -> [ { name: A, num: NA1 } ], [ { name: A, num: NA2 }, { name: B, num: NB } ]
        pp = [ [ { 'name': q, 'num': src[i][k] } for k,q in enumerate(p) ] for i,p in enumerate(proto) ]

        p = Params().create_params(pp)

        # append rotation degree of freedom
        for grp,n in group.items():
            if n > 1:
                species = [ grp ]
                coords = [ [[ 0.0, 0.0, 0.0 ]] for i in range(n) ]
                for sublat in p['config']['defect_structure']:
                    for g in sublat['groups']:
                        if g['name'] == grp:
                            g['species'] = species
                            g['coords'] = coords
        # print(p)

        # create parameter set
        params = DFTConfigParams.from_dict(p['config'])

        # create config
        config = defect_config(params)

        for iter in range(nrepeat):
            # backup
            config_latgas_bak = [
                copy.deepcopy(sublat.latgas_rep) for sublat in config.defect_sublattices
            ]

            # try
            stat, info = model._try_rotate(config.defect_sublattices[choice])
            print("--- stat =", stat, info)

            # new state
            config_latgas = [
                copy.deepcopy(sublat.latgas_rep) for sublat in config.defect_sublattices
            ]

            # check if move is successful/unsatisfied as expected
            #XXX self.assertEqual(stat, expect_stat)

            if stat:
                # check number of particles
                x = [ [ len([ i for i,v in enumerate(subl.latgas_rep) if v[0] == grp ])
                        for grp in proto[k] ]
                      for k,subl in enumerate(config.defect_sublattices) ]

                # print("--- dst =", x)
                self.assertEqual(x, src)

                # # check info to find if the changes are consistent with the info
                rot_id = info['rot_id']
                rot_group = info['rot_group']
                rot_from = info['from']
                rot_to = info['to']

                for lat,(rep0,rep1) in enumerate(zip(config_latgas_bak, config_latgas)):
                    if lat == choice:
                        for rep_idx,(v0,v1) in enumerate(zip(rep0, rep1)):
                            if rep_idx == rot_id:
                                self.assertEqual(v0[0], v1[0])
                                self.assertEqual(v0[0], rot_group)
                                self.assertEqual(v0[1], rot_from)
                                self.assertEqual(v1[1], rot_to)
                            else:
                                self.assertEqual(v0[0], v1[0])
                                self.assertEqual(v0[1], v1[1])
                    else:
                        for rep_idx,(v0,v1) in enumerate(zip(rep0, rep1)):
                            self.assertEqual(v0[0], v1[0])
                            self.assertEqual(v0[1], v1[1])
            else:
                # check info to find if the changes are consistent with the info
                for lat,(rep0,rep1) in enumerate(zip(config_latgas_bak, config_latgas)):
                    for rep_idx,(v0,v1) in enumerate(zip(rep0, rep1)):
                        self.assertTrue(v0[0] == v1[0])
                        self.assertTrue(v0[1] == v1[1])

        pass

    #----------------
    def test_rotate_1species(self):
        proto = [ ['A'] ]
        group = { 'A': 4 }
        N = 24

        if True:
        # # if False:
            self.try_rotate(proto, group, 0, [[N]])
            self.try_rotate(proto, group, 0, [[N-1]])
            self.try_rotate(proto, group, 0, [[0]])

    #----------------
    def test_rotate_2species(self):
        proto = [ ['A','B'] ]
        group = { 'A': 4 }
        N = 24

        if True:
        # # if False:
            self.try_rotate(proto, group, 0, [[N//2,N//2]])
            self.try_rotate(proto, group, 0, [[N,0]])
            self.try_rotate(proto, group, 0, [[0,N]])
            self.try_rotate(proto, group, 0, [[0,0]])

    #----------------
    def test_exchange_1species(self):
        proto = [ ['A'] ]
        N = 24

        if True:
        # if False:
            self.try_exchange(proto, 0, [[N]])
            self.try_exchange(proto, 0, [[N-1]])
            self.try_exchange(proto, 0, [[0]])

    #----------------
    def test_exchange_2species(self):
        proto = [ ['A', 'B'], ['A', 'B'] ]
        N = 12

        if True:
        # if False:
            self.try_exchange(proto, 1, [[N//2,N//2],[0,0]])
            self.try_exchange(proto, 1, [[N//2,N//2],[N,0]])
            self.try_exchange(proto, 1, [[N//2,N//2],[N//3,N//3]])
            self.try_exchange(proto, 1, [[N//2,N//2],[N//2,N//2]])

    #----------------
    def test_addremove_1sublat_1species(self):
        chem_tbl = [
            { 'species': [ 'A' ], 'mu': 0.1 },
        ]

        proto = [ ['A'] ]
        N = 24

        if True:
        # if False:
            # try add
            self.try_add(proto, chem_tbl, 0, [[N]],   False, [ [[-1]] ])
            self.try_add(proto, chem_tbl, 0, [[N-1]], True,  [ [[N]] ])
            self.try_add(proto, chem_tbl, 0, [[0]],   True,  [ [[1]] ])

            # try remove
            self.try_remove(proto, chem_tbl, 0, [[0]], False, [ [[-1]] ])
            self.try_remove(proto, chem_tbl, 0, [[1]], True,  [ [[0]] ])
            self.try_remove(proto, chem_tbl, 0, [[N]], True,  [ [[N-1]] ])

    #----------------
    def test_addremove_1sublat_2species(self):
        chem_tbl = [
            { 'species': [ 'A' ], 'mu': 0.1 },
            { 'species': [ 'B' ], 'mu': 0.2 },
            { 'species': [ 'A', 'B' ], 'mu': 0.3 },
        ]

        proto = [ ['A', 'B'] ]
        N = 24

        if True:
        # if False:
            # try add A
            self.try_add(proto, chem_tbl, 0, [[N, 0]],    False, [ [[-1, -1]] ])
            self.try_add(proto, chem_tbl, 0, [[N-1, 1]],  False, [ [[-1, -1]] ])
            self.try_add(proto, chem_tbl, 0, [[N-1, 0]],  True,  [ [[N, 0]] ])
            self.try_add(proto, chem_tbl, 0, [[0, N//2]], True,  [ [[1, N//2]] ])

            # try remove A
            self.try_remove(proto, chem_tbl, 0, [[0, 0]], False, [ [[-1, -1]] ])
            self.try_remove(proto, chem_tbl, 0, [[1, 0]], True,  [ [[0, 0]] ])
            self.try_remove(proto, chem_tbl, 0, [[N, 0]], True,  [ [[N-1, 0]] ])

            # try add A,B pair
            self.try_add(proto, chem_tbl, 2, [[N, 0]],    False, [ [[-1, -1]] ])
            self.try_add(proto, chem_tbl, 2, [[N-1, 0]],  False, [ [[-1, -1]] ])
            self.try_add(proto, chem_tbl, 2, [[N//2-1, N//2]],    False, [ [[-1, -1]] ])
            self.try_add(proto, chem_tbl, 2, [[N-2, 0]],  True,  [ [[N-1, 1]] ])
            self.try_add(proto, chem_tbl, 2, [[N//2-1, N//2-1]],  True,  [ [[N//2, N//2]] ])
            self.try_add(proto, chem_tbl, 2, [[0, 0]],    True,  [ [[1, 1]] ])

            # try remove A,B pair
            self.try_remove(proto, chem_tbl, 2, [[0, 0]], False, [ [[-1, -1]] ])
            self.try_remove(proto, chem_tbl, 2, [[1, 0]], False, [ [[-1, -1]] ])
            self.try_remove(proto, chem_tbl, 2, [[1, 1]], True,  [ [[0, 0]] ])
            self.try_remove(proto, chem_tbl, 2, [[N//2, N//2]], True,  [ [[N//2-1, N//2-1]] ])

    #----------------
    def test_addremove_2sublat_2species_disjoint(self):
        chem_tbl = [
            { 'species': [ 'A' ], 'mu': 0.1 },
            { 'species': [ 'B' ], 'mu': 0.2 },
            { 'species': [ 'A', 'B' ], 'mu': 0.3 },
        ]

        proto = [ ['A'], ['B'] ]
        N = 12

        if True:
        # if False:
            # try add A
            self.try_add(proto, chem_tbl, 0, [[N],[N]],    False, [])
            self.try_add(proto, chem_tbl, 0, [[N-1],[N]],  True,  [ [[N],[N]] ])
            self.try_add(proto, chem_tbl, 0, [[0],[N//2]], True,  [ [[1],[N//2]] ])

            # try remove A
            self.try_remove(proto, chem_tbl, 0, [[0],[N]], False, [])
            self.try_remove(proto, chem_tbl, 0, [[1],[N]], True,  [ [[0],[N]] ])
            self.try_remove(proto, chem_tbl, 0, [[N],[N]], True,  [ [[N-1],[N]] ])

            # try add A,B pair
            self.try_add(proto, chem_tbl, 2, [[N],[N]],    False, [])
            self.try_add(proto, chem_tbl, 2, [[N-1],[N]],  False, [])
            self.try_add(proto, chem_tbl, 2, [[N],[N-1]],  False, [])
            self.try_add(proto, chem_tbl, 2, [[N-1],[N-1]],  True, [ [[N],[N]] ])
            self.try_add(proto, chem_tbl, 2, [[N//2],[N//2]],  True, [ [[N//2+1],[N//2+1]] ])

            # try remove A,B pair
            self.try_remove(proto, chem_tbl, 2, [[0],[0]],    False, [])
            self.try_remove(proto, chem_tbl, 2, [[1],[0]],    False, [])
            self.try_remove(proto, chem_tbl, 2, [[0],[1]],    False, [])
            self.try_remove(proto, chem_tbl, 2, [[1],[1]],    True, [ [[0],[0]] ])
            self.try_remove(proto, chem_tbl, 2, [[N//2],[N//2]],   True, [ [[N//2-1],[N//2-1]] ])

    #----------------
    def test_addremove_2sublat_2species_multilayer(self):
        chem_tbl = [
            { 'species': [ 'A' ], 'mu': 0.1 },
            { 'species': [ 'B' ], 'mu': 0.2 },
            { 'species': [ 'A', 'B' ], 'mu': 0.3 },
        ]

        proto = [ ['A'], ['A', 'B'] ]
        N = 12

        if True:
        # if False:
            # try add A
            self.try_add(proto, chem_tbl, 0, [[N], [0,N]],   False, [])
            self.try_add(proto, chem_tbl, 0, [[N], [1,N-1]], False, [])
            self.try_add(proto, chem_tbl, 0, [[N], [N,0]],   False, [])
            self.try_add(proto, chem_tbl, 0, [[N], [0,N-1]], True,  [ [[N],[1,N-1]] ])
            self.try_add(proto, chem_tbl, 0, [[N], [1,N-2]], True,  [ [[N],[2,N-2]] ])
            self.try_add(proto, chem_tbl, 0, [[N-1], [0,N]], True,  [ [[N],[0,N]] ])
            self.try_add(proto, chem_tbl, 0, [[N-1], [0,N-1]], True,  [ [[N],[0,N-1]], [[N-1],[1,N-1]] ])

            # try remove A
            self.try_remove(proto, chem_tbl, 0, [[0], [0,N]],   False, [])
            self.try_remove(proto, chem_tbl, 0, [[1], [0,N]],   True,  [ [[0],[0,N]] ])
            self.try_remove(proto, chem_tbl, 0, [[0], [1,N-1]], True,  [ [[0],[0,N-1]] ])
            self.try_remove(proto, chem_tbl, 0, [[1], [1,N-1]], True,  [ [[0],[1,N-1]], [[1],[0,N-1]] ])

            # try add B
            self.try_add(proto, chem_tbl, 1, [[N], [0,N]],         False, [])
            self.try_add(proto, chem_tbl, 1, [[N], [1,N-1]],       False, [])
            self.try_add(proto, chem_tbl, 1, [[N], [N//2,N//2]],   False, [])
            self.try_add(proto, chem_tbl, 1, [[N], [0,N-1]],       True,  [ [[N],[0,N]] ])
            self.try_add(proto, chem_tbl, 1, [[N], [N//2,N//2-1]], True,  [ [[N],[N//2,N//2]] ])

            # try remove B
            self.try_remove(proto, chem_tbl, 1, [[N], [0,0]],   False, [])
            self.try_remove(proto, chem_tbl, 1, [[N], [N,0]],   False, [])
            self.try_remove(proto, chem_tbl, 1, [[N], [0,1]],   True, [ [[N],[0,0]] ])
            self.try_remove(proto, chem_tbl, 1, [[N], [0,N]],   True, [ [[N],[0,N-1]] ])

            # try add A,B pair
            self.try_add(proto, chem_tbl, 2, [[0], [0,N]],         False, [])
            self.try_add(proto, chem_tbl, 2, [[0], [1,N-1]],       False, [])
            self.try_add(proto, chem_tbl, 2, [[N], [0,N-1]],       False, [])
            self.try_add(proto, chem_tbl, 2, [[N], [1,N-2]],       False, [])
            self.try_add(proto, chem_tbl, 2, [[N], [1,N-2]],       False, [])

            self.try_add(proto, chem_tbl, 2, [[N], [0,N-2]],       True, [ [[N],[1,N-1]] ])
            self.try_add(proto, chem_tbl, 2, [[N-1], [0,N-1]],     True, [ [[N],[0,N]] ])
            self.try_add(proto, chem_tbl, 2, [[N//2], [0,N//2]],   True, [ [[N//2+1],[0,N//2+1]], [[N//2],[1,N//2+1]] ])

            # try remove A,B pair
            self.try_remove(proto, chem_tbl, 2, [[0], [0,0]], False, [])
            self.try_remove(proto, chem_tbl, 2, [[1], [0,1]], True, [ [[0],[0,0]] ])
            self.try_remove(proto, chem_tbl, 2, [[0], [1,1]], True, [ [[0],[0,0]] ])
            self.try_remove(proto, chem_tbl, 2, [[N//2], [N//2,N//2]], True, [ [[N//2-1],[N//2,N//2-1]], [[N//2],[N//2-1,N//2-1]] ])

    #----------------
    def test_input(self):
        rootdir = os.path.dirname(__file__)
        datadir = os.path.join(rootdir, "data", "gc")

        input_file = os.path.join(datadir, "input.toml")
        input_param = toml.load(input_file)

        nsite = 24
        nAl = 12
        nMg = 8
        nvac = nsite - nAl - nMg

        # create parameter set from input file
        params = DFTConfigParams.from_dict(input_param["config"])

        self.assertEqual(params.supercell, [1,1,1])
        self.assertEqual(len(params.defect_sublattices), 1)

        self.assertEqual(len(params.num_defects), 1)
        self.assertEqual(len(params.num_defects[0]), 2)  # [ Al, Mg ]
        self.assertEqual(params.num_defects[0]['Al'], nAl)
        self.assertEqual(params.num_defects[0]['Mg'], nMg)

        # create config
        config = defect_config(params)

        self.assertEqual(len(config.defect_sublattices), 1)
        self.assertEqual(len(config.defect_sublattices[0].groups), 3)  # [ Al, Mg, vac ]
        self.assertEqual(len(config.defect_sublattices[0].site_centers_sc), nsite)

        self.assertEqual(len(config.defect_sublattices[0].groups), 3)
        self.assertEqual(config.defect_sublattices[0].vac_group, '@')

        self.assertEqual(config.defect_sublattices[0].groups[0].name, 'Al')
        self.assertEqual(len(config.defect_sublattices[0].groups[0].species), 1)
        self.assertEqual(config.defect_sublattices[0].groups[0].species[0], 'Al')
        self.assertEqual(config.defect_sublattices[0].groups[0].orientations, 1)

        self.assertEqual(config.defect_sublattices[0].groups[1].name, 'Mg')
        self.assertEqual(len(config.defect_sublattices[0].groups[1].species), 1)
        self.assertEqual(config.defect_sublattices[0].groups[1].species[0], 'Mg')
        self.assertEqual(config.defect_sublattices[0].groups[1].orientations, 1)

        num_Al = len(match_latgas_group(config.defect_sublattices[0].latgas_rep, config.defect_sublattices[0].groups[0]))
        num_Mg = len(match_latgas_group(config.defect_sublattices[0].latgas_rep, config.defect_sublattices[0].groups[1]))
        num_vac = len(match_latgas_group(config.defect_sublattices[0].latgas_rep, config.defect_sublattices[0].groups[2]))

        self.assertEqual(num_Al, nAl)
        self.assertEqual(num_Mg, nMg)
        self.assertEqual(num_vac, nvac)

        self.assertEqual(len(config.chemical_potential), 3)
        self.assertEqual(config.chemical_potential[0]['species'], ['Al'])
        self.assertEqual(config.chemical_potential[0]['mu'], 0.1)
        self.assertEqual(config.chemical_potential[1]['species'], ['Mg'])
        self.assertEqual(config.chemical_potential[1]['mu'], 0.2)
        self.assertEqual(config.chemical_potential[2]['species'], ['Al', 'Mg'])
        self.assertEqual(config.chemical_potential[2]['mu'], 0.3)

        pass
