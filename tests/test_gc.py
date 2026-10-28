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
    def create_params(self, sublat, chem=None, gc_move=None):

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

        if gc_move:
            params['config']['grandcanonical_move'] = gc_move

        return params

class TestGrandCanonical(unittest.TestCase):
    def setUp(self):
        print()
        print("initialize random number sequence")
        random.seed(random_seed)
        pass

    def tearDown(self):
        pass

    def try_add_remove(self, model_desc, mode, choice, src, expect_stat, expect_dst):

        # examples of arguments:
        #   model_desc -- model description
        #     'proto' -- list of species on each sublattice
        #       [ [ 'A' ], [ 'A', 'B' ] ]
        #     'chem_tbl' -- table of chemical potentials
        #       [ { 'species': ['A'], 'mu': 0.1 }, { 'species': ['A','B'], 'mu': 0.2 } ]
        #     'move_tbl' -- table of grandcanonical moves
        #       [ { 'from': ['A'], 'to': [] }, { 'species': ['A','B'] } ]
        #   mode -- run mode
        #     "add" or "remove"
        #   choice -- move choice from chem table, specify id
        #     0 (for ['A'])
        #   src -- input state
        #     [ [ 4 ], [ 2, 2 ] ]
        #   exect_stat -- expected return status
        #     True/False
        #   expect_dst -- expected result
        #     [ [[ 5 ], [ 2, 2 ]], [[ 4 ], [ 3, 2 ]] ]

        # create model
        model = DFTLatticeGas(StabCalculator(),
                              save_history = False,
                              enable_grandcanonical = True,
                              gc_ratio = 1.0,
                              debug = True)

        proto = model_desc.get('proto', [])
        
        # proto = [ ['A'], ['A','B'] ]
        # -> [ { name: A, num: NA1 } ], [ { name: A, num: NA2 }, { name: B, num: NB } ]
        pp = [ [ { 'name': q, 'num': src[i][k] } for k,q in enumerate(p) ] for i,p in enumerate(proto) ]

        chem_tbl = model_desc.get('chem_tbl', None)
        move_tbl = model_desc.get('move_tbl', None)

        p = Params().create_params(pp, chem_tbl, move_tbl)
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
        move = config.gc_move[choice]  # from: [A] -> to: []
        if mode == "add":
            spec = (move['to'],move['from'])
            self.assertTrue(len(spec[0]) == 0)
        else:
            spec = (move['from'],move['to'])
            self.assertTrue(len(spec[1]) == 0)

        stat, info = model._try_add_remove(config, spec, mode=="add")
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
                    for grp,choices in info['choice']:
                        for ch in choices:
                            sublat_id, idx = ch
                            if lat == sublat_id and rep_idx == idx:
                                if mode == "add":
                                    self.assertTrue(v0[0] == '@' and v1[0] == grp)
                                else:
                                    self.assertTrue(v0[0] == grp and v1[0] == '@')
                                break
                        else:
                            continue
                        break
                    else:
                        self.assertTrue(v0[0] == v1[0])
        else:
            # check info to find if the changes are consistent with the info
            for lat,(rep0,rep1) in enumerate(zip(config_latgas_bak, config_latgas)):
                for rep_idx,(v0,v1) in enumerate(zip(rep0, rep1)):
                    self.assertTrue(v0[0] == v1[0])

        pass

    def try_add(self, model, choice, src, expect_stat, expect_dst):
        return self.try_add_remove(model, "add", choice, src, expect_stat, expect_dst)

    def try_remove(self, model, choice, src, expect_stat, expect_dst):
        #return self.try_add_remove(model, "remove", choice, src, expect_stat, expect_dst)
        return True

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

    def try_replace(self, choice, src, expect_stat, expect_dst, args={}):

        # examples of arguments:
        #   src -- input state
        #     [ [ 4 ], [ 2, 2 ] ]
        #   expect_stat -- expected return status
        #     True/False
        #   expect_dst  -- expected state
        #     [ [[ 5 ], [ 2, 2 ]], [[ 4 ], [ 3, 2 ]] ]

        # create model
        model = DFTLatticeGas(StabCalculator(),
                              save_history = False,
                              enable_grandcanonical = True,
                              gc_ratio = 1.0,
                              debug = True)

        proto = args.get('proto', [ ['A', 'B'] ])

        # proto = [ ['A'], ['A','B'] ] and src = [ [NA1], [NA2,NB] ]
        # -> [ { name: A, num: NA1 } ], [ { name: A, num: NA2 }, { name: B, num: NB } ]
        pp = [ [ { 'name': q, 'num': src[i][k] } for k,q in enumerate(p) ] for i,p in enumerate(proto) ]

        chem_tbl = args.get('chem_tbl', [
            { 'species': 'A', 'mu': 0.1 },
            { 'species': 'B', 'mu': 0.1 },
        ])

        move_tbl = args.get('move_tbl', [
            { 'from': 'A', 'to': 'B', 'reverse': True },
        ])
        
        p = Params().create_params(pp, chem_tbl, move_tbl)
        #print(p)

        # create parameter set
        params = DFTConfigParams.from_dict(p['config'])

        # create config
        config = defect_config(params)

        # backup
        config_latgas_bak = [
            copy.deepcopy(sublat.latgas_rep) for sublat in config.defect_sublattices
        ]

        # try
        move = config.gc_move[choice]
        stat, info = model._try_replace(config, (move['from'],move['to']))
        print("--- stat =", stat, info)

        # new state
        config_latgas = [
            copy.deepcopy(sublat.latgas_rep) for sublat in config.defect_sublattices
        ]

        # check if move is successful/unsatisfied as expected
        self.assertEqual(stat, expect_stat, "return status unexpected")

        if stat:
            # check number of particles
            x = [ [ len([ i for i,v in enumerate(subl.latgas_rep) if v[0] == grp ])
                    for grp in proto[k] ]
                  for k,subl in enumerate(config.defect_sublattices) ]

            print("--- dst =", x)
            self.assertTrue(x in expect_dst)

        else:
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
    def exchange_stats(self, proto, src, nrepeat):
        # run _try_exchange nrepeat times on sublattice 0 and
        # return the list of (stat, info)
        model = DFTLatticeGas(StabCalculator(),
                              save_history = False,
                              enable_grandcanonical = False,
                              gc_ratio = 0.0,
                              debug = True)
        pp = [ [ { 'name': q, 'num': src[i][k] } for k,q in enumerate(p) ] for i,p in enumerate(proto) ]
        params = DFTConfigParams.from_dict(Params().create_params(pp)['config'])
        config = defect_config(params)
        return [ model._try_exchange(config.defect_sublattices[0]) for _ in range(nrepeat) ]

    #----------------
    def test_exchange_skips_unoccupied_groups(self):
        # fully occupied sublattice: the vacancy group has zero sites and
        # must never be chosen, so every trial is a real exchange
        N = 24
        for proto, src in [ ([['A','B']], [[N//2, N//2]]),
                            ([['A','B','C']], [[N//3, N//3, N//3]]),
                            ([['A','B']], [[N, 0]]) ]:
            results = self.exchange_stats(proto, src, 200)
            if proto == [['A','B']] and src == [[N, 0]]:
                # only one occupied group: nothing to exchange
                self.assertTrue(all(stat is False for stat, info in results))
            else:
                self.assertTrue(all(stat is True for stat, info in results))
                chosen = set()
                for stat, info in results:
                    chosen.add(info['ex1_group'])
                    chosen.add(info['ex2_group'])
                self.assertEqual(chosen, set(proto[0]))

    #----------------
    def test_exchange_uses_vacancy_when_present(self):
        # with a vacancy the vacancy group is a regular exchange partner
        N = 24
        results = self.exchange_stats([['A','B']], [[N//2, N//2 - 1]], 300)
        self.assertTrue(all(stat is True for stat, info in results))
        chosen = set()
        for stat, info in results:
            chosen.add(info['ex1_group'])
            chosen.add(info['ex2_group'])
        self.assertEqual(chosen, {'A', 'B', '@'})

    #----------------
    def test_trialstep_rotate_only_when_single_group_occupied(self):
        # a sublattice filled by one orientable group must only try rotations
        N = 24
        model = DFTLatticeGas(StabCalculator(),
                              save_history = False,
                              enable_grandcanonical = False,
                              gc_ratio = 0.0,
                              debug = True)
        pp = [ [ { 'name': 'A', 'num': N } ] ]
        p = Params().create_params(pp)
        # give group A two orientations
        p['config']['defect_structure'][0]['groups'][0]['coords'] = [[[0,0,0]], [[0,0,0.1]]]
        params = DFTConfigParams.from_dict(p['config'])
        config = defect_config(params)
        model.energy(config)
        for _ in range(50):
            dconfig, dE, dW, info = model.trialstep(config, config.energy)
            self.assertEqual(info['mode'], 'rotate')
            self.assertEqual(info['result'], 'ok')
            self.assertEqual(dW, 0.0)

    #----------------
    def make_canonical_model_config(self, proto, src, orientations={}):
        # build a canonical model and a config; orientations maps a group
        # name to the number of orientations it should have
        model = DFTLatticeGas(StabCalculator(),
                              save_history = False,
                              enable_grandcanonical = False,
                              gc_ratio = 0.0,
                              debug = True)
        pp = [ [ { 'name': q, 'num': src[i][k] } for k,q in enumerate(p) ] for i,p in enumerate(proto) ]
        p = Params().create_params(pp)
        for ds in p['config']['defect_structure']:
            for g in ds['groups']:
                if g['name'] in orientations:
                    g['coords'] = [ [[0, 0, 0.1 * i]] for i in range(orientations[g['name']]) ]
        params = DFTConfigParams.from_dict(p['config'])
        config = defect_config(params)
        model.energy(config)
        return model, config

    #----------------
    def test_rotate_skips_unoccupied_orientable_group(self):
        # A has orientations but no sites; B and C are not orientable.
        # Rotation must not be proposed, every trial is an exchange.
        N = 24
        model, config = self.make_canonical_model_config(
            [['A', 'B', 'C']], [[0, N//2, N//2]], orientations={'A': 3})
        for _ in range(100):
            dconfig, dE, dW, info = model.trialstep(config, config.energy)
            self.assertEqual(info['mode'], 'exchange')
            self.assertEqual(info['result'], 'ok')
            self.assertEqual(dW, 0.0)

    #----------------
    def test_trialstep_skips_immobile_sublattice(self):
        # sublattice 0 is filled by one non-orientable group and cannot move;
        # sublattice 1 can exchange. Only sublattice 1 must be chosen.
        # (the template has 24 sites in total, split evenly over sublattices)
        N = 12
        model, config = self.make_canonical_model_config(
            [['A'], ['A', 'B']], [[N], [N//2, N//2]])
        for _ in range(100):
            dconfig, dE, dW, info = model.trialstep(config, config.energy)
            self.assertEqual(info['sublat_id'], 1)
            self.assertEqual(info['result'], 'ok')
            self.assertEqual(dW, 0.0)

    #----------------
    def test_trialstep_raises_when_nothing_can_move(self):
        # canonical sampling with no movable degree of freedom is an input error
        N = 24
        model, config = self.make_canonical_model_config([['A']], [[N]])
        from abics.exception import InputError
        with self.assertRaises(InputError):
            model.trialstep(config, config.energy)

    #----------------
    def test_canonical_mc_counts_failed_trials(self):
        # a model whose trialstep cannot propose a move (dW = NaN) must be
        # counted in nfail and never accepted by CanonicalMonteCarlo
        from abics.sampling.mc import CanonicalMonteCarlo

        class NullModel:
            def energy(self, config):
                return 0.0
            def trialstep(self, config, energy_now):
                return None, 0.0, float("nan")
            def newconfig(self, config, dconfig):
                raise AssertionError("failed trial must not be accepted")

        mc = CanonicalMonteCarlo(NullModel(), 1.0, config=None)
        mc.energy = 0.0
        mc.MCstep(10)
        self.assertEqual(mc.ntrials, 10)
        self.assertEqual(mc.naccepted, 0)
        self.assertEqual(mc.nfail, 10)

    #----------------
    def test_canonical_mc_accepts_two_value_trialstep(self):
        # models that return only (dconfig, dE) keep working
        from abics.sampling.mc import CanonicalMonteCarlo

        class TwoValueModel:
            def energy(self, config):
                return 0.0
            def trialstep(self, config, energy_now):
                return "new", -1.0
            def newconfig(self, config, dconfig):
                return dconfig

        mc = CanonicalMonteCarlo(TwoValueModel(), 1.0, config="old")
        mc.energy = 0.0
        mc.MCstep(5)
        self.assertEqual(mc.ntrials, 5)
        self.assertEqual(mc.naccepted, 5)
        self.assertEqual(mc.nfail, 0)
        self.assertEqual(mc.config, "new")

    #----------------
    def test_gc_failed_trial_returns_nan_weight(self):
        # grandcanonical: removing a species that is absent cannot be proposed
        model = DFTLatticeGas(StabCalculator(),
                              save_history = False,
                              enable_grandcanonical = True,
                              gc_ratio = 1.0,
                              debug = True)
        N = 24
        pp = [ [ { 'name': 'A', 'num': 0 } ] ]
        chem_tbl = [ { 'species': ['A'], 'mu': 0.1 } ]
        move_tbl = [ { 'from': ['A'], 'to': [] } ]
        p = Params().create_params(pp, chem_tbl, move_tbl)
        params = DFTConfigParams.from_dict(p['config'])
        config = defect_config(params)
        model.energy(config)
        results = set()
        for _ in range(50):
            dconfig, dE, dW, info = model.trialstep(config, config.energy)
            results.add((info['mode'], info['result'], np.isnan(dW)))
        # the reversible move picks add or remove; remove must fail with NaN
        self.assertIn(('remove', 'nochoice', True), results)
        self.assertIn(('add', 'ok', False), results)
        self.assertTrue(all(is_nan == (res == 'nochoice') for mode, res, is_nan in results))

    #----------------
    def test_write_acceptance_ratio(self):
        import tempfile
        from abics.sampling.mc import write_acceptance_ratio
        with tempfile.TemporaryDirectory() as d:
            fname = os.path.join(d, "acceptance_ratio.dat")
            write_acceptance_ratio(fname, [100.0, 200.0], [10, 0], [4, 0], [2, 0])
            rows = [ line.split() for line in open(fname) if not line.startswith("#") ]
        self.assertEqual(rows[0], ["100.0", "0.4", "10", "4", "2"])
        self.assertEqual(rows[1][0], "200.0")
        self.assertEqual(rows[1][1], "nan")
        self.assertEqual(rows[1][2:], ["0", "0", "0"])

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
        model = {
            'proto': [ ['A'] ],
            'chem_tbl': [
                { 'species': [ 'A' ], 'mu': 0.1 },
            ],
            'move_tbl': [
                { 'species': 'A' },
            ],
        }
        N = 24

        if True:
        # if False:
            # try add
            self.try_add(model, 0, [[N]],   False, [ [[-1]] ])
            self.try_add(model, 0, [[N-1]], True,  [ [[N]] ])
            self.try_add(model, 0, [[0]],   True,  [ [[1]] ])

            # try remove
            self.try_remove(model, 0, [[0]], False, [ [[-1]] ])
            self.try_remove(model, 0, [[1]], True,  [ [[0]] ])
            self.try_remove(model, 0, [[N]], True,  [ [[N-1]] ])

    #----------------
    def test_addremove_1sublat_2species(self):
        model = {
            'proto': [ ['A', 'B'] ],
            'chem_tbl': [
                { 'species': [ 'A' ], 'mu': 0.1 },
                { 'species': [ 'B' ], 'mu': 0.2 },
                { 'species': [ 'A', 'B' ], 'mu': 0.3 },
            ],
        }
        N = 24

        if True:
        # if False:
            # try add A
            self.try_add(model, 0, [[N, 0]],    False, [ [[-1, -1]] ])
            self.try_add(model, 0, [[N-1, 1]],  False, [ [[-1, -1]] ])
            self.try_add(model, 0, [[N-1, 0]],  True,  [ [[N, 0]] ])
            self.try_add(model, 0, [[0, N//2]], True,  [ [[1, N//2]] ])

            # try remove A
            self.try_remove(model, 0, [[0, 0]], False, [ [[-1, -1]] ])
            self.try_remove(model, 0, [[1, 0]], True,  [ [[0, 0]] ])
            self.try_remove(model, 0, [[N, 0]], True,  [ [[N-1, 0]] ])

            # try add A,B pair
            self.try_add(model, 2, [[N, 0]],    False, [ [[-1, -1]] ])
            self.try_add(model, 2, [[N-1, 0]],  False, [ [[-1, -1]] ])
            self.try_add(model, 2, [[N//2-1, N//2]],    False, [ [[-1, -1]] ])
            self.try_add(model, 2, [[N-2, 0]],  True,  [ [[N-1, 1]] ])
            self.try_add(model, 2, [[N//2-1, N//2-1]],  True,  [ [[N//2, N//2]] ])
            self.try_add(model, 2, [[0, 0]],    True,  [ [[1, 1]] ])

            # try remove A,B pair
            self.try_remove(model, 2, [[0, 0]], False, [ [[-1, -1]] ])
            self.try_remove(model, 2, [[1, 0]], False, [ [[-1, -1]] ])
            self.try_remove(model, 2, [[1, 1]], True,  [ [[0, 0]] ])
            self.try_remove(model, 2, [[N//2, N//2]], True,  [ [[N//2-1, N//2-1]] ])

    #----------------
    def test_addremove_2sublat_2species_disjoint(self):
        model = {
            'proto': [ ['A'], ['B'] ],
            'chem_tbl': [
                { 'species': [ 'A' ], 'mu': 0.1 },
                { 'species': [ 'B' ], 'mu': 0.2 },
                { 'species': [ 'A', 'B' ], 'mu': 0.3 },
            ],
        }
        N = 12

        if True:
        # if False:
            # try add A
            self.try_add(model, 0, [[N],[N]],    False, [])
            self.try_add(model, 0, [[N-1],[N]],  True,  [ [[N],[N]] ])
            self.try_add(model, 0, [[0],[N//2]], True,  [ [[1],[N//2]] ])

            # try remove A
            self.try_remove(model, 0, [[0],[N]], False, [])
            self.try_remove(model, 0, [[1],[N]], True,  [ [[0],[N]] ])
            self.try_remove(model, 0, [[N],[N]], True,  [ [[N-1],[N]] ])

            # try add A,B pair
            self.try_add(model, 2, [[N],[N]],    False, [])
            self.try_add(model, 2, [[N-1],[N]],  False, [])
            self.try_add(model, 2, [[N],[N-1]],  False, [])
            self.try_add(model, 2, [[N-1],[N-1]],  True, [ [[N],[N]] ])
            self.try_add(model, 2, [[N//2],[N//2]],  True, [ [[N//2+1],[N//2+1]] ])

            # try remove A,B pair
            self.try_remove(model, 2, [[0],[0]],    False, [])
            self.try_remove(model, 2, [[1],[0]],    False, [])
            self.try_remove(model, 2, [[0],[1]],    False, [])
            self.try_remove(model, 2, [[1],[1]],    True, [ [[0],[0]] ])
            self.try_remove(model, 2, [[N//2],[N//2]],   True, [ [[N//2-1],[N//2-1]] ])

    #----------------
    def test_addremove_2sublat_2species_multilayer(self):
        model = {
            'proto': [ ['A'], ['A', 'B'] ],
            'chem_tbl': [
                { 'species': [ 'A' ], 'mu': 0.1 },
                { 'species': [ 'B' ], 'mu': 0.2 },
                { 'species': [ 'A', 'B' ], 'mu': 0.3 },
            ],
        }
        N = 12

        if True:
        # if False:
            # try add A
            self.try_add(model, 0, [[N], [0,N]],   False, [])
            self.try_add(model, 0, [[N], [1,N-1]], False, [])
            self.try_add(model, 0, [[N], [N,0]],   False, [])
            self.try_add(model, 0, [[N], [0,N-1]], True,  [ [[N],[1,N-1]] ])
            self.try_add(model, 0, [[N], [1,N-2]], True,  [ [[N],[2,N-2]] ])
            self.try_add(model, 0, [[N-1], [0,N]], True,  [ [[N],[0,N]] ])
            self.try_add(model, 0, [[N-1], [0,N-1]], True,  [ [[N],[0,N-1]], [[N-1],[1,N-1]] ])

            # try remove A
            self.try_remove(model, 0, [[0], [0,N]],   False, [])
            self.try_remove(model, 0, [[1], [0,N]],   True,  [ [[0],[0,N]] ])
            self.try_remove(model, 0, [[0], [1,N-1]], True,  [ [[0],[0,N-1]] ])
            self.try_remove(model, 0, [[1], [1,N-1]], True,  [ [[0],[1,N-1]], [[1],[0,N-1]] ])

            # try add B
            self.try_add(model, 1, [[N], [0,N]],         False, [])
            self.try_add(model, 1, [[N], [1,N-1]],       False, [])
            self.try_add(model, 1, [[N], [N//2,N//2]],   False, [])
            self.try_add(model, 1, [[N], [0,N-1]],       True,  [ [[N],[0,N]] ])
            self.try_add(model, 1, [[N], [N//2,N//2-1]], True,  [ [[N],[N//2,N//2]] ])

            # try remove B
            self.try_remove(model, 1, [[N], [0,0]],   False, [])
            self.try_remove(model, 1, [[N], [N,0]],   False, [])
            self.try_remove(model, 1, [[N], [0,1]],   True, [ [[N],[0,0]] ])
            self.try_remove(model, 1, [[N], [0,N]],   True, [ [[N],[0,N-1]] ])

            # try add A,B pair
            self.try_add(model, 2, [[0], [0,N]],         False, [])
            self.try_add(model, 2, [[0], [1,N-1]],       False, [])
            self.try_add(model, 2, [[N], [0,N-1]],       False, [])
            self.try_add(model, 2, [[N], [1,N-2]],       False, [])
            self.try_add(model, 2, [[N], [1,N-2]],       False, [])

            self.try_add(model, 2, [[N], [0,N-2]],       True, [ [[N],[1,N-1]] ])
            self.try_add(model, 2, [[N-1], [0,N-1]],     True, [ [[N],[0,N]] ])
            self.try_add(model, 2, [[N//2], [0,N//2]],   True, [ [[N//2+1],[0,N//2+1]], [[N//2],[1,N//2+1]] ])

            # try remove A,B pair
            self.try_remove(model, 2, [[0], [0,0]], False, [])
            self.try_remove(model, 2, [[1], [0,1]], True, [ [[0],[0,0]] ])
            self.try_remove(model, 2, [[0], [1,1]], True, [ [[0],[0,0]] ])
            self.try_remove(model, 2, [[N//2], [N//2,N//2]], True, [ [[N//2-1],[N//2,N//2-1]], [[N//2],[N//2-1,N//2-1]] ])

    #----------------
    def test_addremove_1sublat_1species_2elem(self):
        model = {
            'proto': [ ['A'] ],
            'chem_tbl': [
                { 'species': [ 'A' ], 'mu': 0.1 },
            ],
            'move_tbl': [
                { 'from': ['A','A'], 'to': [] },
            ],
        }
        N = 24

        if True:
        # if False:
            # try add
            self.try_add(model, 0, [[N]],   False, [ [[-1]] ])
            self.try_add(model, 0, [[N-1]], False, [ [[-1]] ])
            self.try_add(model, 0, [[N-2]], True,  [ [[N]] ])
            self.try_add(model, 0, [[0]],   True,  [ [[2]] ])

            # try remove
            self.try_remove(model, 0, [[0]], False, [ [[-1]] ])
            self.try_remove(model, 0, [[1]], False, [ [[-1]] ])
            self.try_remove(model, 0, [[2]], True,  [ [[0]] ])
            self.try_remove(model, 0, [[N]], True,  [ [[N-2]] ])

    #----------------
    def test_addremove_1sublat_2species_2elem(self):
        model = {
            'proto': [ ['A','B'] ],
            'chem_tbl': [
                { 'species': [ 'A' ], 'mu': 0.1 },
                { 'species': [ 'B' ], 'mu': 0.1 },
            ],
            'move_tbl': [
                { 'from': ['A','A','B'], 'to': [] },
            ],
        }
        N = 24

        if True:
        # if False:
            # try add
            self.try_add(model, 0, [[N,0]],   False, [ [[-1]] ])
            self.try_add(model, 0, [[N-1,0]], False, [ [[-1]] ])
            self.try_add(model, 0, [[N-3,0]], True,  [ [[N-1,1]] ])
            self.try_add(model, 0, [[0,0]],   True,  [ [[2,1]] ])

            # try remove
            self.try_remove(model, 0, [[0,0]], False, [ [[-1]] ])
            self.try_remove(model, 0, [[1,1]], False, [ [[-1]] ])
            self.try_remove(model, 0, [[2,1]], True,  [ [[0,0]] ])
            self.try_remove(model, 0, [[N,0]], False, [ [[-1]] ])
            self.try_remove(model, 0, [[N-1,1]], True,  [ [[N-3,0]] ])

    #----------------
    def test_replace(self):
        N = 24

        model = {
            'proto': [ ['A','B','C'] ],
            'chem_tbl': [
                {'species': 'A', 'mu': 0.1},
                {'species': 'B', 'mu': 0.1},
                {'species': 'C', 'mu': 0.1},
            ],
            'move_tbl': [
                {'from': 'A', 'to': 'C'},
            ],
        }

        self.try_replace(0, [[0,0,N]], False, [], model)
        self.try_replace(0, [[1,0,N-1]], True, [ [[0,0,N]] ], model)
        self.try_replace(0, [[N,0,0]], True, [ [[N-1,0,1]] ], model)

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

        self.assertEqual(len(config.gc_move), 2)
        self.assertEqual(config.gc_move[0]['from'], [('Mg',1)])
        self.assertEqual(config.gc_move[0]['to'], [])
        self.assertEqual(config.gc_move[0]['reverse'], True)
        self.assertEqual(config.gc_move[1]['from'], [('Al',1)])
        self.assertEqual(config.gc_move[1]['to'], [('Mg',1)])
        self.assertEqual(config.gc_move[1]['reverse'], True)

        pass
