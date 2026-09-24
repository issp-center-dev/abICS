import copy
import unittest

import numpy as np
from pymatgen.core import Structure

from abics.applications.latgas_abinitio_interface.model_setup import (
    Config,
    DFTLatticeGas,
    DefectSublattice,
    Group,
    _latgas_index,
    _patch_indices,
)


class DummyCalculator:
    def submit(self, structure, output_dir):
        return 0.0, structure


def make_config(cellsize=(2, 2, 2)):
    base = Structure(
        np.diag(cellsize),
        [],
        [],
    )
    first = DefectSublattice(
        [[0.0, 0.0, 0.0]],
        [
            Group("A", ["A"], coords=[[[0.0, 0.0, 0.0]]]),
            Group("B", ["B"], coords=[[[0.0, 0.0, 0.0]]]),
        ],
    )
    second = DefectSublattice(
        [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]],
        [
            Group("C", ["C"], coords=[[[0.0, 0.0, 0.0]]]),
            Group("D", ["D"], coords=[[[0.0, 0.0, 0.0]]]),
        ],
    )
    return Config(
        base,
        [first, second],
        [{"A": 4, "B": 4}, {"C": 8, "D": 8}],
        cellsize=list(cellsize),
    )


class TestPatchExchange(unittest.TestCase):
    def make_model(self, shapes=((1, 1, 1),), debug=True, mode="single_sublattice"):
        return DFTLatticeGas(
            DummyCalculator(),
            save_history=False,
            patch_exchange_enable=True,
            patch_exchange_ratio=1.0,
            patch_exchange_shapes=shapes,
            patch_exchange_mode=mode,
            debug=debug,
        )

    def test_flat_index_and_wrapped_order(self):
        cellsize = (2, 2, 2)
        ids = _patch_indices((1, 1, 1), (2, 2, 2), cellsize, 2)
        expected = [
            _latgas_index(ix, iy, iz, ib, cellsize, 2)
            for dx in range(2)
            for dy in range(2)
            for dz in range(2)
            for ib in range(2)
            for ix, iy, iz in [((1 + dx) % 2, (1 + dy) % 2, (1 + dz) % 2)]
        ]
        self.assertEqual(ids, expected)

    def test_displacements_are_symmetric(self):
        model = self.make_model()
        config = make_config()
        model._prepare_patch_geometry(config.cellsize)
        valid = set(model._patch_displacements[(2, 2, 2)][(1, 1, 1)])
        for displacement in valid:
            inverse = tuple((-value) % 2 for value in displacement)
            self.assertIn(inverse, valid)

    def test_patch_preserves_composition_and_second_sublattice(self):
        model = self.make_model()
        config = make_config()
        before = [
            copy.deepcopy(sublat.latgas_rep)
            for sublat in config.defect_sublattices
        ]
        for _ in range(20):
            model._try_patch_exchange(config)
        for before_sublat, sublat in zip(before, config.defect_sublattices):
            for group in set(entry[0] for entry in before_sublat):
                self.assertEqual(
                    sum(entry[0] == group for entry in before_sublat),
                    sum(entry[0] == group for entry in sublat.latgas_rep),
                )

    def test_fixed_exchange_is_an_involution(self):
        config = make_config()
        sublat = config.defect_sublattices[0]
        shape = (1, 1, 1)
        origin1 = (0, 0, 0)
        origin2 = (1, 0, 0)
        ids1 = _patch_indices(origin1, shape, config.cellsize, 1)
        ids2 = _patch_indices(origin2, shape, config.cellsize, 1)
        before = copy.deepcopy(sublat.latgas_rep)
        for index1, index2 in zip(ids1, ids2):
            sublat.latgas_rep[index1], sublat.latgas_rep[index2] = (
                sublat.latgas_rep[index2],
                sublat.latgas_rep[index1],
            )
        for index1, index2 in zip(ids1, ids2):
            sublat.latgas_rep[index1], sublat.latgas_rep[index2] = (
                sublat.latgas_rep[index2],
                sublat.latgas_rep[index1],
            )
        self.assertEqual(before, sublat.latgas_rep)

    def test_null_patch_does_not_change_state(self):
        config = make_config()
        for entry in config.defect_sublattices[0].latgas_rep:
            entry[:] = ["A", 0]
        for entry in config.defect_sublattices[1].latgas_rep:
            entry[:] = ["C", 0]
        model = self.make_model(shapes=((1, 1, 1),))
        before = [
            copy.deepcopy(sublat.latgas_rep)
            for sublat in config.defect_sublattices
        ]
        status, info = model._try_patch_exchange(config)
        self.assertFalse(status)
        self.assertEqual(info["result"], "null")
        self.assertEqual(
            before,
            [sublat.latgas_rep for sublat in config.defect_sublattices],
        )


    def test_all_sublattices_preserves_composition_per_sublattice(self):
        model = self.make_model(mode="all_sublattices")
        config = make_config()
        before = [
            copy.deepcopy(sublat.latgas_rep)
            for sublat in config.defect_sublattices
        ]
        for _ in range(20):
            model._try_patch_exchange(config)
        for before_sublat, sublat in zip(before, config.defect_sublattices):
            for group in set(entry[0] for entry in before_sublat):
                self.assertEqual(
                    sum(entry[0] == group for entry in before_sublat),
                    sum(entry[0] == group for entry in sublat.latgas_rep),
                )

    def test_all_sublattices_is_an_involution(self):
        config = make_config()
        shape = (1, 1, 1)
        origin1 = (0, 0, 0)
        origin2 = (1, 0, 0)
        before = [
            copy.deepcopy(sublat.latgas_rep)
            for sublat in config.defect_sublattices
        ]
        for _ in range(2):
            for sublat in config.defect_sublattices:
                nbasis = len(sublat.site_centers)
                ids1 = _patch_indices(origin1, shape, config.cellsize, nbasis)
                ids2 = _patch_indices(origin2, shape, config.cellsize, nbasis)
                rep = sublat.latgas_rep
                for index1, index2 in zip(ids1, ids2):
                    rep[index1], rep[index2] = rep[index2], rep[index1]
        self.assertEqual(
            before,
            [sublat.latgas_rep for sublat in config.defect_sublattices],
        )

    def test_all_sublattices_null_requires_every_sublattice_identical(self):
        config = make_config()
        # Make only the first sublattice's candidate patches identical;
        # the second sublattice still differs, so the combined move is not null.
        for entry in config.defect_sublattices[0].latgas_rep:
            entry[:] = ["A", 0]
        model = self.make_model(shapes=((1, 1, 1),), mode="all_sublattices")
        found_ok = False
        for _ in range(50):
            before = [
                copy.deepcopy(sublat.latgas_rep)
                for sublat in config.defect_sublattices
            ]
            status, info = model._try_patch_exchange(config)
            if status:
                found_ok = True
                self.assertEqual(info["result"], "ok")
                break
            self.assertEqual(
                before,
                [sublat.latgas_rep for sublat in config.defect_sublattices],
            )
        self.assertTrue(found_ok, "expected at least one non-null all_sublattices trial")

    def test_all_sublattices_null_when_every_sublattice_identical(self):
        config = make_config()
        for entry in config.defect_sublattices[0].latgas_rep:
            entry[:] = ["A", 0]
        for entry in config.defect_sublattices[1].latgas_rep:
            entry[:] = ["C", 0]
        model = self.make_model(shapes=((1, 1, 1),), mode="all_sublattices")
        before = [
            copy.deepcopy(sublat.latgas_rep)
            for sublat in config.defect_sublattices
        ]
        status, info = model._try_patch_exchange(config)
        self.assertFalse(status)
        self.assertEqual(info["result"], "null")
        self.assertEqual(
            before,
            [sublat.latgas_rep for sublat in config.defect_sublattices],
        )


if __name__ == "__main__":
    unittest.main()
