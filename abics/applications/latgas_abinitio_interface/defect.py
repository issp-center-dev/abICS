from pymatgen import Lattice

from .model_setup import config, defect_sublattice, base_structure

from abics.util import read_coords


class DFTConfigParams:
    def __init__(self):
        pass

    @classmethod
    def from_dict(cls, d):
        if "config" in d:
            d = d["config"]
        params = cls()
        params.lat = Lattice(read_coords(d["unitcell"]))
        params.supercell = d.get("supercell", [1, 1, 1])
        params.base_structure = base_structure(params.lat, d["base_structure"])
        params.defect_sublattices = [
            defect_sublattice.from_dict(ds) for ds in d["defect_structure"]
        ]
        params.num_defects = [
            {g["name"]: g["num"] for g in ds["groups"]} for ds in d["defect_structure"]
        ]
        return params

    @classmethod
    def from_toml(cls, f):
        import toml

        return cls.from_dict(toml.load(f))


def defect_config(cparam: DFTConfigParams):
    spinel_config = config(
        cparam.base_structure,
        cparam.defect_sublattices,
        cparam.num_defects,
        cparam.supercell,
    )
    spinel_config.shuffle()
    return spinel_config
