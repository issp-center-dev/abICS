import toml
import os, sys
import argparse
from pymatgen import Structure


def main():

    # Input parser
    parser = argparse.ArgumentParser(
        description="Prepare abICS config from structure file"
    )

    parser.add_argument("inputfi", help="toml input file for %(prog)s")
    parser.add_argument(
        "structurefi",
        help="Structure file that can be read"
        + " by pymatgen Structure.from_file() method",
    )
    parser.add_argument(
        "outfi",
        nargs="?",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file to be used as abics input." + " Defaults to standard output",
    )

    args = parser.parse_args()
    abics_input_dict = {}

    # Read input files
    st2abics_params = toml.load(args.inputfi)
    input_structure = Structure.from_file(args.structurefi)
    input_species = list(input_structure.symbol_set)

    # Make supercell
    input_structure.make_supercell(st2abics_params["supercell"])

    # Make base structure
    base_params = st2abics_params["base_structure"]
    base_species = base_params["species"]
    fix = base_params["fix"]

    base_structure = input_structure.copy()
    remove_species = [sp for sp in input_species if sp not in base_species]
    base_structure.remove_species(remove_species)
    abics_input_dict["config"] = {
        "unitcell": base_structure.lattice.matrix.tolist(),
        "supercell": [1, 1, 1],
    }
    if len(base_structure.symbol_set) == 0:
        abics_input_dict["config"]["base_structure"] = [{}]
    else:
        abics_input_dict["config"]["base_structure"] = []
    for sp in base_structure.symbol_set:
        if fix:
            relaxations = [[False, False, False]] * int(base_structure.composition[sp])
        else:
            relaxations = [[True, True, True]] * int(base_structure.composition[sp])
        abics_input_dict["config"]["base_structure"].append(
            {
                "type": sp,
                "coords": base_structure.frac_coords[
                    base_structure.indices_from_symbol(sp), :
                ].tolist(),
                "relaxation": relaxations,
            }
        )

    # Make defect structures
    abics_input_dict["config"]["defect_structure"] = []
    defect_params = st2abics_params["defect_structure"]
    for defect_sublattice in defect_params:
        site_center_species = defect_sublattice["site_center_species"]
        defect_structure = input_structure.copy()
        remove_species = [sp for sp in input_species if sp not in site_center_species]
        defect_structure.remove_species(remove_species)
        coords = defect_structure.frac_coords.tolist()
        groups = []
        if "groups" in defect_sublattice.keys():
            for group in defect_sublattice["groups"]:
                groups.append(
                    {
                        "name": group["name"],
                        "species": group.get("species", [group["name"]]),
                        "coords": group.get("coords", [[[0, 0, 0]]]),
                        "num": int(group["num"]),
                    }
                )
        else:
            for sp in site_center_species:
                groups.append(
                    {
                        "name": sp,
                        "species": [sp],
                        "coords": [[[0, 0, 0]]],
                        "num": int(defect_structure.composition[sp]),
                    }
                )
        num = 0
        for group in groups:
            num += group["num"]
        assert num == len(coords)

        abics_input_dict["config"]["defect_structure"].append(
            {"coords": coords, "groups": groups}
        )

    config = ["[config]"]
    config.append("unitcell = [")
    for i in range(3):
        r = abics_input_dict["config"]["unitcell"][i]
        config.append("[{}, {}, {}],".format(r[0], r[1], r[2]))
    config.append("]")
    config.append("supercell = [1, 1, 1]")
    config.append("")
    if abics_input_dict["config"]["base_structure"] == [{}]:
        config.append("[[config.base_structure]]")
    else:
        for base in abics_input_dict["config"]["base_structure"]:
            config.append("[[config.base_structure]]")
            config.append('type = "{}"'.format(base["type"]))
            config.append('coords = """')
            for r in base["coords"]:
                config.append("{} {} {}".format(r[0], r[1], r[2]))
            config.append('"""')
            config.append('relaxation = """')
            for r in base["relaxation"]:
                r = list(map(lambda x: str(x).lower(), r))
                config.append("{} {} {}".format(r[0], r[1], r[2]))
            config.append('"""')
            config.append("")
    for defect in abics_input_dict["config"]["defect_structure"]:
        config.append("[[config.defect_structure]]")
        config.append('coords = """')
        for r in defect["coords"]:
            config.append("{} {} {}".format(r[0], r[1], r[2]))
        config.append('"""')
        config.append("")
        for group in defect["groups"]:
            config.append("[[config.defect_structure.groups]]")
            config.append('name = "{}"'.format(group["name"]))
            config.append("species = {}".format(group["species"]))
            config.append("coords = [")
            for coord in group["coords"]:
                if len(coord) == 1:
                    c = coord[0]
                    config.append('"{} {} {}"'.format(c[0], c[1], c[2]))
                else:
                    config.append('"""')
                    for c in coord:
                        config.append("{} {} {}".format(c[0], c[1], c[2]))
                    config.append('""",')
            config.append("]")
            config.append("num = {}".format(group["num"]))
            config.append("")

    # config = toml.dumps(abics_input_dict)
    with open(
        os.path.join(os.path.dirname(__file__), "input_template.toml"), "r"
    ) as fi:
        template = fi.read()
    # args.outfi.write(template + config)
    args.outfi.write(template)
    for c in config:
        args.outfi.write(c)
        args.outfi.write("\n")


if __name__ == "__main__":
    main()
