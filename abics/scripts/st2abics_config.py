import toml
import os, sys
import argparse
from pymatgen import Structure


def main():

    # Input parser
    parser = argparse.ArgumentParser(
        description='Prepare abICS config from structure file')

    parser.add_argument('inputfi',
                        help='toml input file for %(prog)s',
    )
    parser.add_argument('structurefi',
                        help='Structure file that can be read' + \
                        ' by pymatgen Structure.from_file() method',
    )
    parser.add_argument('outfi', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Output file to be used as abics input.' + \
                        ' Defaults to standard output',
    )
                        
    args = parser.parse_args()
    abics_input_dict = {}

    # Read input files
    st2abics_params = toml.load(args.inputfi)
    input_structure = Structure.from_file(args.structurefi)
    input_species = list(input_structure.symbol_set)

    # Make supercell
    input_structure.make_supercell(st2abics_params['supercell'])

    # Make base structure
    base_params = st2abics_params['base_structure']
    base_species = base_params['species']
    fix = base_params['fix']

    base_structure = input_structure.copy()
    remove_species = [sp for sp in input_species if sp not in base_species]
    base_structure.remove_species(remove_species)
    abics_input_dict['config'] = {
        'unitcell': base_structure.lattice.matrix.tolist(),
        'supercell': [1,1,1],
    }
    for sp in base_structure.symbol_set:
        if fix:
            relaxations = [[False, False, False]]*int(base_structure.composition[sp])
        else:
            relaxations = [[True, True, True]]*int(base_structure.composition[sp])
        abics_input_dict['config']['base_structure'] = [
            {
            'type': sp,
            'coords': base_structure.frac_coords[
                base_structure.indices_from_symbol(sp),:
                ].tolist(),
            'relaxation': relaxations
            }
        ]

    # Make defect structures
    abics_input_dict['config']['defect_structure'] = []
    defect_params = st2abics_params['defect_structure']
    for defect_sublattice in defect_params:
        site_center_species = defect_sublattice['site_center_species']
        defect_structure = input_structure.copy()
        remove_species = [
            sp for sp in input_species if sp not in site_center_species
            ]
        defect_structure.remove_species(remove_species)
        coords = defect_structure.frac_coords.tolist()
        groups = []
        if 'groups' in defect_sublattice.keys():
            for group in defect_sublattice['groups']:
                groups.append(
                    {'name': group['name'],
                     'species': group.get('species', group['name']),
                     'coords': group.get('coords', [[[0,0,0]]]),
                     'num': int(group['num']),
                     }
                )
        else:
            for sp in site_center_species:
                groups.append(
                    {'name': sp,
                     'species': [sp],
                     'coords': [[[0,0,0]]],
                     'num': int(defect_structure.composition[sp]),
                     }
                )
        num = 0
        for group in groups:
            num += group['num']
        assert(num == len(coords))

        abics_input_dict['config']['defect_structure'].append(
            {'coords': coords,
             'groups': groups,
             }
        )

    config = toml.dumps(abics_input_dict)
    with open(os.path.join(os.path.dirname(__file__), 'input_template.toml'), 'r') as fi:
        template = fi.read()
    args.outfi.write(template + config)

if __name__ == '__main__':
    main()
    
