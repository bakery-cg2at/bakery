"""
Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>

This file is part of Backmapper.

Backmapper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import collections
import cPickle
import networkx
import xml.etree.ElementTree as etree

import files_io
import lammps

__doc__ = 'Convert LAMMPS to GROMACS'


def read_settings(input_file):
    tree = etree.parse(input_file)
    root = tree.getroot()
    type2chain = {}
    output_type = collections.namedtuple('output_type', ['chain_name', 'type_name'])
    for e in root.find('type2chain').text.split():
        type_id, chain_name, type_name = e.split(':')
        type2chain[int(type_id)] = output_type(chain_name.strip(), type_name.strip())
    name_sequence = {}
    for e in root.findall('name_seq'):
        name_sequence[e.attrib['chain_name']] = e.text.split()

    molecule = root.find('molecule_type')
    output_tuple = collections.namedtuple(
        'Output', ['name_seq', 'type2chain', 'mol_name', 'nrexcl', 'nrmols'])
    return output_tuple(
        name_sequence,
        type2chain,
        molecule.attrib['name'],
        molecule.attrib['nrexcl'],
        molecule.attrib['nrmols'])


def lammps2gromacs_ff(lmp_input, settings):
    """Converts bonded parameters from LAMMPS input file to GROMACS input format. Currently only
    limited subset of LAMMPS *_style is supported.

    Args:
        lmp_input: The LAMMPSReader object.

    Returns:
        The dictionary with topology.
    """

    def kcal2kJ(v):
        """Converts kcal -> kJ"""
        return v*4.184

    output_ff = {'bonds': {}, 'angles': {}, 'dihedrals': {}, 'atomtypes': {}}
    if 'bond_style' in lmp_input.force_field:
        bond_style = lmp_input.force_field['bond_style'][0]
        if bond_style == 'harmonic':
            coeff = {k: map(float, v) for k, v in lmp_input.force_field['bond'].iteritems()}
            for k in coeff:
                K = kcal2kJ(coeff[k][0]) * 100.0
                r = coeff[k][1] / 10.0
                output_ff['bonds'][k] = [1, r, K, '; cg term']
        elif bond_style == 'table':
            for k, v in lmp_input.force_field['bond'].iteritems():
                output_ff['bonds'][k] = [8, k, 0.0, '; cg term', v]
        else:
            raise RuntimeError('bond_style {} not supported'.format(bond_style))

    if 'angle_style' in lmp_input.force_field:
        angle_style = lmp_input.force_field['angle_style'][0]
        if angle_style == 'harmonic':
            coeff = {k: map(float, v) for k, v in lmp_input.force_field['angle'].iteritems()}
            for k in coeff:
                K = kcal2kJ(coeff[k][0])
                output_ff['angles'][k] = [1, coeff[k][1], K, '; cg term']
        elif angle_style == 'table':
            for k, v in lmp_input.force_field['angle'].iteritems():
                output_ff['angles'][k] = [8, k, 0.0, '; cg term', v]
        else:
            raise RuntimeError('angle_style {} not supported'.format(angle_style))

    if 'dihedral_style' in lmp_input.force_field:
        dihedral_style = lmp_input.force_field['dihedral_style'][0]
        if dihedral_style == 'harmonic':
            coeff = {k: map(float, v) for k, v in lmp_input.force_field['dihedral'].iteritems()}
            for k in coeff:
                K = kcal2kJ(coeff[k][0])
                output_ff['dihedrals'][k] = [1, coeff[k][1], K, '; cg term']
        elif dihedral_style == 'table':
            for k, v in lmp_input.force_field['dihedral'].iteritems():
                output_ff['dihedrals'][k] = [8, k, 0.0, '; cg term', v]
        else:
            raise RuntimeError('dihedral_style {} not supported'.format(dihedral_style))

    if 'pair_style' in lmp_input.force_field:
        pair_style = lmp_input.force_field['pair_style'][0]
        if pair_style == 'table':

            for pair_type, output_type in settings.type2chain.items():
                output_ff['atomtypes'][output_type.type_name] = {
                    'name': output_type.type_name,
                    'mass': lmp_input._mass_type[pair_type],
                    'charge': lmp_input.atom_charges[pair_type],
                    'type': 'A',
                    'sigma': 1.0,
                    'epsilon': 1.0
                }
        else:
            raise RuntimeError('pair_style {} not supported'.format(pair_style))

    return output_ff


def build_gromacs(lr, settings, output_file):
    """Prepares GROMASC topology file. Not everything is supported!"""

    output = files_io.GROMACSTopologyFile(output_file)
    output.init()

    # Write defaults
    output.defaults = {
        'nbfunc': 1,
        'comb-rule': 1,
        'gen-pairs': 'no',
        'fudgeLJ': 0.0,
        'fudgeQQ': 0.0
    }
    print('Warning, [ defaults ] set to {}'.format(output.defaults))

    # Write atoms.
    # Stores sequence index for chains.
    seq_idx = {k: 0 for k in settings.name_seq}
    for atid in sorted(lr.atoms):
        at_data = lr.atoms[atid]
        chain_name, at_type = settings.type2chain[at_data['atom_type']]
        at_seq = settings.name_seq[chain_name]
        chain_len = len(at_seq)
        at_name = at_seq[seq_idx[chain_name] % chain_len]
        output.atoms[atid] = files_io.TopoAtom(
            atom_id=atid,
            atom_type=at_type,
            chain_idx=at_data['res_id'],
            chain_name=chain_name,
            name=at_name,
            cgnr=atid,
            charge=at_data['charge'] if at_data['charge'] else 0.0,
            mass=at_data['mass']
        )
        seq_idx[chain_name] += 1
    # Create bonded lists.
    lammps_ff = lammps2gromacs_ff(lr, settings)

    # Write atomtypes
    output.atomtypes = lammps_ff['atomtypes']

    for bond_type, bond_list in lr.topology['bonds'].items():
        for bp in sorted(bond_list):
            output.bonds[bp] = lammps_ff['bonds'][bond_type][:3]
    for angle_type, angle_list in lr.topology['angles'].items():
        for bp in sorted(angle_list):
            output.angles[bp] = lammps_ff['angles'][angle_type][:3]
    for dihedral_type, dihedral_list in lr.topology['dihedrals'].items():
        for bp in sorted(dihedral_list):
            output.dihedrals[bp] = lammps_ff['dihedrals'][dihedral_type][:3]
    for dihedral_type, dihedral_list in lr.topology['impropers'].items():
        for bp in sorted(dihedral_list):
            output.improper_dihedrals[bp] = lammps_ff['dihedrals'][dihedral_type][:3]

    output.moleculetype = {
        'name': settings.mol_name,
        'nrexcl': settings.nrexcl
    }
    output.molecules = {
        'name': settings.mol_name,
        'mol': settings.nrmols
    }
    output.system_name = settings.mol_name
    return output


def prepare_coordinate(lr, file_name, input_topology):
    out_coordinate = files_io.GROFile(file_name)
    out_coordinate.box = (lr.box['x'], lr.box['y'], lr.box['z'])
    for at_id, at_data in input_topology.atoms.items():
        out_coordinate.atoms[at_id] = files_io.Atom(
            atom_id=at_id,
            name=at_data.name,
            chain_name=at_data.chain_name,
            chain_idx=at_data.chain_idx,
            position=lr.atoms[at_id]['position']
        )
    out_coordinate.write(force=True)


def _args():
    parser = argparse.ArgumentParser(
        description='Convert LAMMPS to GROMACS topology',
        add_help=True)

    parser.add_argument('--data', help='LAMMPS data file', required=True, dest='lammps_data')
    parser.add_argument('--in', help='LAMMPS in file', required=True, dest='lammps_in')
    parser.add_argument('--out', help='GROMACS out file', required=True)
    parser.add_argument('--out_coordinate', help='.gro file', required=True)
    parser.add_argument('--out_graph', help='Write NetworkX graph', required=False)
    parser.add_argument('--out_lammps', help='Write cPickled LammpsReader object')
    parser.add_argument('--options', help='Options file', required=True)

    return parser


def main():
    args = _args().parse_args()
    lammps_reader = files_io.LammpsReader()
    lammps_reader.read_input(args.lammps_in)
    lammps_reader.read_data(args.lammps_data, update=True)
    settings = read_settings(args.options)
    if args.out_graph:
        out_graph = lammps_reader.get_graph(settings)
        networkx.write_gpickle(out_graph, args.out_graph)
    if args.out_lammps:
        cPickle.dump({'atoms': lammps_reader.atoms, 'topology': lammps_reader.topology},
                     open(args.out_lammps, 'wb'))
    output_topology = build_gromacs(lammps_reader, settings, args.out)
    output_topology.write()

    # Write coordinate file.
    prepare_coordinate(lammps_reader, args.out_coordinate, output_topology)


if __name__ == '__main__':
    main()
