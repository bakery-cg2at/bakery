"""
Copyright (C) 2015-2016 Jakub Krajniak <jkrajniak@gmail.com>

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__doc__ = """Set of function that are used to process LAMMPS input files."""


def generate_cg_bonded_terms(lmp_input, ff_definition, output_topology, plain=False):
    """Generates CG bonded terms (LAMMPS)

    Warning: This function has side effect, it will modify `output_topology` object.

    Args:
        lmp_input: LAMMPSReader object.
        ff_definition: The definition of force-field.
        output_topology: GROMACSTopology object that will be updates.
        plain: If set to True then instead of cross_ terms it will generate regular one.
    """
    b_prefix = '' if plain else 'cross_'

    for bond_name, v in lmp_input.topology.items():
        for b_type, b_lists in v.items():
            for l in b_lists:
                l_new = list(map(output_topology.cg_old_new_id.get, l))
                assert None not in l_new
                properties = ff_definition[bond_name][b_type]
                output_topology.new_data[
                    '{}{}'.format(b_prefix, bond_name)][tuple(l_new)] = properties


def lammps2gromacs_ff(lmp_input):
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

    output_ff = {'bonds': {}, 'angles': {}, 'dihedrals': {}}
    if 'bond' in lmp_input.force_field:
        bond_style = lmp_input.force_field['bond_style'][0]
        if bond_style == 'harmonic':
            coeff = {k: list(map(float, v)) for k, v in lmp_input.force_field['bond'].items()}
            for k in coeff:
                K = kcal2kJ(coeff[k][0]) * 100.0
                r = coeff[k][1] / 10.0
                output_ff['bonds'][k] = [1, r, K, '; cg term']
        elif bond_style == 'table':
            for k, v in lmp_input.force_field['bond'].items():
                output_ff['bonds'][k] = [8, k, 0.0, '; cg term', v]
        else:
            raise RuntimeError('bond_style {} not supported'.format(bond_style))

    if 'angle' in lmp_input.force_field:
        angle_style = lmp_input.force_field['angle_style'][0]
        if angle_style == 'harmonic':
            coeff = {k: list(map(float, v)) for k, v in lmp_input.force_field['angle'].items()}
            for k in coeff:
                K = kcal2kJ(coeff[k][0])
                output_ff['angles'][k] = [1, coeff[k][1], K, '; cg term']
        elif angle_style == 'table':
            for k, v in lmp_input.force_field['angle'].items():
                output_ff['angles'][k] = [8, k, 0.0, '; cg term', v]
        else:
            raise RuntimeError('angle_style {} not supported'.format(angle_style))

    if 'dihedral' in lmp_input.force_field:
        dihedral_style = lmp_input.force_field['dihedral_style'][0]
        if dihedral_style == 'harmonic':
            coeff = {k: list(map(float, v)) for k, v in lmp_input.force_field['dihedral'].items()}
            for k in coeff:
                K = kcal2kJ(coeff[k][0])
                output_ff['dihedrals'][k] = [1, coeff[k][1], K, '; cg term']
        elif dihedral_style == 'table':
            for k, v in lmp_input.force_field['dihedral'].items():
                output_ff['dihedrals'][k] = [8, k, 0.0, '; cg term', v]
        else:
            raise RuntimeError('dihedral_style {} not supported'.format(dihedral_style))
    return output_ff
