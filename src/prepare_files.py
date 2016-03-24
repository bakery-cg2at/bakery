#!/usr/bin/env python
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

import argparse
import backmapping
import files_io
import gromacs
import lammps
import structures
import tools

__doc__ = 'Prepare step of bakery'


def _args():
    parser = argparse.ArgumentParser(
        description='Prepares hybrid coordinate and topology files.',
        add_help=True)

    parser.add_argument('--options', help='XML options file', required=True)
    parser.add_argument('--plain', help='Prepare non-hybrid configuration files.',
                        default=False, choices=('CG', 'AA'))
    parser.add_argument('--prefix', default='', help='Prefix for output files')

    return parser


def main():
    args = _args().parse_args()

    bck_settings = structures.BackmapperSettings(args.options)
    cg_com, cg_aa = backmapping.calculate_com_fragments(bck_settings)
    is_lammps = False
    if bck_settings.cg_configuration['format'] == 'LAMMPS':
        is_lammps = True
        print('Processing LAMMPS input files')
        lr = files_io.LammpsReader()
        for in_f in bck_settings.cg_configuration['input_files']:
            print('Reading {}'.format(in_f))
            lr.read_input(in_f)
        for idx, in_d in enumerate(bck_settings.cg_configuration['data_files']):
            print('Reading {}'.format(in_d))
            lr.read_data(in_d, update=idx > 0)
        cg_graph = lr.get_graph(bck_settings)
        force_field = lammps.lammps2gromacs_ff(lr)
    elif bck_settings.cg_configuration['format'] == 'GROMACS':
        print('Processing GROMACS input files')
        cg_graph = tools.get_graph(bck_settings)
    else:
        raise RuntimeError('Unknown CG format {}'.format(
            bck_settings.cg_configuration['format']))

    # Processing flow:
    # 1. Prepares hybrid coordinates and hybrid topology file.
    hyb_file, output_topology = backmapping.prepare_hyb_coordinates(
        bck_settings, cg_graph, cg_com, cg_aa, plain=args.plain)
    # 2. Renumerate bonds because of new ids of atoms.
    backmapping.update_atomistic_bonds(output_topology, bck_settings, plain=args.plain)
    # 3. Generate bonded terms for CG particles.
    if not args.plain or args.plain == 'CG':
        if is_lammps:
            lammps.generate_cg_bonded_terms(lr, force_field, output_topology, args.plain == 'CG')
        else:
            gromacs.generate_cg_bonded_terms(
                bck_settings, cg_graph, output_topology,
                args.plain == 'CG')
    # 4. Generate cross-link terms for atomistic particles
    if not args.plain or args.plain == 'AA':
        backmapping.generate_crosslink_at_terms(bck_settings, output_topology, plain=args.plain)
    # 5. Postprocess topology.
    backmapping.postprocess_topology(bck_settings, output_topology)
    # 6. Write output.
    print('Writing {}{}'.format(args.prefix, bck_settings.hybrid_topology['file']))
    output_topology.write('{}{}'.format(args.prefix, bck_settings.hybrid_topology['file']))
    print('Writing {}{}'.format(args.prefix, bck_settings.hybrid_configuration['file']))
    hyb_file.write(
        '{}{}'.format(args.prefix, bck_settings.hybrid_configuration['file']), force=True)


if __name__ == '__main__':
    main()
