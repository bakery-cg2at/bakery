#!/usr/bin/env python
"""
Copyright (C) 2015 Jakub Krajniak <jkrajniak@gmail.com>

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
import math


def _args():
    parser = argparse.ArgumentParser(
        description='Converts LAMMPS tabulated potentials to Espresso++ format')
    parser.add_argument('input_file')
    parser.add_argument('output_file')
    parser.add_argument('--table_type', choices=('pair', 'bond', 'angle', 'dihedral'),
                        default='pair')

    return parser.parse_args()


def _pair_convert(input_f, output_f):
    for line in input_f:
        l = line.strip()
        if not l.startswith('#') and not l.startswith('N') and l:
            sl = l.split()
            if len(sl) == 4:
                output_f.write('{} {} {}\n'.format(
                    sl[1], float(sl[2])*4.184, float(sl[3])*41.84))
            elif len(sl) == 3:
                output_f.write('{} {} {}\n'.format(
                    sl[0], float(sl[1])*4.184, float(sl[2])*41.84))


def _bond_convert(input_f, output_f):
    for line in input_f:
        l = line.strip()
        if not l.startswith('#') and not l.startswith('N') and l:
            sl = l.split()
            if len(sl) < 3:
                continue
            output_f.write('{} {} {}\n'.format(
                sl[0], float(sl[1])*4.184, float(sl[2])*41.84))


def _angle_convert(input_f, output_f):
    for line in input_f:
        l = line.strip()
        if not l.startswith('#') and not l.startswith('N') and l:
            sl = l.split()
            if len(sl) < 3:
                continue
            output_f.write('{} {} {}\n'.format(
                math.radians(float(sl[0])), float(sl[1])*4.184,
                float(sl[2])*4.184*180.0/math.pi))


def _dihedral_convert(input_f, output_f):
    degrees = True
    nof = False
    data = []
    for line in input_f:
        l = line.strip()
        if not l.startswith('#') and not l.startswith('N') and l:
            sl = l.split()
            if len(sl) < 2:
                continue
            phi = float(sl[1])
            if degrees:
                phi = math.radians(phi) - math.pi
            if nof:
                data.append([phi, float(sl[2])*4.184, 0.0])
            else:
                output_f.write('{} {} {}\n'.format(
                    phi, float(sl[2])*4.184, float(sl[3])*4.184*180.0/math.pi))
        elif l.startswith('N'):
            degrees = 'RADIANS' not in l
            nof = 'NOF' in l
    # Calculate force and then write a file.
    if data:
        for idx in range(0, len(data)-1):
            data[idx][2] = (data[idx+1][1] - data[idx][1])/(data[idx+1][0] - data[idx][0])
            output_f.write('{} {} {}\n'.format(*data[idx]))


def main():
    args = _args()
    input_f = open(args.input_file, 'r')
    output_f = open(args.output_file, 'w')

    table_type2func = {
        'pair': _pair_convert,
        'bond': _bond_convert,
        'angle': _angle_convert,
        'dihedral': _dihedral_convert
    }

    table_type2func[args.table_type](input_f, output_f)

    output_f.close()
    input_f.close()

if __name__ == '__main__':
    main()
