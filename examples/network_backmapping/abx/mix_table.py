#! /usr/bin/env python
#  Copyright (C) 2016
#      Jakub Krajniak (jkrajniak at gmail.com)
#
#  This file is part of ChemLab.
#
#  ChemLab is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ChemLab is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import datetime
import numpy as np

from chemlab import gromacs_topology

parser = argparse.ArgumentParser('Mix table')
parser.add_argument('--conversion', help='Scalling factor', type=float, default=0.5)
parser.add_argument('--table_mono')
parser.add_argument('--table_poly')
parser.add_argument('--out')

args = parser.parse_args()


def mix_arithmetic(tab_mono, tab_poly, coupling):
    max_length = 0
    if tab_poly.shape[0] != tab_mono.shape[0]:
        if tab_poly.shape[0] <= tab_mono.shape[0]:
            out_tab = np.array(tab_poly)
        else:
            out_tab = np.array(tab_mono)
        if (tab_poly[:, 0][:out_tab.shape[0]] != tab_mono[:, 0][:out_tab.shape[0]]).all():
            raise RuntimeError('Both r columns should be the same')
    else:
        out_tab = np.array(tab_poly)
    max_length = out_tab.shape[0]
    if max_length == 0:
        raise RuntimeError('The length of output table is zero???')
    out_tab[:, 1] = coupling*tab_poly[:max_length, 1] + (1.0-coupling)*tab_mono[:max_length, 1]
    out_tab[:, 2] = coupling*tab_poly[:max_length, 2] + (1.0-coupling)*tab_mono[:max_length, 2]
    return out_tab


mono_tab = np.loadtxt(args.table_mono)
poly_tab = np.loadtxt(args.table_poly)
mixed_table = mix_arithmetic(mono_tab, poly_tab, args.conversion)
print('Saved {}'.format(out_name))
np.savetxt(
    args.out,
    mixed_table,
    header='Mixed of {} and {} at {}'.format(
        t1, t2, datetime.datetime.now()),
    fmt='%2.9e')
