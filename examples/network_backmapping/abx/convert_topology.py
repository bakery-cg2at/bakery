#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#


import argparse
from md_libs import files_io

parser = argparse.ArgumentParser('Convert CG topology for AB2 system')
parser.add_argument('in_top', help='Input topology')
parser.add_argument('out_top', help='Output topology')

res_type2name = {
    ('MA', 'ML', 'MB'): ('CMS', ('MA', 'ML', 'MB')),
    ('PA', 'PL', 'RB'): ('FOU', ('PA1', 'PL1', 'RB1')),
    ('DA', 'PL', 'RB'): ('LCU', ('DA2', 'PL2', 'RB2')),
    ('RA', 'PL', 'RB'): ('DEU', ('RA3', 'PL3', 'RB3')),
    ('DA', 'PL', 'PB'): ('TEU', ('DA4', 'PL4', 'PB4')),
    ('RA', 'PL', 'PB'): ('LEU', ('RA5', 'PL5', 'PB5'))
}


window_size = [len(x) for x in res_type2name]
if len(window_size) != window_size.count(window_size[0]):
    raise RuntimeError('Wrong res_type2name, window size has to be the same')
print('Window size: {}'.format(window_size))
window_size = window_size[0]

args = parser.parse_args()

in_top = files_io.GROMACSTopologyFile(args.in_top)
in_top.read()

at_ids = sorted(in_top.atoms.keys())

for i in range(0, len(at_ids), window_size):
    chunk_at_ids = at_ids[i:i+window_size]
    chunk_types = tuple(in_top.atoms[x].atom_type for x in chunk_at_ids)
    new_res_name, new_atom_names = res_type2name[chunk_types]
    print((chunk_types, new_res_name))
    for ai, at_id in enumerate(chunk_at_ids):
        in_top.atoms[at_id].chain_name = new_res_name
        in_top.atoms[at_id].name = new_atom_names[ai]

in_top.write(args.out_top)
