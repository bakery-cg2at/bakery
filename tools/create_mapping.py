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
from xml.etree import ElementTree
from xml.dom import minidom


def prettify(elem):
    """Return a pretty-printed XML string for the Element."""
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


def _args():
    parser = argparse.ArgumentParser('Creates xml mapping file for long polymer chains.')
    parser.add_argument('--config', required=True)
    parser.add_argument('--out', required=True)

    return parser.parse_args()


def _prepare_cg_beads(cg_beads, config):
    for bead in config['CG_BEADS']:
        bead_item = ElementTree.SubElement(cg_beads, 'cg_bead')
        bead_name = ElementTree.SubElement(bead_item, 'name')
        bead_name.text = bead['name']
        bead_type = ElementTree.SubElement(bead_item, 'type')
        bead_type.text = bead['type']
        bead_mapping = ElementTree.SubElement(bead_item, 'mapping')
        bead_mapping.text = bead['mapping']
        bead_beads = ElementTree.SubElement(bead_item, 'beads')
        bead_beads.text = bead['beads']


def _prepare_cg_bonded(cg_bonded, config):
    config_xml = {'BONDS': 'bond', 'ANGLES': 'angle', 'DIHEDRALS': 'dihedral'}
    for cname, xml_name in config_xml.iteritems():
        for bond_name in config.get(cname, {}):
            bond_item = ElementTree.SubElement(cg_bonded, xml_name)
            bond_name_item = ElementTree.SubElement(bond_item, 'name')
            bond_name_item.text = bond_name
            bond_beads = ElementTree.SubElement(bond_item, 'beads')
            bond_beads.text = '\n'.join(config[cname][bond_name])


def _prepare_cg_mapping(maps, config):
    for map_name in config.get('MAPPING', {}):
        map_item = ElementTree.SubElement(maps, 'map')
        map_name_item = ElementTree.SubElement(map_item, 'name')
        map_name_item.text = map_name
        map_weights = ElementTree.SubElement(map_item, 'weights')
        map_weights.text = config['MAPPING'][map_name]


def main():
    args = _args()

    config = {}
    execfile(args.config, {}, config)

    root = ElementTree.Element('cg_molecule')
    mol_name = ElementTree.SubElement(root, 'name')
    mol_name.text = config['NAME']
    mol_ident = ElementTree.SubElement(root, 'ident')
    mol_ident.text = config['IDENT']

    topology = ElementTree.SubElement(root, 'topology')

    cg_beads = ElementTree.SubElement(topology, 'cg_beads')
    _prepare_cg_beads(cg_beads, config)

    cg_bonded = ElementTree.SubElement(topology, 'cg_bonded')
    _prepare_cg_bonded(cg_bonded, config)

    maps = ElementTree.SubElement(root, 'maps')
    _prepare_cg_mapping(maps, config)

    with open(args.out, 'w') as out:
        out.write(prettify(root))

if __name__ == '__main__':
    main()
