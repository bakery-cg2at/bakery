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

import collections
import files_io
import xml.etree.ElementTree as etree

__doc__ = "Other data structures."""

BeadID = collections.namedtuple('BeadID', ['name', 'degree'])

XMLMap = collections.namedtuple(
    'XMLMap',
    ['ident',  # Identity name of molecules.
     'name',   # Name of molecules.
     'mass_map',  # Map with mass of CG beads.
     'source_coordinates',  # File name of atomistic fragments.
     'source_topology',  # File name of topology for atomistic fragments.
     'molecule_beads',  # Returns molecule CG beads.
     'molecule_topology'  # Topology at the level of CG bead.
     ])


class BackmapperSettings(object):
    def __init__(self, input_xml):
        tree = etree.parse(input_xml)
        self.root = tree.getroot()
        self._parse()

    def _parse(self):
        cg_molecules = [self._parse_cg_molecule(r) for r in self.root.findall('cg_molecule')]
        self.cg_molecules = {r.name: r for r in cg_molecules}
        self.hybrid_topology = self._parse_hybrid_topology()

        cg_configuration = self.root.find('cg_configuration')
        self.cg_configuration = {
            'format': cg_configuration.find('format').text.strip()
            }
        if self.cg_configuration['format'] == 'LAMMPS':
            self.cg_configuration.update(
                {'data_files': cg_configuration.find('data_files').text.strip().split(),
                 'input_files': cg_configuration.find('input_files').text.strip().split()})
            self.cg_configuration['LAMMPS'] = {
                'name_seq': {},
                'type2chain': {}
                }
            name_seqs = cg_configuration.findall('name_seq')
            for ns in name_seqs:
                cn = ns.attrib['chain_name']
                self.cg_configuration['LAMMPS']['name_seq'][cn] = ns.text.strip().split()
            for k in cg_configuration.find('type2chain').text.strip().split():
                kt, kn = k.strip().split(':')
                self.cg_configuration['LAMMPS']['type2chain'][int(kt)] = kn
        elif self.cg_configuration['format'] == 'GROMACS':
            self.cg_configuration['file'] = cg_configuration.find('file').text.strip()
            try:
                self.cg_configuration['topology'] = cg_configuration.find('topology').text.strip()
            except AttributeError:
                self.cg_configuration['topology'] = None
        hybrid_configuration = self.root.find('hybrid_configuration')
        self.hybrid_configuration = {
            'file': hybrid_configuration.find('file').text.strip(),
            'format': hybrid_configuration.find('format').text.strip()
        }

    def _parse_cg_molecule(self, root):
        source_file = root.find('source_file').text.strip()
        source_topology = root.find('source_topology').text.strip()
        chain_name = root.find('name').text.strip()
        ident_name = root.find('ident').text.strip()

        # For each of cg bead and different degree holds a list of weights in the same
        # order as the <beads> section
        mapping2mass = {
            x.find('name').text: map(float, x.find('weights').text.strip().split())
            for x in root.iter('map')
        }

        # Prepares cg_bead_mass and cg_beads structures.
        cg_bead_mass = {}
        cg_beads = {}
        for x in root.iter('cg_bead'):
            name = x.find('name').text.strip()
            atom_type = x.find('type').text.strip()
            try:
                degree = int(x.find('degree').text.strip())
            except:
                degree = '*'
            chain_name = chain_name
            mapping = x.find('mapping').text.strip()
            beads = x.find('beads').text.strip().split()
            active_site = x.find('active_site')
            if active_site is not None:
                active_site = active_site.text.strip()
            cg_bead_mass[BeadID(name, degree)] = {
                b_name: mapping2mass[mapping][i] for i, b_name in enumerate(beads)
            }
            cg_beads[BeadID(name, degree)] = files_io.TopoAtom(
                atom_type=atom_type,
                name=name,
                chain_name=chain_name,
                active_site=active_site,
                mass=sum(cg_bead_mass[(name, degree)].values()))

        # Prepares topological information.
        cg_bonded = {
            'bond': {},
            'angle': {},
            'dihedral': {}
        }
        num = {'bond': 2, 'angle': 3, 'dihedral': 4}
        for name in cg_bonded:
            for x in root.iter(name):
                term_name = x.find('name').text.strip()
                bonds = x.find('beads').text.strip().split()
                params = x.find('params')
                if params is not None:
                    params = params.text.strip()
                tuple_size = num[name]
                cg_bonded[name][term_name] = {
                    'list': [bonds[i:i+tuple_size] for i in range(0, len(bonds), tuple_size)],
                    'params': params
                }
        return XMLMap(
            ident=ident_name,
            name=chain_name,
            mass_map=cg_bead_mass,
            source_coordinates=source_file,
            source_topology=source_topology,
            molecule_beads=cg_beads,
            molecule_topology=cg_bonded)

    def _parse_hybrid_topology(self):
        cg_bonded = {
            'hybrid_bonds': {},
            'output_type': 'single',
            'file': None,
            'include': ''
        }
        root = self.root.find('hybrid_topology')
        cg_bonded['file'] = root.find('file').text.strip()
        cg_bonded['moleculetype'] = {
            'name': root.find('molecule_type').find('name').text,
            'nrexcl': root.find('molecule_type').find('exclusion').text
            }
        cg_bonded['system'] = root.find('system').text
        include_section = root.find('include')
        if include_section is not None:
            include_text = '\n'.join(map(str.strip, include_section.text.strip().split('\n')))
            cg_bonded['include'] = include_text + '\n\n'

        for x in root.findall('bond'):
            beads = x.find('beads').text.strip().split()
            bonds = [beads[i:i+2] for i in range(0, len(beads), 2)]
            bond_params = x.find('bond_params')
            if bond_params is not None:
                bond_params = bond_params.text.strip().split()
            angle_params = x.find('angle_params')
            if angle_params is not None:
                angle_params = angle_params.text.strip().split()
            dih_params = x.find('dihedral_params')
            if dih_params is not None:
                dih_params = dih_params.text.strip().split()
            pair_params = x.find('pair_params')
            if pair_params is not None:
                pair_params = pair_params.text.strip().split()
            for b in bonds:
                cg_bonded['hybrid_bonds'][tuple(b)] = {
                    'bond_params': bond_params,
                    'angle_params': angle_params,
                    'dihedral_params': dih_params,
                    'pair_params': pair_params
                }
        return cg_bonded
