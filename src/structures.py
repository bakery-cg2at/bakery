"""
Copyright (C) 2015-2016 Jakub Krajniak <jkrajniak@gmail.com>

This file is part of Backmapper.

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
import copy
import files_io
import tools
import xml.etree.ElementTree as etree

import numpy as np
import random
import networkx
import logging
import sys

__doc__ = "Data structures."""

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

CGMolecule = collections.namedtuple(
    'CGMolecule',
    ['name', 'ident', 'source_file', 'source_topology']
)

logger = logging.getLogger()

class CGFragment:
    def __init__(self, cg_molecule, fragment_list, active_sites=None, charge_map=None, as_remove=None):
        self.topology = cg_molecule.source_topology
        self.coordinate = cg_molecule.source_file
        self.fragment_list = fragment_list

        self.cg_molecule = cg_molecule

        self.charge_map = charge_map
        self.active_sites_remove_map = as_remove

        # Select the set of atomistic coordinates.
        chains = self.coordinate.chains[cg_molecule.ident]
        # Select random chain from the ensemble of chains.
        atoms = chains[random.sample(chains.keys(), 1)[0]]
        self.atomparams = {
            k: v[0] for k, v in self.topology.chain_atom_names[cg_molecule.ident].items()
            }
        self.active_sites = {}

        self.com = np.zeros(3)
        total_mass = 0.0
        tmp_atom_list = []
        for bead in self.fragment_list:
            atom_name = bead.split(':')[2]
            atom = atoms[atom_name]._replace(chain_name=cg_molecule.ident)
            atom_mass = self.atomparams[atom_name].mass
            self.com += atom_mass * atom.position
            total_mass += atom_mass
            tmp_atom_list.append(atom)

        self.com /= total_mass
        self.cg_mass = total_mass
        # Move the atoms in fragments to the origin, always by substracting the com of the fragment.
        self.atom_in_fragments = [
            files_io.Atom(x[0], x[1], x[2], x[3], x[4] - self.com)
            for x in tmp_atom_list
            ]
        if active_sites:
            for at_as in active_sites:
                t = at_as.split(':')
                self.active_sites[t[1]] = int(t[2])


class BackmapperSettings2:
    def __init__(self, input_xml):
        self.res2atom = collections.defaultdict(list)
        self.cg_active_sites = collections.defaultdict(list)
        self.atom_id2fragment = {}
        self.atom_ids = []
        self.cg_old_new_id = {}
        self.cross_prefix = 'cross_'
        self.cg_new_id_old = {}
        self.atom2cg = {}
        self.fragments = {}  # key-> molecule name, val-> dict
        self.mol_atomid_map = collections.defaultdict(dict)  # Map between old id and new id, divided into molecules.
        self.mol_atomname_map = collections.defaultdict(dict)  # Map between old id and new id, divided into molecules.
        self.cg2atom = collections.defaultdict(list)

        self.charge_transfer = {}   # Map with charge transfer.

        self.global_graph = networkx.Graph()

        tree = etree.parse(input_xml)
        self.root = tree.getroot()
        self._parse()

    def _parse(self):
        cg_configuration = self.root.find('cg_configuration')
        self.cg_topology = files_io.GROMACSTopologyFile(
            cg_configuration.find('topology').text.strip())
        self.cg_topology.read()

        self.cg_topology.replicate()

        self.cg_graph = self.cg_topology.get_graph()
        self.cg_coordinate = files_io.read_coordinates(
            cg_configuration.find('file').text.strip()
        )

        hyb_configuration = self.root.find('hybrid_configuration')
        self.hybrid_configuration = {
            'format': hyb_configuration.find('format').text.strip(),
            'file': files_io.GROFile(hyb_configuration.find('file').text.strip())
        }

        hyb_topology = self.root.find('hybrid_topology')
        self.hyb_topology = files_io.GROMACSTopologyFile(hyb_topology.find('file').text.strip())
        self.hyb_topology.init()
        self.hyb_topology.new_data.update({
            '{}{}'.format(self.cross_prefix, k): {} for k in self.hyb_topology.new_data
        })
        self.hyb_topology.header_section.append('; GROMACS like hybrid topology file.\n; ===================== \n\n')
        if hyb_topology.find('include') is not None:
            for include_entry in hyb_topology.find('include').text.split():
                ie = include_entry.strip()
                if ie.startswith(';'):
                    self.hyb_topology.header_section.append(';#include {}\n'.format(ie.replace(';', '')))
                else:
                    self.hyb_topology.header_section.append('#include {}\n'.format(ie))
            self.hyb_topology.header_section.append('\n')

        self.hyb_topology.atomtypes = {}
        # Change atomtypes to 'V' for CG.
        for k, d in self.cg_topology.atomtypes.items():
            d['type'] = 'V'
            self.hyb_topology.atomtypes[k] = d
        self.hyb_topology.moleculetype = {
            'name': hyb_topology.find('molecule_type').find('name').text.strip(),
            'nrexcl': hyb_topology.find('molecule_type').find('exclusion').text.strip()
        }
        self.hyb_topology.molecules = {
            'name': hyb_topology.find('system').text.strip(),
            'mol': 1
        }
        self.hyb_topology.bondtypes = self.cg_topology.bondtypes
        self.hyb_topology.angletypes = self.cg_topology.angletypes
        self.hyb_topology.dihedraltypes = self.cg_topology.dihedraltypes

        # Reads parameters of bonded terms.
        self.bond_params = collections.defaultdict(dict)
        self.angle_params = collections.defaultdict(dict)
        self.dihedral_params = collections.defaultdict(dict)
        bonded_terms = hyb_topology.findall('bonds')
        for bond_term in bonded_terms:
            params = bond_term.attrib['params'].split()
            bead_list = bond_term.text.strip().split()
            if len(bead_list) % 2 != 0:
                raise RuntimeError('Wrong pairs in <bonds> section, found {}'.format(bead_list))
            pair_list = zip(bead_list[::2], bead_list[1::2])
            for p1, p2 in pair_list:
                self.bond_params[p1][p2] = params
                self.bond_params[p2][p1] = params
        angle_terms = hyb_topology.findall('angles')
        for angle_term in angle_terms:
            params = angle_term.attrib['params'].split()
            bead_list = angle_term.text.strip().split()
            if len(bead_list) % 3 != 0:
                raise RuntimeError('Wrong triplets in <angles> section, found {}'.format(bead_list))
            for idx in range(0, len(bead_list)-3, 3):
                p1, p2, p3 = bead_list[idx], bead_list[idx+1], bead_list[idx+2]
                if p1 not in self.angle_params:
                    self.angle_params[p1] = {}
                if p3 not in self.angle_params:
                    self.angle_params[p3] = {}
                if p2 not in self.angle_params[p1]:
                    self.angle_params[p1][p2] = {}
                if p2 not in self.angle_params[p3]:
                    self.angle_params[p3][p2] = {}
                self.angle_params[p1][p2][p3] = params
                self.angle_params[p3][p2][p1] = params
        dihedral_terms = hyb_topology.findall('dihedrals')
        for dihedral_term in dihedral_terms:
            params = dihedral_term.attrib['params'].split()
            bead_list = dihedral_term.text.strip().split()
            if len(bead_list) % 4 != 0:
                raise RuntimeError('Wrong quadruplets in <dihedrals> section, found {}'.format(bead_list))
            for idx in range(0, len(bead_list) - 4, 4):
                p1, p2, p3, p4 = bead_list[idx], bead_list[idx + 1], bead_list[idx + 2], bead_list[idx + 3]
                if p2 not in self.dihedral_params[p1]:
                    self.dihedral_params[p1][p2] = {}
                if p3 not in self.dihedral_params[p1][p2]:
                    self.dihedral_params[p1][p2][p3] = {}
                if p4 in self.dihedral_params[p1][p2][p3]:
                    raise RuntimeError('Parameters for triplet: {} already defined {}'.format(
                        (p1, p2, p3, p4), params
                    ))
                self.dihedral_params[p1][p2][p3][p4] = params
                if p3 not in self.dihedral_params[p4]:
                    self.dihedral_params[p4][p3] = {}
                if p2 not in self.dihedral_params[p4][p3]:
                    self.dihedral_params[p4][p3][p2] = {}
                if p1 in self.dihedral_params[p4][p3][p2]:
                    raise RuntimeError('Parameters for triplet: {} already defined {}'.format(
                        (p4, p3, p2, p1), params
                    ))
                self.dihedral_params[p4][p3][p2][p1] = params

        # For testing purposes, write the cg_topol.top as a graph structure.
        networkx.write_gpickle(self.cg_graph, 'cg_topol.top.pck')
        self._parse_cg_molecules()

    def _parse_cg_molecules(self):
        """Parse the cg_molecule section and prepares data structures."""
        for r in self.root.findall('cg_molecule'):
            cg_molecule = CGMolecule(
                r.find('name').text,
                r.find('ident').text,
                files_io.read_coordinates(r.find('source_file').text),
                files_io.GROMACSTopologyFile(r.find('source_topology').text))
            cg_molecule.source_topology.read()
            self.fragments[cg_molecule.name] = {
                'cg_molecule': cg_molecule
            }
            # Generate the beads.
            for cg_bead in r.findall('cg_bead'):
                name = cg_bead.find('name').text
                self.fragments[cg_molecule.name][name] = {}
                for beads in cg_bead.findall('beads'):
                    degree = beads.attrib.get('degree', '*')
                    active_sites = beads.attrib.get('active_site')
                    as_remove = {}
                    if active_sites:
                        active_sites = map(str.strip, active_sites.split())
                        tmp_as = [':'.join(x.split(':')[:2]) for x in active_sites]
                        # Remove with active sites. So whenever the active site will be used then this list
                        # of atoms will be removed.
                        removes = beads.findall('remove')
                        for asr in removes:
                            as_name = asr.attrib['active_site']
                            if as_name not in tmp_as:
                                raise RuntimeError(
                                    'Defined <remove> entry for active_site {} but it\'s not on the list'.format(
                                        as_name))
                            beads_to_remove = asr.text
                            if not beads_to_remove:
                                raise RuntimeError('No beads to remove in active_site {}'.format(as_name))
                            as_remove[as_name] = map(str.strip, asr.text.split())
                    bead_list = beads.text.split()

                    cm = beads.find('charge_map')
                    charge_map = None
                    if cm is not None:
                        charge_map = cm.text.split()
                        if len(charge_map) != len(bead_list):
                            raise RuntimeError(
                                'Number of entries in charge_map {} does not match number of beads {}'.format(
                                    len(charge_map), len(bead_list)))
                    cg_fragment = CGFragment(cg_molecule, bead_list, active_sites, charge_map, as_remove)
                    self.fragments[cg_molecule.name][name][degree] = cg_fragment

            # Read charge transfer.
            # <charge_transfer on="IPD:N1:2" from="IPD:H8" to="EPO:C23#H25,EPO:C41#H66,HDD:C21#H43,HDD:C32#H44" />
            charge_transfers = r.find('charge_transfers')
            if charge_transfers is not None:
                for ct in charge_transfers.findall('charge_transfer'):
                    t = ct.attrib['on'].split(':')
                    on_key = '{}:{}'.format(t[0], t[1])
                    on_deg = int(t[2])
                    if on_key not in self.charge_transfer:
                        self.charge_transfer[on_key] = {}
                    if on_deg in self.charge_transfer[on_key]:
                        raise RuntimeError('Mistake in charge transfer, key {} already defined'.format(t))
                    transfer_from = ct.attrib['from']
                    transfer_to_map = {}
                    self.charge_transfer[on_key][on_deg] = {'from': transfer_from.split(':'), 'to': transfer_to_map}
                    transfer_to = ct.attrib['to'].split(',')
                    print transfer_to
                    for tt in transfer_to:
                        tt_to_on, tt_to = tt.split('#')
                        transfer_to_map[tt_to_on] = tt_to


    def prepare_hybrid(self):
        """Creates hybrid files."""
        outfile = self.hybrid_configuration['file']
        new_at_id = 1
        new_res_id = 0
        last_res_id = 0
        cg_ids = sorted(self.cg_graph.nodes())

        for cg_id in cg_ids:
            cg_atom = self.cg_coordinate.atoms[cg_id]
            cg_bead = self.cg_graph.node[cg_id]
            cg_bead_id = new_at_id
            if cg_bead['res_id'] != last_res_id:
                last_res_id = cg_bead['res_id']
                new_res_id += 1

            cg_degree = str(cg_bead['degree'])
            fragments = self.fragments[cg_bead['chain_name']][cg_bead['name']]
            # Choose based on the degree of node.
            cg_fragment = fragments.get(cg_degree, fragments.get('*'))
            if not cg_fragment:
                raise RuntimeError('Problem with getting atomistic fragments')

            outfile.atoms[cg_bead_id] = files_io.Atom(
                atom_id=cg_bead_id,
                name=cg_bead['name'],
                chain_name=cg_fragment.cg_molecule.ident, # Change chain name to one from <ident> tag in the cg_molecule
                chain_idx=new_res_id,
                position=cg_atom.position
            )
            self.cg_new_id_old[cg_bead_id] = cg_id
            self.cg_old_new_id[cg_id] = cg_bead_id

            self.atom_id2fragment[new_at_id] = cg_fragment

            self.global_graph.add_node(cg_bead_id, **cg_bead)

            # Change the mass of CG bead
            topol_atom = copy.copy(self.cg_topology.atoms[cg_id])
            topol_atom.mass = cg_fragment.cg_mass
            topol_atom.atom_id = cg_bead_id
            topol_atom.chain_idx = new_res_id
            topol_atom.chain_name = cg_fragment.cg_molecule.ident # Change chain name to one from <ident> tag in the cg_molecule
            topol_atom.cgnr = cg_bead_id
            self.hyb_topology.atoms[new_at_id] = topol_atom
            if cg_fragment.cg_molecule.ident not in self.hyb_topology.chains:
                self.hyb_topology.chains[cg_fragment.cg_molecule.ident] = {}
            if new_res_id not in self.hyb_topology.chains[cg_fragment.cg_molecule.ident]:
                self.hyb_topology.chains[cg_fragment.cg_molecule.ident][new_res_id] = {}
            if topol_atom.name in self.hyb_topology.chains[cg_fragment.cg_molecule.ident][new_res_id]:
                raise RuntimeError(
                    '{} already defined, please make sure that atom names are unique'.format(topol_atom.name))
            self.hyb_topology.chains[cg_fragment.cg_molecule.ident][new_res_id][topol_atom.name] = topol_atom

            # Set the atomistic coordinates for this fragment, put it in the
            # output coordinate file.
            cg_com = cg_atom.position
            new_at_id += 1
            for idx, at in enumerate(cg_fragment.atom_in_fragments):
                new_at_atom = files_io.Atom(
                    new_at_id,
                    at.name,
                    cg_fragment.cg_molecule.ident,  # Change chain name to one from <ident> tag in the cg_molecule
                    new_res_id,
                    at.position + cg_com
                )
                outfile.atoms[new_at_id] = new_at_atom
                self.global_graph.add_node(
                    new_at_id,
                    name=at.name,
                    res_id=new_res_id,
                    position=at.position + cg_com,
                    chain_name=at.chain_name)

                # Set topology atom
                topol_atom = copy.copy(cg_fragment.topology.chain_atom_names[at.chain_name][at.name][0])
                if new_res_id not in self.mol_atomid_map[at.chain_name]:
                    self.mol_atomid_map[at.chain_name][new_res_id] = {}
                    self.mol_atomname_map[at.chain_name][new_res_id] = {}
                self.mol_atomid_map[at.chain_name][new_res_id][topol_atom.atom_id] = new_at_id
                self.mol_atomname_map[at.chain_name][new_res_id][topol_atom.name] = new_at_id

                if cg_fragment.charge_map:
                    topol_atom.charge = cg_fragment.charge_map[idx]
                topol_atom.atom_id = new_at_id
                topol_atom.chain_idx = new_res_id
                topol_atom.chain_name = cg_fragment.cg_molecule.ident
                topol_atom.cgnr = new_at_id
                self.hyb_topology.atoms[new_at_id] = topol_atom
                if cg_fragment.cg_molecule.ident not in self.hyb_topology.chains:
                    self.hyb_topology.chains[cg_fragment.cg_molecule.ident] = {}
                if new_res_id not in self.hyb_topology.chains[cg_fragment.cg_molecule.ident]:
                    self.hyb_topology.chains[cg_fragment.cg_molecule.ident][new_res_id] = {}
                if topol_atom.name in self.hyb_topology.chains[cg_fragment.cg_molecule.ident][new_res_id]:
                    raise RuntimeError(
                        '{} already defined, please make sure that atom names are unique'.format(topol_atom.name))
                self.hyb_topology.chains[cg_fragment.cg_molecule.ident][new_res_id][topol_atom.name] = topol_atom

                # Set active sites.
                if at.name in cg_fragment.active_sites:
                    self.cg_active_sites[cg_bead_id].append((new_at_atom, cg_fragment.active_sites[at.name]))

                self.atom_id2fragment[new_at_id] = cg_fragment
                self.atom_ids.append(new_at_id)
                self.atom2cg[new_at_id] = cg_bead_id
                self.cg2atom[cg_bead_id].append(new_at_id)
                self.res2atom[new_res_id].append(new_at_id)

                new_at_id += 1

        outfile.box = self.cg_coordinate.box

        # Rebuild hybrid topology.
        self.rebuild_hybrid_topology()
        self.hyb_topology.write()
        # Write the hybrid coordinate file.
        outfile.write(force=True)

    def rebuild_hybrid_topology(self):
        """Regenerate the hybrid topology based on the new particle ids."""
        # First build coarse-grained topology.
        def generate_cg_b_list(old_list):
            if not old_list:
                return {}
            return {tuple((map(self.cg_old_new_id.get, k))): v + [' ; cg_bonded'] for k, v in old_list.items()}

        def generate_at_b_list(at):
            """Generates atomistic list of bonds"""
            topology = self.atom_id2fragment[at.atom_id].topology
            bonded_lists = [
                (2, topology.bonds, 'bonds'),
                (3, topology.angles, 'angles'),
                (4, topology.dihedrals, 'dihedrals'),
                (4, topology.improper_dihedrals, 'dihedrals'),
                (2, topology.pairs, 'pairs')
            ]
            for ternary, top_list, output_name in bonded_lists:
                for p, params in top_list.items():
                    old2new_id = self.mol_atomid_map[at.chain_name][at.chain_idx]
                    new_tuple = tuple(map(old2new_id.get, p))

                    # Skip the terms that involves the missing atoms. The atoms can be missing because of
                    # the degree dependent atomistic fragments.
                    if None not in new_tuple and at.atom_id in new_tuple:  # correct pair
                        cg_beads = map(self.atom2cg.get, new_tuple)
                        if cg_beads.count(cg_beads[0]) == len(cg_beads):  # Atoms in the same CG bead
                            out_name = output_name
                        else:
                            out_name = '{}{}'.format(self.cross_prefix, output_name)
                        self.hyb_topology.new_data[out_name][new_tuple] = params

        print('Renumering atomistic and coarse-grained bonds')
        self.hyb_topology.new_data['{}bonds'.format(self.cross_prefix)] = generate_cg_b_list(
            self.cg_topology.bonds)
        self.hyb_topology.new_data['{}angles'.format(self.cross_prefix)] = generate_cg_b_list(
            self.cg_topology.angles)
        self.hyb_topology.new_data['{}dihedrals'.format(self.cross_prefix)] = generate_cg_b_list(
            self.cg_topology.dihedrals)
        self.hyb_topology.new_data['{}dihedrals'.format(self.cross_prefix)].update(
            generate_cg_b_list(self.cg_topology.improper_dihedrals))

        # Store exclusion list based on the boned-terms from renumerated cg topology
        ex_cg_list = open('exclusion_list_cg.dat', 'w')
        for b in self.hyb_topology.new_data['{}bonds'.format(self.cross_prefix)]:
            ex_cg_list.write('{} {}\n'.format(*b))
        for a in self.hyb_topology.new_data['{}angles'.format(self.cross_prefix)]:
            ex_cg_list.write('{} {}\n'.format(a[0], a[2]))
        for d in self.hyb_topology.new_data['{}dihedrals'.format(self.cross_prefix)]:
            ex_cg_list.write('{} {}\n'.format(d[0], d[2]))
            ex_cg_list.write('{} {}\n'.format(d[1], d[3]))
            ex_cg_list.write('{} {}\n'.format(d[0], d[3]))
        ex_cg_list.close()

        # Create the atomistic topology.
        for atid in self.atom_ids:
            generate_at_b_list(self.hyb_topology.atoms[atid])

        for (b1, b2), params in self.hyb_topology.new_data['{}bonds'.format(self.cross_prefix)].items():
            self.global_graph.add_edge(b1, b2)
        for (b1, b2), params in self.hyb_topology.new_data['bonds'].items():
            self.global_graph.add_edge(b1, b2)

        # Generate list of CG bonds that are not defined at AT level
        cg_cross_bonds = set()
        for b1, b2 in self.global_graph.edges():
            if b1 in self.cg_new_id_old and b2 in self.cg_new_id_old:  # Look only on CG bonds
                ats_b1 = self.cg2atom[b1]
                ats_b2 = set(self.cg2atom[b2])
                is_connected = False
                for at1 in ats_b1:
                    edges = set(self.global_graph.edge[at1])
                    intersect = ats_b2.intersection(edges)
                    if intersect:
                        is_connected = True
                        break
                if not is_connected:
                    cg_cross_bonds.add(tuple(sorted([b1, b2])))

        print('Found {} bonds between coarse-grained beads'.format(len(cg_cross_bonds)))
        print('Generating atomistic cross-bonds between coarse-grained beads; It will take a while...')
        # Create the atomistic bonds across the coarse-grained beads.
        # We iterate over bonds in cg_graph and then generate the bonds.
        # At the level of topology, CG bonds are already defined.
        at_cross_bonds = []
        atoms_to_remove = []
        charge_to_transfer = []
        progress_indc = 0.0
        progress_indc_total = len(cg_cross_bonds)
        global_degree = self.global_graph.degree()
        for b1, b2 in cg_cross_bonds:
            n1 = self.global_graph.node[b1]
            n2 = self.global_graph.node[b2]
            if n1['res_id'] != n2['res_id']:  # Cross bond between beads in different chains.
                # Look for active sites on both CG molecules.
                ats1, ats2 = None, None
                for at1, max_d1 in self.cg_active_sites[b1]:
                    for at2, max_d2 in self.cg_active_sites[b2]:
                        b1_key = '{}:{}'.format(at1.chain_name, at1.name)
                        b2_key = '{}:{}'.format(at2.chain_name, at2.name)
                        test_bond = self.bond_params.get(b1_key, {}).get(b2_key, None)
                        if (test_bond and self.global_graph.has_node(at1.atom_id)
                            and self.global_graph.has_node(at2.atom_id)):
                            at1_deg = global_degree[at1.atom_id]
                            at2_deg = global_degree[at2.atom_id]
                            at_remove1 = self.atom_id2fragment[b1].active_sites_remove_map.get(b1_key)
                            at_remove2 = self.atom_id2fragment[b2].active_sites_remove_map.get(b2_key)
                            tmp_atoms_to_remove = []
                            valid = True
                            # Before we check the degree, let's first try to remove atoms (if they are defined to
                            # remove and after that compare the degree.
                            if at_remove1:
                                for atr1 in at_remove1:
                                    atr_chain_name, atr_name = atr1.split(':')[1], atr1.split(':')[2]
                                    atr1_id = self.mol_atomname_map[atr_chain_name][at1.chain_idx][atr_name]
                                    if not self.global_graph.has_node(atr1_id):
                                        valid = False
                                        break
                                    tmp_atoms_to_remove.append(atr1_id)
                                    if self.global_graph.has_edge(atr1_id, at1.atom_id):
                                        at1_deg -= 1
                            if at_remove2 and valid:
                                for atr2 in at_remove2:
                                    atr_chain_name, atr_name = atr2.split(':')[1], atr2.split(':')[2]
                                    atr2_id = self.mol_atomname_map[atr_chain_name][at2.chain_idx][atr_name]
                                    if not self.global_graph.has_node(atr2_id):
                                        valid = False
                                        break
                                    tmp_atoms_to_remove.append(atr2_id)
                                    if self.global_graph.has_edge(atr2_id, at2.atom_id):
                                        at2_deg -= 1

                            # Check the degree after update with virtual removing of atoms.
                            if at1_deg < max_d1 and at2_deg < max_d2 and valid:
                                # Found correct pair of active sites. Remove the atoms that were
                                # connected to this atom if in settings the set of atoms to remove were defined.
                                ats1, ats2 = at1, at2   # Pair of selected atoms.
                                atoms_to_remove.extend(tmp_atoms_to_remove)
                                for ai in tmp_atoms_to_remove:
                                    self.global_graph.remove_node(ai)
                                global_degree = self.global_graph.degree()
                                break
                    if ats1 is not None and ats2 is not None:
                        break

                if ats1 is not None and ats2 is not None:
                    # Selected proper active sites.
                    at_cross_bonds.append((ats1.atom_id, ats2.atom_id))
                    self.global_graph.add_edge(ats1.atom_id, ats2.atom_id)
                    global_degree[ats1.atom_id] += 1
                    global_degree[ats2.atom_id] += 1
                    # Look for charge to transfer.
                    deg1 = global_degree[ats1.atom_id]
                    deg2 = global_degree[ats2.atom_id]
                    b_key1 = '{}:{}'.format(ats1.chain_name, ats1.name)
                    b_key2 = '{}:{}'.format(ats2.chain_name, ats2.name)
                    charge_transfer_cfg = None
                    target_atom = None
                    target_atom_key = None
                    if b_key1 in self.charge_transfer:
                        charge_transfer_cfg = self.charge_transfer[b_key1][deg1]
                        target_atom_key = b_key2
                        target_atom = ats2
                    elif b_key2 in self.charge_transfer:
                        charge_transfer_cfg = self.charge_transfer[b_key2][deg2]
                        target_atom_key = b_key1
                        target_atom = ats1

                    if charge_transfer_cfg is not None:
                        at_to_transfer = charge_transfer_cfg['to'][target_atom_key]
                        # Gets the charge from the source topology.
                        from_chain_name, from_atom_name = charge_transfer_cfg['from']
                        at_from_transfer = (
                            self.fragments[from_chain_name]['cg_molecule']
                                .source_topology.chain_atom_names[from_chain_name][from_atom_name][0])
                        # Now update directly hybrid topology with new charge.
                        self.hyb_topology.chains[target_atom.chain_name][target_atom.chain_idx][at_to_transfer].charge = at_from_transfer.charge
                else:
                    print b1, b2, n1, n2
                    print self.cg_active_sites[b1]
                    print self.cg_active_sites[b2]
                    raise RuntimeError('Something is really wrong!')
            sys.stdout.write('{} %\r'.format(100.0*(progress_indc/progress_indc_total)))
            progress_indc += 1.0

        # Generate entries for AT cross bonds.
        print('Found {} atomistic cross bonds'.format(len(at_cross_bonds)))
        self._generate_atomistic_bonds(at_cross_bonds)
        self._remove_atomistic_particles(set(atoms_to_remove))
        # Generate exclusion list.
        self._generate_exclusion_lists()

    def _generate_atomistic_bonds(self, at_cross_bonds):
        """Generates parameters for atomistic bonds."""
        fout_filename = 'missing_definitions.txt'
        fout = open(fout_filename, 'w')
        missing_definitions = set()

        bond_key = '{}bonds'.format(self.cross_prefix)
        angle_key = '{}angles'.format(self.cross_prefix)
        dihedral_key = '{}dihedrals'.format(self.cross_prefix)

        for b1, b2 in at_cross_bonds:
            n1 = self.global_graph.node[b1]
            n2 = self.global_graph.node[b2]
            b1_key = n1_key = '{}:{}'.format(n1['chain_name'], n1['name'])
            b2_key = n2_key = '{}:{}'.format(n2['chain_name'], n2['name'])
            if ((b1, b2) in self.hyb_topology.new_data[bond_key] or
                    (b2, b2) in self.hyb_topology.new_data[bond_key]):
                raise RuntimeError('Bond {}-{} already defined in the topology, wtf?'.format(b1, b2))
            try:
                self.hyb_topology.new_data[bond_key][(b1, b2)] = self.bond_params[n1_key][n2_key]
            except KeyError as ex:
                print('Missing definition of bond: {}-{}'.format(n1_key, n2_key))
                raise ex
            # Generate angles.
            triplets = tools.gen_bonded_tuples(self.global_graph, 3, (b1, b2))
            quadruplets = tools.gen_bonded_tuples(self.global_graph, 4, (b1, b2))
            for triplet in triplets:
                n1_key, n2_key, n3_key = [
                    '{}:{}'.format(x['chain_name'], x['name'])
                    for x in map(self.global_graph.node.get, triplet)]
                try:
                    if (triplet in self.hyb_topology.new_data[angle_key] or
                            tuple(reversed(triplet)) in self.hyb_topology.new_data[angle_key]):
                        raise RuntimeError('Angle {} already defined in the topology, wtf?'.format(triplet))
                    self.hyb_topology.new_data[angle_key][triplet] = (
                        self.angle_params[n1_key][n2_key][n3_key])
                except KeyError as ex:
                    print('Missing definition of angle: {}-{}-{} (bond {}-{})'.format(
                        n1_key, n2_key, n3_key, b1_key, b2_key))
                    missing_definitions.add((n1_key, n2_key, n3_key))
            # Generate dihedrals.
            for quadruplet in quadruplets:
                n1_key, n2_key, n3_key, n4_key = [
                    '{}:{}'.format(x['chain_name'], x['name'])
                    for x in map(self.global_graph.node.get, quadruplet)]
                try:
                    if (quadruplet in self.hyb_topology.new_data[dihedral_key] or
                            tuple(reversed(quadruplet)) in self.hyb_topology.new_data[dihedral_key]):
                        raise RuntimeError('Dihedral {} already defined in the topology, wtf?'.format(quadruplet))
                    self.hyb_topology.new_data[dihedral_key][quadruplet] = (
                        self.dihedral_params[n1_key][n2_key][n3_key][n4_key])
                except KeyError:
                    print('Missing definition of dihedral: {}-{}-{}-{} (bond {}-{})'.format(
                        n1_key, n2_key, n3_key, n4_key, b1_key, b2_key))
                    missing_definitions.add((n1_key, n2_key, n3_key, n4_key))
        fout.writelines('\n'.join([' '.join(x) for x in sorted(missing_definitions)]))
        print('Wrote missing definitions in {}'.format(fout_filename))
        fout.close()

    def _remove_atomistic_particles(self, atoms_to_remove):
        """Update coordinate and topology file by removing atoms and renumbering"""
        if atoms_to_remove:
            print('Clean up atomistic particles after creating bonds, atoms to remove: {}'.format(len(atoms_to_remove)))
            self.hyb_topology.remove_atoms(atoms_to_remove, renumber=False)
            self.hybrid_configuration['file'].remove_atoms(atoms_to_remove, renumber=False)
            # Renumber data files
            self.hyb_topology.renumber()
            self.hybrid_configuration['file'].renumber()

    def _generate_exclusion_lists(self):
        # Collect all bonds and greate global_graph again...
        print('Generating exclusion lists, nrexcl={}'.format(self.hyb_topology.moleculetype['nrexcl']))
        g = networkx.Graph()
        bonds = [b for b in self.hyb_topology.bonds.keys()]
        for k in self.hyb_topology.new_data:
            if 'bonds' in k:
                bonds.extend([b for b in self.hyb_topology.new_data[k]])
        g.add_edges_from(bonds)
        paths = networkx.all_pairs_shortest_path(g, int(self.hyb_topology.moleculetype['nrexcl']))
        exclusions = set()
        for l in paths.values():
            for p in l.values():
                if len(p) > 1:
                    exclusions.add(tuple(sorted([p[0], p[-1]])))
        output_filename = 'exclusion_{}.list'.format(self.hyb_topology.file_name.split('.')[0])
        out_file = open(output_filename, 'w')
        out_file.writelines('\n'.join(['{} {}'.format(*d) for d in sorted(exclusions)]))
        out_file.close()
        print('Generated {} exclusions, writen to {}'.format(len(exclusions), output_filename))