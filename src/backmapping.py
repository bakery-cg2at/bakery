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
import copy
from . import files_io
import numpy as np
from . import structures
from . import tools

__doc__ = 'Set of functions that are use for backmapping.'


def calculate_com_fragments(backmapper_settings):
    """For each of the bead calculate the center of mass based on the all-atom coordinate file.

    Args:
        backmapper_settings: The parsed configuration file.

    Returns:
        The center of mass of each of CG beads for molecules and
        position of atomistic fragments translated so that COM of fragment matches COM of CG bead.
    """
    print('Calculate COM of atomistic fragments')
    cg_com = {}
    cg_aa_beads = {}
    for mol_name, cg_mol in backmapper_settings.cg_molecules.items():
        aa_file = files_io.GROFile(cg_mol.source_coordinates)
        aa_file.read()
        cg_com_coords = {}
        cg_aa = collections.defaultdict(list)
        for cg_bead, aa_beads in cg_mol.mass_map.items():
            c_mp = [0.0, 0.0, 0.0]
            tot_mass = 0.0
            for aa_bead, aa_mass in aa_beads.items():
                chain_id, chain_name, atom_name = aa_bead.split(':')
                pos = aa_file.fragments[chain_name][atom_name].position
                tot_mass += aa_mass
                cg_aa[cg_bead].append(aa_file.fragments[chain_name][atom_name])
                for i in range(3):
                    c_mp[i] += pos[i] * aa_mass
            for i in range(3):
                c_mp[i] /= tot_mass
            cg_com_coords[cg_bead] = np.array(c_mp)

        for ca in list(cg_aa.values()):
            ca.sort(key=lambda x: x.atom_id)
        cg_com[cg_mol.name] = cg_com_coords
        cg_aa_beads[cg_mol.name] = cg_aa
    return cg_com, cg_aa_beads


def prepare_hyb_coordinates(settings, cg_graph, cg_com, cg_aa, plain=False):
    """Prepares hybrid coordinates file and topology.

    Args:
        settings: The BackmapperSettings object.
        cg_graph: The networkx.Graph object with the CG structure.
        cg_com: Indexed by molecule name and bead index, position of COM of CG beads.
        cg_aa: Dictionary with atomistic fragments.
        plain: If set to 'CG' or 'AA' then prepare only configuration for that kind
            of systems.

    Returns:
        hybrid coordinate file object and topology object for further process.
    """
    print('Preparing hybrid configuration files')
    cg_ids = sorted(cg_graph.nodes())
    hyb_file = files_io.GROFile(settings.hybrid_configuration['file'])
    input_aa_top = {}
    atom_id_old2new = {}
    # Reads AT topologies
    for mol_name, cg_mol in settings.cg_molecules.items():
        input_aa_top[mol_name] = files_io.GROMACSTopologyFile(cg_mol.source_topology)
        input_aa_top[mol_name].read()
        atom_id_old2new[mol_name] = collections.defaultdict(dict)

    output_topology = files_io.GROMACSTopologyFile(settings.hybrid_topology['file'])
    atom_types = {}
    atom_id2cg_bead_id = {}
    cg_old_new_id = {}
    chain2atoms = collections.defaultdict(dict)

    new_at_id = 1
    if plain == 'AA':
        for cg_id in cg_ids:
            cg_bead = cg_graph.node[cg_id]
            cg_bead_id = cg_id
            name = cg_bead['chain_name']
            cg_key = structures.BeadID(cg_bead['name'], cg_bead['degree'])
            if cg_key not in settings.cg_molecules[name].molecule_beads:
                cg_key = structures.BeadID(cg_bead['name'], '*')
            for aa_bead in cg_aa[name][cg_key]:
                hyb_file.atoms[new_at_id] = files_io.Atom(
                    atom_id=new_at_id,
                    name=aa_bead.name,
                    chain_name=aa_bead.chain_name,
                    chain_idx=cg_bead['res_id'],
                    position=(aa_bead.position +
                              (np.array(cg_bead['position']) -
                               np.array(cg_com[cg_bead['chain_name']][cg_key]))))
                topo_at = copy.copy(
                    input_aa_top[aa_bead.chain_name]
                    .chain_atom_names[aa_bead.chain_name][aa_bead.name][0])
                atom_id_old2new[aa_bead.chain_name][cg_bead['res_id']][topo_at.atom_id] = new_at_id
                topo_at.chain_idx = cg_bead['res_id']
                topo_at.cgnr = cg_bead['res_id']
                topo_at.atom_id = new_at_id
                output_topology.atoms[new_at_id] = topo_at
                atom_id2cg_bead_id[new_at_id] = cg_bead_id
                chain2atoms[cg_bead['res_id']]['{}:{}'.format(
                    aa_bead.chain_name, aa_bead.name)] = new_at_id
                new_at_id += 1
    elif plain == 'CG':
        for cg_id in cg_ids:
            cg_bead = cg_graph.node[cg_id]
            cg_bead_id = new_at_id
            name = cg_bead['chain_name']
            hyb_file.atoms[new_at_id] = files_io.Atom(
                atom_id=cg_bead_id,
                name=cg_bead['name'],
                chain_name=cg_bead['chain_name'],
                chain_idx=cg_bead['res_id'],
                position=cg_bead['position'])
            try:
                cg_key = structures.BeadID(cg_bead['name'], cg_bead['degree'])
                topo_at = copy.copy(
                    settings.cg_molecules[cg_bead['chain_name']]
                    .molecule_beads[cg_key])
            except KeyError:
                cg_key = structures.BeadID(cg_bead['name'], '*')
                topo_at = copy.copy(
                    settings.cg_molecules[cg_bead['chain_name']]
                    .molecule_beads[cg_key])
            if topo_at.atom_type not in atom_types:
                atom_types[topo_at.atom_type] = {
                    'name': topo_at.atom_type,
                    'mass': topo_at.mass,
                    'charge': 0.00,
                    'type': 'V',
                    'sigma': 1.0,
                    'epsilon': 1.0
                }
            cg_old_new_id[cg_id] = new_at_id
            topo_at.chain_idx = cg_bead['res_id']
            topo_at.cgnr = cg_bead['res_id']
            topo_at.atom_id = new_at_id
            output_topology.atoms[new_at_id] = topo_at
            if new_at_id not in chain2atoms:
                chain2atoms[new_at_id] = {}
            new_at_id += 1
    else:
        for cg_id in cg_ids:
            cg_bead = cg_graph.node[cg_id]
            cg_bead_id = new_at_id
            name = cg_bead['chain_name']
            hyb_file.atoms[new_at_id] = files_io.Atom(
                atom_id=cg_bead_id,
                name=cg_bead['name'],
                chain_name=cg_bead['chain_name'],
                chain_idx=cg_bead['res_id'],
                position=cg_bead['position'])
            try:
                cg_key = structures.BeadID(cg_bead['name'], cg_bead['degree'])
                topo_at = copy.copy(
                    settings.cg_molecules[cg_bead['chain_name']]
                    .molecule_beads[cg_key])
            except KeyError:
                cg_key = structures.BeadID(cg_bead['name'], '*')
                topo_at = copy.copy(
                    settings.cg_molecules[cg_bead['chain_name']]
                    .molecule_beads[cg_key])
            if topo_at.atom_type not in atom_types:
                atom_types[topo_at.atom_type] = {
                    'name': topo_at.atom_type,
                    'mass': topo_at.mass,
                    'charge': 0.00,
                    'type': 'V',
                    'sigma': 1.0,
                    'epsilon': 1.0
                }
            cg_old_new_id[cg_id] = new_at_id
            topo_at.chain_idx = cg_bead['res_id']
            topo_at.cgnr = cg_bead['res_id']
            topo_at.atom_id = new_at_id
            output_topology.atoms[new_at_id] = topo_at
            if new_at_id not in chain2atoms:
                chain2atoms[new_at_id] = {}
            new_at_id += 1

            # Insert atomistic part
            for aa_bead in cg_aa[name][cg_key]:
                hyb_file.atoms[new_at_id] = files_io.Atom(
                    atom_id=new_at_id,
                    name=aa_bead.name,
                    chain_name=aa_bead.chain_name,
                    chain_idx=cg_bead['res_id'],
                    position=(aa_bead.position +
                              (np.array(cg_bead['position']) -
                               np.array(cg_com[cg_bead['chain_name']][cg_key]))))
                topo_at = copy.copy(
                    input_aa_top[aa_bead.chain_name]
                    .chain_atom_names[aa_bead.chain_name][aa_bead.name][0])
                atom_id_old2new[aa_bead.chain_name][cg_bead['res_id']][topo_at.atom_id] = new_at_id
                topo_at.chain_idx = cg_bead['res_id']
                topo_at.cgnr = cg_bead['res_id']
                topo_at.atom_id = new_at_id
                output_topology.atoms[new_at_id] = topo_at
                atom_id2cg_bead_id[new_at_id] = cg_bead_id
                chain2atoms[cg_bead_id]['{}:{}'.format(
                    aa_bead.chain_name, aa_bead.name)] = new_at_id
                new_at_id += 1

    output_topology.atomtypes.update(atom_types)
    output_topology.atom_id_old2new = atom_id_old2new
    output_topology.atom_id2cg_bead_id = atom_id2cg_bead_id
    output_topology.cg_old_new_id = cg_old_new_id
    output_topology.chain2atoms = chain2atoms
    hyb_file.box = cg_graph.graph['box']
    return hyb_file, output_topology


def update_atomistic_bonds(output_topology, settings, plain=False):
    """Insert atomistic bonds from atomistic topology to hybrid topology with reindexing atom ids.

    Args:
        output_topology: The GROMACSTopology object.
        settings: The backmapping settings object.
        plain: Do not generate cross_terms.
    """
    print('Update atomistic bonded terms')
    input_aa_top = {}
    bonded_terms = ['bonds', 'angles', 'dihedrals', 'pairs']
    for mol_name, cg_mol in settings.cg_molecules.items():
        input_aa_top[mol_name] = files_io.GROMACSTopologyFile(cg_mol.source_topology)
        input_aa_top[mol_name].read()

    # Prepares output data structures in output_topology.
    for b in bonded_terms:
        output_topology.new_data[b] = {}
        if not plain:
            output_topology.new_data['cross_%s' % b] = {}

    # Replicate bonded lists by taking list from atomistic topology and
    # replicate those bonded terms which are involved in current atom set.
    for mol_name, chains in output_topology.atom_id_old2new.items():
        aa_top = input_aa_top[mol_name]
        for chain_idx in chains:
            old2new = chains[chain_idx]
            for bt in bonded_terms:
                bonded_list = getattr(aa_top, bt)
                for at_old, at_new in old2new.items():
                    at_cg_bead_id = output_topology.atom_id2cg_bead_id[at_new]
                    for b_def, b_pref in bonded_list.items():
                        if at_old in b_def:  # Found bonded term.
                            new_b_def = tuple(map(old2new.get, b_def))
                            # If one of the term in bonded tuple, triplet, quadruplet
                            # is not present then ignore that entry. This can happen
                            # for example when we use different fragments with respect to
                            # degree of CG beads.
                            if None in new_b_def:
                                continue
                            cg_bead_ids = list(map(output_topology.atom_id2cg_bead_id.get, new_b_def))
                            if plain == 'AA':
                                output_topology.new_data[bt][new_b_def] = b_pref
                            elif len(cg_bead_ids) != cg_bead_ids.count(at_cg_bead_id):
                                output_topology.new_data['cross_%s' % bt][new_b_def] = b_pref
                            else:
                                output_topology.new_data[bt][new_b_def] = b_pref


def generate_crosslink_at_terms(settings, output_topology, plain=False):
    """Generates cross-link atomistic terms in topology.

    Args:
        settings: The BackmapperSettings object.
        output_topology: The output topology.
        plain: If set to True then move everything to 'bonded' section from cross_.
    """
    print('Generate crosslinking terms at atomistic scale')
    b_prefix = '' if plain in ('AA', 'CG') else 'cross_'
    cg_ids = set(output_topology.cg_old_new_id.values())
    hybrid_bonds = settings.hybrid_topology['hybrid_bonds']
    new_bonds = {}
    for (b0, b1) in output_topology.new_data['{}bonds'.format(b_prefix)]:
        # Look only on CG bonds
        if b0 in cg_ids and b1 in cg_ids:
            at0 = output_topology.atoms[b0]
            at1 = output_topology.atoms[b1]
            if at0.active_site is not None and at1.active_site is not None:
                as0 = output_topology.atoms[
                    output_topology.chain2atoms[at0.atom_id][at0.active_site]]
                as1 = output_topology.atoms[
                    output_topology.chain2atoms[at1.atom_id][at1.active_site]]
                key0 = ('{}:{}'.format(at0.chain_name, at0.name),
                        '{}:{}'.format(at1.chain_name, at1.name))
                key1 = tuple(reversed(key0))
                hyb_bond_def = hybrid_bonds.get(key0, hybrid_bonds.get(key1))
                if hyb_bond_def is not None:
                    new_bonds[(as0.atom_id, as1.atom_id)] = hyb_bond_def
    output_topology.new_data['{}bonds'.format(b_prefix)].update(
        {k: v['bond_params'] for k, v in new_bonds.items()})
    # Search for angles and dihedrals if angle_params and dihedral_params present.
    g = output_topology.get_graph()
    for pair, hyb_params in new_bonds.items():
        gen_angle_params = hyb_params.get('angle_params')
        gen_dih_params = hyb_params.get('dihedral_params')
        gen_pair_params = hyb_params.get('pair_params')
        if gen_angle_params:
            output_topology.new_data['{}angles'.format(b_prefix)].update(
                {x: gen_angle_params for x in tools.gen_bonded_tuples(g, 3, pair)})
        dih_list = None
        if gen_dih_params or gen_pair_params:
            dih_list = tools.gen_bonded_tuples(g, 4, pair)
        if gen_dih_params:
            output_topology.new_data['{}dihedrals'.format(b_prefix)].update(
                {x: gen_dih_params for x in dih_list})
        if gen_pair_params:
            output_topology.new_data['{}pairs'.format(b_prefix)].update(
                {(x[0], x[3]): gen_pair_params for x in dih_list})


def postprocess_topology(settings, output_topology):
    print('Postprocess output topology file')
    if settings.hybrid_topology['include']:
        output_topology.header_section.append(settings.hybrid_topology['include'])
    output_topology.moleculetype = settings.hybrid_topology['moleculetype']
    output_topology.system_name = settings.hybrid_topology['system']
    output_topology.molecules = {
        'name': output_topology.system_name,
        'mol': 1
    }
