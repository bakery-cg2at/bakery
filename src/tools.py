"""
Copyright (C) 2015-2017 Jakub Krajniak <jkrajniak@gmail.com>

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
import files_io
import networkx as nx
import sys

__doc__ = "Tool functions."


class MyArgParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(MyArgParser, self).__init__(*args, **kwargs)

    def convert_arg_line_to_args(self, line):
        for arg in line.split():
            t = arg.strip()
            if not t:
                continue
            if t.startswith('#'):
                break
            if not t.startswith('--'):
                t = '--{}'.format(t)
            yield t

    @staticmethod
    def save_to_file(output_file, namespace):
        """Saves arguments to file so it can be read again.

        Args:
            output_file: The string with the name of output file.
            namespace: The namespace with arguements.
        """
        with open(output_file, "w") as of:
            for k, v in namespace.__dict__.iteritems():
                if v is not None:
                    of.write('{}={}\n'.format(k, v))


def gen_bonded_tuples(g, num, bond_pair):
    """Generates tuples of different size, based on the graph and input edge.

    Args:
        g: The networkx Graph object.
        num: The length of the tuple.
        bond_pair: The edge which has to be included in all tuples.

    Returns:
        The set of all tuples of defined length from graph `g`.
    """
    b0, b1 = bond_pair
    paths = []
    if num > 3:
        for nb0 in g.edge[b0]:
            paths.extend(nx.single_source_shortest_path(g, nb0, num-1).values())
        for nb1 in g.edge[b1]:
            paths.extend(nx.single_source_shortest_path(g, nb1, num-1).values())

    paths.extend(nx.single_source_shortest_path(g, b0, num-1).values())
    paths.extend(nx.single_source_shortest_path(g, b1, num-1).values())
    output = set()
    for b in paths:
        if len(b) == num and b0 in b and b1 in b:
            if tuple(reversed(b)) not in output:
                output.add(tuple(b))
    return output


def get_graph(settings):
    """Build graph based on settings file. Useful for GROMACS."""
    gro = files_io.GROFile(settings.cg_configuration['file'])
    gro.read()

    g = nx.Graph(box=gro.box)
    for at_id, at in gro.atoms.iteritems():
        g.add_node(
            at_id,
            name=at.name,
            res_id=at.chain_idx,
            position=at.position,
            chain_name=at.chain_name)

    # Adding edges
    for mol in gro.chains:
        try:
            cg_bonds = settings.cg_molecules[mol].molecule_topology.get('bond')
        except KeyError:
            print(('\nError:\nMolecule \'{}\' not found in input CG trajectory'
                   '(valid molecule\' names: {})\nExit, nothing to do.'
                   ).format(mol, settings.cg_molecules.keys()))
            sys.exit(1)
        if cg_bonds:
            for chain_idx in gro.chains[mol]:
                for bond_name, bond_def in cg_bonds.iteritems():
                    for b1, b2 in bond_def['list']:
                        a1 = gro.chains[mol][chain_idx][b1]
                        a2 = gro.chains[mol][chain_idx][b2]
                        g.add_edge(a1.atom_id, a2.atom_id, params=bond_def['params'])
    # Update degree
    for n_id in g.node:
        g.node[n_id]['degree'] = g.degree(n_id)

    return g


def generate_list(input_list, id_map):
    output = {}
    for k, v in input_list.items():
        try:
            new_k = tuple(map(lambda x: id_map[x], k))
            output[new_k] = v
        except KeyError:
            continue
    return output


def get_atomistic_topology(in_top):
    """Returns atomistic topology from hybrid topology.

    Args:
        in_top: Input hybrid topology.
        coord: Input coordinate file.
    Returns:
        atomistic topology.
    """
    virtual_atomtypes = {k for k, v in in_top.atomtypes.items() if v['type'] == 'V'}

    # Map topol id -> atom_id
    topol_old2new = {}
    new_topol_atoms = {}
    new_id = 1
    for old_id in sorted(in_top.atoms):
        at = in_top.atoms[old_id]
        if at.atom_type in virtual_atomtypes:
            continue
        at.atom_id = new_id
        at.cgnr = at.chain_idx
        new_topol_atoms[new_id] = at
        topol_old2new[old_id] = new_id
        new_id += 1

    in_top.atomtypes = {k: v for k, v in in_top.atomtypes.items() if v['type'] != 'V'}

    bondtypes = {}
    for i in in_top.bondtypes:
        for j, params in in_top.bondtypes[i].items():
            if not {i, j}.intersection(virtual_atomtypes):
                if i not in bondtypes:
                    bondtypes[i] = {}
                bondtypes[i][j] = params
    in_top.bondtypes = bondtypes
    angletypes = {}
    for i in in_top.angletypes:
        for j in in_top.angletypes[i]:
            for k, params in in_top.angletypes[i][j].items():
                if not {i, j, k}.intersection(virtual_atomtypes):
                    if i not in angletypes:
                        angletypes[i] = {j: {k: params}}
                    elif j not in angletypes[i]:
                        angletypes[i][j] = {k: params}
                    else:
                        angletypes[i][j][k] = params
    in_top.angletypes = angletypes
    dihedraltypes = {}
    for i in in_top.dihedraltypes:
        for j in in_top.dihedraltypes[i]:
            for k in in_top.dihedraltypes[i][j]:
                for l in in_top.dihedraltypes[i][j][k]:
                    if not {i, j, k, l}.intersection(virtual_atomtypes):
                        if i not in dihedraltypes:
                            dihedraltypes[i] = {j: {k: {l: params}}}
                        elif j not in dihedraltypes[i]:
                            dihedraltypes[i][j] = {k: {l: params}}
                        elif k not in dihedraltypes[i][j]:
                            dihedraltypes[i][j][k] = {l: params}
                        else:
                            dihedraltypes[i][j][k][l] = params
    in_top.dihedraltypes = dihedraltypes

    in_top.header_section.insert(0, '; input_topol: {}\n'.format(in_top.file_name))
    in_top.header_section.insert(1, '; clean: {}\n; remove_cross: {}\n'.format(True, True))

    in_top.atoms = new_topol_atoms

    new_bonds = generate_list(in_top.bonds, topol_old2new)
    new_angles = generate_list(in_top.angles, topol_old2new)
    new_dihs = generate_list(in_top.dihedrals, topol_old2new)
    new_pairs = generate_list(in_top.pairs, topol_old2new)
    new_cr_bonds = generate_list(in_top.cross_bonds, topol_old2new)
    new_cr_angles = generate_list(in_top.cross_angles, topol_old2new)
    new_cr_dihs = generate_list(in_top.cross_dihedrals, topol_old2new)
    new_cr_pairs = generate_list(in_top.cross_pairs, topol_old2new)

    in_top.bonds = new_bonds
    in_top.angles = new_angles
    in_top.dihedrals = new_dihs
    in_top.pairs = new_pairs

    in_top.bonds.update(new_cr_bonds)
    in_top.angles.update(new_cr_angles)
    in_top.dihedrals.update(new_cr_dihs)
    in_top.pairs.update(new_cr_pairs)
    in_top.cross_bonds = {}
    in_top.cross_angles = {}
    in_top.cross_dihedrals = {}
    in_top.cross_pairs = {}
    in_top.skip_cross = True

    return in_top
