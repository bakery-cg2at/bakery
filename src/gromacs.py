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

from . import tools

__doc__ = "Helper functions for processing GROMACS input files."


def generate_cg_bonded_terms(settings, cg_graph, output_topology, plain=False):
    """Generates CG bonded term based on the topology on XML file.

    Args:
        settings: The XML settings file.
        cg_graph: The networkx.Graph object.
        output_topology: The GROMACS topology object.
        plain: If set to True then cg_terms will not be in cross_ sections.
    """
    b_prefix = '' if plain else 'cross_'
    bond_label = '{}bonds'.format(b_prefix)
    ang_label = '{}angles'.format(b_prefix)
    dih_label = '{}dihedrals'.format(b_prefix)
    pair_label = '{}pairs'.format(b_prefix)

    # Prepare ds
    ff = {}
    for mol in settings.cg_molecules:
        ff[mol] = {'bonds': {}, 'angles': {}, 'dihedrals': {}, 'pairs': {}}
        for k in ff[mol]:
            k0 = k[:-1]
            top_term = settings.cg_molecules[mol].molecule_topology.get(k0, {})
            for b, b_def in top_term.items():
                for l in b_def['list']:
                    ff[mol][k][tuple(l)] = b_def['params'].split()
                    ff[mol][k][tuple(l)].append(b)
    cg_edges = set()
    for b1, b2 in cg_graph.edges():
        n1 = cg_graph.node[b1]
        n2 = cg_graph.node[b2]
        if n1['chain_name'] == n2['chain_name']:
            key = (n1['name'], n2['name'])
            bond_params = ff[n1['chain_name']]['bonds'].get(
                key, ff[n1['chain_name']]['bonds'].get((n2['name'], n1['name'])))
            k_new = list(map(output_topology.cg_old_new_id.get, (b1, b2)))
            assert None not in k_new
            if (tuple(k_new) not in output_topology.new_data[bond_label] and
                    (k_new[1], k_new[0]) not in output_topology.new_data[bond_label]):
                cg_edges.add(tuple(k_new))
                output_topology.new_data[bond_label][tuple(k_new)] = bond_params

    # Generates triplets and quadruplets. 2-nd pass around new topology file.
    new_g = output_topology.get_graph()
    for b1, b2 in cg_edges:
        n1, n2 = output_topology.atoms[b1], output_topology.atoms[b2]
        chain_name = n1.chain_name
        if chain_name == n2.chain_name:
            key = (n1.name, n2.name)
            params = ff[chain_name]
            angle_params = params.get('angles')
            dihedral_params = params.get('dihedrals')
            pairs_params = params.get('pairs')
            # Generates angles.
            if angle_params:
                triplets = {
                    tuple([output_topology.atoms[x].name for x in z]): tuple(z)
                    for z in tools.gen_bonded_tuples(new_g, 3, (b1, b2))
                }
                for tr, tr_ids in triplets.items():
                    a_params = angle_params.get(tr, angle_params.get(tuple(reversed(tr))))
                    if (a_params and tr_ids not in output_topology.new_data[ang_label] and
                            tuple(reversed(tr_ids)) not in output_topology.new_data[ang_label]):
                        output_topology.new_data[ang_label][tr_ids] = a_params
            # Generates dihedrals.
            if dihedral_params or pairs_params:
                quadruplets = {
                    tuple([output_topology.atoms[x].name for x in z]): tuple(z)
                    for z in tools.gen_bonded_tuples(new_g, 4, (b1, b2))
                }
                for q, q_ids in quadruplets.items():
                    d_params = dihedral_params.get(q, dihedral_params.get(tuple(reversed(q))))
                    p_params = pairs_params.get((q[0], q[3]), pairs_params.get((q[3], q[0])))
                    rq_ids = tuple(reversed(q_ids))
                    if (d_params and
                        q_ids not in output_topology.new_data[dih_label] and
                            rq_ids not in output_topology.new_data[dih_label]):
                        output_topology.new_data[dih_label][q_ids] = d_params
                    # Generates pairs.
                    if (p_params and
                        (q_ids[0], q_ids[3]) not in output_topology.new_data[pair_label] and
                            (q_ids[3], q_ids[0]) not in output_topology.new_data[pair_label]):
                        output_topology.new_data[pair_label][(q_ids[0], q_ids[3])] = p_params
