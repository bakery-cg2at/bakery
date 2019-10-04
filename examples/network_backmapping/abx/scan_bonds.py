from __future__ import print_function

import networkx as nx
import pickle
import random
from collections import defaultdict
import logging

import sys

sys.setrecursionlimit(5000)


with open('graph_before_cross_bonds.pck', 'rb') as inputf:
    gb = pickle.load(inputf)


cg_cross_bonds = []
res_graph = nx.Graph()
cg_graph = nx.Graph()
for n in gb.nodes:
    n1 = gb.node[n]
    if n1['bead_type'] == 'CG':
        res_graph.add_node(n1['res_id'], chain_name=n1['chain_name'])
        cg_graph.add_node(n, chain_name=n1['chain_name'], res_id=n1['res_id'])


for e1, e2 in gb.edges:
    n1 = gb.node[e1]
    n2 = gb.node[e2]
    if n1['bead_type'] == 'CG' and n2['bead_type'] == 'CG':
        res_graph.add_node(n1['res_id'], chain_name=n1['chain_name'])
        res_graph.add_node(n2['res_id'], chain_name=n2['chain_name'])
        cg_graph.add_edge(e1, e2)
        if n1['res_id'] != n2['res_id']:
            cg_cross_bonds.append((e1, e2))
            res_graph.add_edge(n1['res_id'], n2['res_id'])


def get_active_sites(g_degree, global_graph, atom_ids):
    max_degree = {'MEA': 2, 'MEB': 3, 'CMB': 2}
    ret = []
    for a in atom_ids:
        n = global_graph.node[a]
        deg_n = g_degree[a]
        if n['name'] in max_degree and deg_n < max_degree[n['name']]:
            ret.append(a)
    return ret


def check_at_cross_bond(g_degree, global_graph, a1, a2):
    max_degree = {'MEA': 2, 'MEB': 3, 'CMB': 2}
    restricted = [('MEA', 'MEA'), ('MEB', 'MEB'), ('CMB', 'CMB'), ('MEB', 'CMB')]
    n1 = global_graph.node[a1]
    n2 = global_graph.node[a2]
    if n1['name'] not in max_degree or n2['name'] not in max_degree:
        logging.debug('{}({})-{}({}) not in max_degree'.format(n1['name'], a1, n2['name'], a2))
        return False
    if (n1['name'], n2['name']) in restricted or (n2['name'], n1['name']) in restricted:
        logging.debug('{}({})-{}({}) is restricted'.format(n1['name'], a1, n2['name'], a2))
        return False
    deg_1 = g_degree[a1]
    deg_2 = g_degree[a2]
    if deg_1 >= max_degree[n1['name']] or deg_2 >= max_degree[n2['name']]:
        logging.debug('{}({}) (deg:{})-{}({}) (deg:{}) >= max_degree'.format(n1['name'], a1, deg_1, n2['name'], a2, deg_2))
        return False
    return True


gb_bak = gb.copy()
global_degree = dict(gb_bak.degree())
global_graph = gb_bak.copy()

end_results = []

def process_graph(local_degree, at_bonds, cross_bonds, cg_index):
    logging.info('process index {}'.format(cg_index))
    if cg_index == len(cross_bonds):
        end_results.append(at_bonds)
        return
    cg_b1, cg_b2 = cross_bonds[cg_index]
    n1, n2 = global_graph.node[cg_b1], global_graph.node[cg_b2]
    possible_as1 = get_active_sites(local_degree, global_graph, n1['atom_ids'])
    possible_as2 = get_active_sites(local_degree, global_graph, n2['atom_ids'])
    cg_index_level = cg_index + 1
    for a1 in possible_as1:
        for a2 in possible_as2:
            is_correct_bond = check_at_cross_bond(local_degree, global_graph, a1, a2)
            logging.debug('{}-{} is_correct_bond={} {}'.format(a1, a2, is_correct_bond, cg_index_level))
            if is_correct_bond:
                new_local_degree = local_degree.copy()
                new_local_degree[a1] += 1
                new_local_degree[a2] += 1
                new_at_bonds = at_bonds[:]
                new_at_bonds.append((a1, a2))
                process_graph(new_local_degree, new_at_bonds, cross_bonds, cg_index_level)


logger = logging.getLogger()
logger.setLevel(logging.INFO)

process_graph(global_degree, [], cg_cross_bonds, 0)

print('Number of solutions: {}'.format(len(end_results)))

end_results.sort(key=lambda l: len(l), reverse=True)
print('Max number of at bonds: {}'.format(len(end_results[0])))

max_end_results = end_results[0]

with open('end_results.dat', 'wb+') as out_data:
    pickle.dump(end_results, out_data)

with open('predefined_active_sites.txt', 'w+') as predef_as:
    for e1, e2 in max_end_results:
        n1 = global_graph.node[e1]
        n2 = global_graph.node[e2]
        print('{} {} {} {}'.format(n1['orig_cg_bead_id'], n2['orig_cg_bead_id'], n1['name'], n2['name']), file=predef_as)
