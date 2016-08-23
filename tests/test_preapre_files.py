import collections
import os
import networkx as nx
import re
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))
print sys.path[0]

import structures


def remove_file(file_name):
    if os.path.exists(file_name):
        os.unlink(file_name)


class TestCaseParseSettings(unittest.TestCase):
    def setUp(self):
        self.settings_file = 'settings.xml'
        self.bck = structures.BackmapperSettings2(self.settings_file)
        self.bck.prepare_hybrid()

    def tearDown(self):
        remove_file('hyb_conf.gro')
        remove_file('hyb_topol.top')

    def test_check_structure(self):
        hyb_topology = self.bck.hyb_topology
        g = nx.Graph()
        atom_names = collections.defaultdict(list)
        for at_id, at in hyb_topology.atoms.items():
            g.add_node(at_id, name=at.name, chain_name=at.chain_name)
            atom_names[at.name].append(at_id)

        g.add_edges_from(hyb_topology.new_data['bonds'])
        g.add_edges_from(hyb_topology.new_data['cross_bonds'])

        deg = g.degree()
        # Match string and the maximum degree
        atsym2deg = [('C..$', 3), ('C.$', 4), ('O.$', 2), ('H.*$', 1), ('N1.*$', 2), ('N2.$', 3)]
        atsym2deg = [(re.compile(x), y) for x, y in atsym2deg]
        for atsym, max_deg in atsym2deg:
            for at, at_ids in atom_names.items():
                if atsym.match(at):
                    at_degs = [deg[i] for i in at_ids]
                    self.assertEqual(
                        len(at_degs),
                        at_degs.count(max_deg),
                        'Degree set of {}: {}, max_deg: {} ({})'.format(at, at_degs, max_deg, atsym.pattern))

    def test_check_removed_atoms(self):
        pass


if __name__ == '__main__':
    unittest.main()