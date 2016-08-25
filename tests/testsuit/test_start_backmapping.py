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
import os
import networkx as nx
import re
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/')))
print sys.path[0]

import espressopp

import structures
import files_io
import tools_sim as tools
import gromacs_topology

def remove_file(file_name):
    if os.path.exists(file_name):
        os.unlink(file_name)


class TestCaseParseSettings(unittest.TestCase):
    def setUp(self):
        self.settings_file = 'settings.xml'
        self.bck = structures.BackmapperSettings2(self.settings_file)
        self.bck.prepare_hybrid()

        self.input_conf = gromacs_topology.read('hyb_conf.gro', 'hyb_topol.top')
        self.system = espressopp.System()
        self.system.rng = espressopp.esutil.RNG()
        self.hyb_topology = self.bck.hyb_topology


    def tearDown(self):
        remove_file('hyb_conf.gro')
        remove_file('hyb_topol.top')

    def test_check_particle_lists(self):
        part_prop, all_particles, adress_tuple = tools.genParticleList(
            self.input_conf, use_velocity=True, adress=True, use_charge=True)
        self.assertEqual(part_prop, ['id', 'type', 'pos', 'res_id', 'mass', 'q', 'adrat', 'lambda_adr', 'vp'])
        self.assertEqual(len(adress_tuple), 6)
        # 6 molecules x CG bead + 6 x 27 atoms - 2 bonds, every bond cause that three atoms are ignored
        self.assertEqual(len(all_particles), 6*28-2*3)

        # Check particle properties.
        self.assertEqual([x.id for x in all_particles if x.adrat == 0],
                         [x.atom_id for x in self.hyb_topology.atoms.values() if x.name == 'A1'])
        self.assertEqual([x.lambda_adr for x in all_particles if x.adrat == 0],
                         [0.0 for x in self.hyb_topology.atoms.values() if x.name == 'A1'])

    def test_potentials(self):
        pass


if __name__ == '__main__':
    unittest.main()