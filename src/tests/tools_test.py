#!/usr/bin/env python
"""
Copyright (C) 2017 Jakub Krajniak <jkrajniak@gmail.com>

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

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import tools
import files_io


class GetAtomisticTopologyTest(unittest.TestCase):
    def setUp(self):
        self.top = files_io.GROMACSTopologyFile('hyb_topol.top')
        self.top.read()

    def tearDown(self):
        self.clean_file('aa_top.top')

    def clean_file(self, filename):
        if os.path.exists(filename):
            os.remove(filename)

    def test_convert_cross(self):
        self.assertTrue(self.top.cross_bonds)
        aa_topol = tools.get_atomistic_topology(self.top)
        self.assertFalse(self.top.atomtypes)
        self.assertEqual({x.atom_type for x in self.top.atoms.values()}, {'CH2', 'CH3'})
        self.assertEqual(self.top.angletypes['CH3']['CH2']['CH2'],
                         {'params': ['2', '1'], 'func': 8})
        self.assertEqual(self.top.bondtypes['CH2']['CH3'],
                         {'params': ['1', '1'], 'func': 1})
        self.assertEqual(self.top.dihedraltypes['CH2']['CH3']['CH2']['CH3'],
                         {'params': ['1', '1'], 'func': 8})
        self.assertFalse(self.top.cross_bonds)
        aa_topol.write('aa_top.top')

        aa = files_io.GROMACSTopologyFile('aa_top.top')
        aa.read()

        # Check if cross section was removed
        for l in aa.content:
            if l.startswith(';') or l.startswith('#'):
                continue
            self.assertFalse('cross' in l)

if __name__ == '__main__':
    unittest.main()

