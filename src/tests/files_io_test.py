"""
Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>

This file is part of Backmapper.

Backmapper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import os
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import files_io

class Storage:
    def getParticle(self, particle_id):
        class P:
            def __init__(self, pid):
                self.pid = pid

            @property
            def pos(self):
                return (self.pid, self.pid, self.pid)

        return P(particle_id)

class System:
    storage = Storage()


class GROFileTestCase(unittest.TestCase):

    def test_create_file(self):
        gro_file = files_io.GROFile('abc.gro')
        gro_file.atoms[1] = files_io.Atom(
            atom_id=1,
            name='A1',
            chain_name='MOL',
            chain_idx=1,
            position=(1,1,1)
        )
        gro_file.title='XXX'
        gro_file.box = (10, 10, 10)
        gro_file.write('abc.gro', force=True)

        self.assertTrue(os.path.exists('abc.gro'))
        os.unlink('abc.gro')

    def test_copy_file(self):
        gro_file = files_io.GROFile('abc.gro')
        gro_file.atoms[1] = files_io.Atom(
            atom_id=1,
            name='A1',
            chain_name='MOL',
            chain_idx=1,
            position=(1,1,1)
        )
        gro_file.atoms[2] = files_io.Atom(
            atom_id=2,
            name='A2',
            chain_name='MOL',
            chain_idx=1,
            position=(2,2,2)
        )
        gro_file.title='XXX'
        gro_file.box = (10, 10, 10)

        new_gro = files_io.GROFile.copy(gro_file)
        self.assertEquals(new_gro.atoms, gro_file.atoms)
        self.assertEqual(new_gro.box, gro_file.box)

    def test_copy_part(self):
        gro_file = files_io.GROFile('abc.gro')
        gro_file.atoms[1] = files_io.Atom(
            atom_id=1,
            name='A1',
            chain_name='MOL',
            chain_idx=1,
            position=(1,1,1)
        )
        gro_file.atoms[2] = files_io.Atom(
            atom_id=2,
            name='A2',
            chain_name='MOL',
            chain_idx=1,
            position=(2,2,2)
        )
        gro_file.atoms[3] = files_io.Atom(
            atom_id=3,
            name='A3',
            chain_name='MOL',
            chain_idx=1,
            position=(2,2,2)
        )
        gro_file.title='XXX'
        gro_file.box = (10, 10, 10)

        new_gro = files_io.GROFile.copy(gro_file, particle_ids=[2], renumber=True)
        self.assertNotEquals(new_gro.atoms, gro_file.atoms)
        self.assertEquals(len(new_gro.atoms), 1)
        self.assertEquals(new_gro.atoms[1].name, 'A2')

        # Copy without renumber
        new_gro = files_io.GROFile.copy(gro_file, particle_ids=[2], renumber=False)
        self.assertNotEquals(new_gro.atoms, gro_file.atoms)
        self.assertEquals(len(new_gro.atoms), 1)
        self.assertEquals(new_gro.atoms[2].name, 'A2')

    def test_update_position(self):
        gro_file = files_io.GROFile('abc.gro')
        gro_file.atoms[1] = files_io.Atom(
            atom_id=1,
            name='A1',
            chain_name='MOL',
            chain_idx=1,
            position=(1,1,1)
        )
        gro_file.atoms[2] = files_io.Atom(
            atom_id=2,
            name='A2',
            chain_name='MOL',
            chain_idx=1,
            position=(2,2,2)
        )
        gro_file.atoms[3] = files_io.Atom(
            atom_id=3,
            name='A3',
            chain_name='MOL',
            chain_idx=1,
            position=(2,2,2)
        )
        gro_file.title='XXX'
        gro_file.box = (10, 10, 10)

        system = System()
        self.assertEquals(gro_file.atoms[3].position, (2, 2, 2))
        gro_file.update_positions(system)
        self.assertEquals(gro_file.atoms[3].position, (3, 3, 3))