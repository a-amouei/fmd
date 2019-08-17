#   test_molframe.py: This file is part of Free Molecular Dynamics
#
#   Copyright (C) 2019 Hossein Ghorbanfekr
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Testing molecular frame"""

import os
import unittest

import numpy as np
from ase.build import bulk
from ase.build.supercells import make_supercell

from pyfmd import Atom, Atoms, Box, MolecularFrame, SelectRegion


def atom1():
    return Atom('Cu', (0.3, 0.4, 0.7), 1, (1.0, 2.0, -3.0))


def atom2():
    return Atom(position=(0.5, 0.6, 0.8), velocity=(-1, -2, -3), symbol='Fe', group_id=2)


class AtomTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_atom_initialization(self):

        self.assertIsInstance(atom1(), Atom)
        self.assertEqual(atom1().symbol, 'Cu')
        self.assertEqual(atom1().x, 0.3)
        self.assertEqual(atom1().y, 0.4)
        self.assertEqual(atom1().z, 0.7)
        self.assertEqual(atom1().vx, 1.0)
        self.assertEqual(atom1().vy, 2.0)
        self.assertEqual(atom1().vz, -3.0)
        self.assertEqual(atom1().group_id, 1)
        self.assertEqual(atom1().mass, 63.546)

        with self.assertRaises(AssertionError):
            Atom(position=('a', 2, 2))

        with self.assertRaises(AssertionError):
            Atom(velocity=('a', 2, 2))

        with self.assertRaises(AssertionError):
            atom = Atom().mass

    def test_atom_position(self):

        atm1 = atom1().set_position((1, 2, 3))
        self.assertEqual(atm1.x, 1)
        self.assertEqual(atm1.y, 2)
        self.assertEqual(atm1.z, 3)

        with self.assertRaises(ValueError):
            atom1().translate((1, 2, 'a'))

    def test_atom_velocity(self):

        atm1 = atom1().set_velocity([1, 2, 3])
        self.assertEqual(atm1.vx, 1)
        self.assertEqual(atm1.vy, 2)
        self.assertEqual(atm1.vz, 3)

        with self.assertRaises(ValueError):
            atom1().set_velocity((1, 2, 'a'))

        atm1 = atom1().add_velocity([1, 2, 3])
        self.assertEqual(atm1.vx, 2.0)
        self.assertEqual(atm1.vy, 4.0)
        self.assertEqual(atm1.vz, 0.0)

        with self.assertRaises(ValueError):
            atom1().add_velocity((2, 'a'))

    def test_atom_translate(self):

        atm1 = atom1().translate((1, 2, 3))
        self.assertEqual(atm1.x, atom1().x + 1)
        self.assertEqual(atm1.y, atom1().y + 2)
        self.assertEqual(atm1.z, atom1().z + 3)

        with self.assertRaises(ValueError):
            atom1().translate((1, 2))


class AtomsTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_atoms_initialization(self):

        self.assertIsInstance(Atoms(), Atoms)
        self.assertEqual(Atoms().name, 'Atoms')

        atoms_section = Atoms(atom1())
        self.assertIsInstance(atoms_section, Atoms)
        self.assertEqual(atoms_section.get_number_of_atoms(), 1)
        self.assertEqual(atoms_section.atoms[0].symbol, atom1().symbol)

        with self.assertRaises(AssertionError):
            Atoms(1)

        with self.assertRaises(AssertionError):
            Atoms([1, 2])

    def test_atoms_add_atom(self):
        atoms_section = Atoms()
        atoms_section.append(atom1())
        self.assertIsInstance(atoms_section, Atoms)
        self.assertEqual(atoms_section.get_number_of_atoms(), 1)
        self.assertEqual(atoms_section.atoms[0].symbol, atom1().symbol)

    def test_atoms_section_list(self):
        atoms_section = Atoms([atom1(), atom2()])
        self.assertIsInstance(atoms_section, Atoms)
        self.assertEqual(atoms_section.get_number_of_atoms(), 2)
        self.assertEqual(atoms_section.atoms[0].symbol, atom1().symbol)
        self.assertEqual(atoms_section.atoms[1].symbol, atom2().symbol)

    def test_atoms_section_equal(self):
        atoms_section = Atoms(atom1())
        new_atoms_section = atoms_section
        self.assertIsInstance(new_atoms_section, Atoms)
        self.assertEqual(new_atoms_section.get_number_of_atoms(), atoms_section.get_number_of_atoms())
        self.assertEqual(new_atoms_section.atoms[0].symbol, atoms_section.atoms[0].symbol)

    def test_atoms_section_plus(self):
        atoms_section1 = Atoms(atom1())
        atoms_section2 = Atoms(atom2())
        new_atoms_section = atoms_section1 + atoms_section2
        self.assertIsInstance(new_atoms_section, Atoms)
        self.assertEqual(new_atoms_section.get_number_of_atoms(),
                         atoms_section1.get_number_of_atoms()+atoms_section2.get_number_of_atoms())
        self.assertEqual(atom1().symbol, new_atoms_section.atoms[0].symbol)
        self.assertEqual(atom2().symbol, new_atoms_section.atoms[1].symbol)

    def test_atoms_getitem(self):
        atoms = Atoms([atom1(), atom2()])
        for d in range(3):
            self.assertEqual(atom1().position[d], atoms[0].position[d])
            self.assertEqual(atom2().position[d], atoms[1].position[d])

        with self.assertRaises(ValueError):
            atoms[12]

    def test_atoms_setitem(self):
        atoms = Atoms([atom1(), atom2()])
        atom = Atom('Li', (15, 0, 0))
        atoms[1] = atom
        for d in range(3):
            self.assertEqual(atom.position[d], atoms[1].position[d])

        with self.assertRaises(AssertionError):
            atoms[0] = 1

        with self.assertRaises(AssertionError):
            atoms[12] = atom


def make_supercell_cu():
    crys = bulk('Cu', 'fcc', a=3.6, orthorhombic=True)
    P = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) * 5
    return make_supercell(crys, P)


def make_supercell_li():
    crys = bulk('Li', 'bcc', a=3.51, orthorhombic=True)
    P = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) * 10
    return make_supercell(crys, P)


class MolecularFrameTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_molecular_frame(self):
        mf = MolecularFrame()
        self.assertIsInstance(mf, MolecularFrame)
        self.assertEqual(mf.name, 'Molecular Frame')

    def test_molecular_frame_import_ase(self):
        mf = MolecularFrame('Cu crystal')
        mf.import_from(package_instance=make_supercell_cu(), package_name="ASE")
        self.assertEqual(mf.get_number_of_atoms(), make_supercell_cu().get_number_of_atoms())
        self.assertEqual(mf.name, 'Cu crystal')
        self.assertAlmostEqual(mf.atoms_section.atoms[1].x, make_supercell_cu().get_positions()[1][0])  # 2nd atom
        self.assertEqual(mf.get_atoms_section().atoms[1].symbol, 'Cu')  # second atom symbol

    def test_molecular_frame_equal(self):
        mf = MolecularFrame('Cu crystal').import_from(make_supercell_cu())
        new_mf = mf  # equal operator
        self.assertEqual(new_mf.get_number_of_atoms(), mf.get_number_of_atoms())
        self.assertEqual(new_mf.name, mf.name)
        self.assertEqual(new_mf.atoms_section.atoms[1].x, mf.atoms_section.atoms[1].x)  # atom x position
        self.assertEqual(new_mf.atoms_section.atoms[1].symbol, mf.atoms_section.atoms[1].symbol) # second atom symbol
        # TODO: limited to atoms section!

    def test_molecular_frame_plus(self):
        mf1 = MolecularFrame('Cu crystal').import_from(make_supercell_cu())
        mf2 = MolecularFrame('Li crystal').import_from(make_supercell_li())
        mf = mf1 + mf2 # plus operator
        self.assertIsInstance(mf, MolecularFrame)
        self.assertTrue(mf1.name in mf.name and mf2.name in mf.name)
        self.assertEqual(mf.get_number_of_atoms(), mf1.get_number_of_atoms() + mf2.get_number_of_atoms())
        # TODO: limited to atoms section!

    def test_molecular_frame_select_region(self):

        def sphere(x, y, z):
            """sphere region"""
            if np.sqrt(x**2+y**2+z**2) < 7.0:
                return True
            else:
                return False

        mf = MolecularFrame('Li crystal').import_from(make_supercell_li())
        sel_mf = mf.select(SelectRegion(sphere))
        self.assertIsInstance(sel_mf, MolecularFrame)
        self.assertEqual(sel_mf.get_number_of_atoms(), 12)
        # TODO: limited to atom sections

    def test_molecular_frame_clone(self):
        mf = MolecularFrame('Li crystal').import_from(make_supercell_li())
        mf_clone = mf.clone()
        mf.atoms_section.atoms[0].x = 1.0  # make a change to the first atom in mf
        self.assertIsInstance(mf_clone, MolecularFrame)
        self.assertTrue(not mf_clone.atoms_section.atoms[0].x == mf.atoms_section.atoms[0].x)
        self.assertEqual(mf_clone.atoms_section.atoms[0].y, mf.atoms_section.atoms[0].y)

    def test_molecular_frame_formatter(self):
        mf1 = MolecularFrame('Cu crystal').import_from(make_supercell_cu())
        filename = "cu_pytest.xyz"
        mf1.write(filename)
        mf2 = MolecularFrame().read(filename)
        self.assertEqual(mf1.get_number_of_atoms(), mf2.get_number_of_atoms())
        os.system("rm -f %s"%filename)

    def test_molecular_frame_adaptor(self):
        mf1 = MolecularFrame('Cu crystal').import_from(make_supercell_cu())
        mf2 = MolecularFrame().import_from(mf1.export())
        self.assertEqual(mf1.get_number_of_atoms(), mf2.get_number_of_atoms())
        self.assertTrue(Atoms().name in mf2.get_sections_name())
        self.assertTrue(Box().name in mf2.get_sections_name())

    def test_temperature(self):
        mf = MolecularFrame('Li crystal').import_from(make_supercell_li())
        mf.temperature = 50.0
        self.assertAlmostEqual(mf.temperature, 50.0)

        with self.assertRaises(AssertionError):
            mf.temperature = 0  # T > 0

    def test_center_of_mass(self):
        mf = MolecularFrame('Li crystal').import_from(make_supercell_li())
        for i in range(3):
            self.assertAlmostEqual(mf.center_of_mass[i], 16.672499, 5)
        # Set COM to center
        mf.center_of_mass = (0, 0, 0)
        for i in range(3):
            self.assertAlmostEqual(mf.center_of_mass[i], 0)


if __name__ == '__main__':
    unittest.main()
