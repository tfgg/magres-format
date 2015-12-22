import math
import numpy
import unittest
import os
from magres.format import MagresFile
from magres.atoms import MagresAtoms

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data")


class AtomsTest(unittest.TestCase):
    cubic = MagresFile.load_json(
        '{"atoms": {"units": [["atom", "Angstrom"], ["lattice", "Angstrom"]], "lattice": [[[1.0, 0.0, 0.0], [0.0, '
        '1.0, 0.0], [0.0, 0.0, 1.0]]], "atom": [{"index": 1, "position": [0.0, 0.0, 0.0], "species": "H", '
        '"label": "H"}]}}')

    orthorhombic = MagresFile.load_json(
        '{"atoms": {"units": [["atom", "Angstrom"], ["lattice", "Angstrom"]], "lattice": [[[1.0, 0.0, 0.0], [0.0, '
        '2.0, 0.0], [0.0, 0.0, 3.0]]], "atom": [{"index": 1, "position": [0.0, 0.0, 0.0], "species": "H", '
        '"label": "H"}]}}')

    species = MagresFile.load_json(
        '{"atoms": {"units": [["atom", "Angstrom"], ["lattice", "Angstrom"]], "lattice": [[[1.0, 0.0, 0.0], [0.0, '
        '1.0, 0.0], [0.0, 0.0, 1.0]]], "atom": [{"index": 1, "position": [0.0, 0.0, 0.0], "species": "H", '
        '"label": "H"},{"index": 1, "position": [0.0, 0.0, 1.0], "species": "C", "label": "C1"},{"index": 1, '
        '"position": [0.0, 0.0, 2.0], "species": "C", "label": "C2"}]}}')

    def test_within_cubic(self):
        atoms = MagresAtoms.load_magres(self.cubic)

        self.assertEqual(len(atoms), 1)

        p = atoms[0].position

        self.assertEqual(len(atoms.within(p, 1.0)), 7)

        num_within_3 = 0
        for atom in atoms.within(atoms[0], 5.0):
            if atom.dist(atoms[0]) <= 3.0:
                num_within_3 += 1

        self.assertEqual(len(atoms.within(p, 3.0)), num_within_3)

        self.assertEqual(len(atoms.within(p, 5.0).within(p, 2.0)), len(atoms.within(p, 2.0).within(p, 5.0)))

    def test_within_orthorhombic(self):
        atoms = MagresAtoms.load_magres(self.orthorhombic)

        self.assertEqual(len(atoms), 1)

        p = atoms[0].position

        self.assertEqual(len(atoms.within(p, 1.0)), 3)

        num_within_3 = 0
        for atom in atoms.within(p, 10.0):
            if atom.dist(atoms[0]) <= 3.0:
                num_within_3 += 1

        self.assertEqual(len(atoms.within(p, 3.0)), num_within_3)
        self.assertEqual(len(atoms.within(p, 5.0).within(p, 2.0)), len(atoms.within(p, 2.0).within(p, 5.0)))

    def test_species(self):
        atoms = MagresAtoms.load_magres(self.species)

        self.assertEqual(len(atoms), 3)
        self.assertEqual(len(atoms.species('H')), 1)
        self.assertEqual(len(atoms.species('C')), 2)

    def test_labels(self):
        atoms = MagresAtoms.load_magres(self.species)

        # self.assertEquals(len(atoms), 3)
        # self.assertEquals(len(atoms.label('H')), 1)
        # self.assertEquals(len(atoms.label('C')), 0)
        # self.assertEquals(len(atoms.label('C1')), 1)
        # self.assertEquals(len(atoms.label('C2')), 1)


if __name__ == "__main__":
    unittest.main()
