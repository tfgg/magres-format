import math
from numpy import mean
import unittest
import os

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data")


class GettingStartedTest(unittest.TestCase):
    def test(self):
        from magres.atoms import MagresAtoms

        atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, 'ethanol-all.magres'))

        self.assertEqual(len(atoms), 9)

        for atom in atoms:
            "{} {}".format(atom, atom.ms.iso)

        atoms.species('H').set_reference(10.0)

        for atom in atoms.species('H'):
            "{} {}".format(atom, atom.ms.cs)

        for atom in atoms.within(atoms.C1, 2.0):
            "{} {}".format(atom, atom.ms.iso, atom.ms.aniso)

        self.assertEqual(len(atoms.species('H').ms.iso), 6)

        self.assertAlmostEqual(mean(atoms.C1.bonded.species('H').ms.iso), 29.9838078506)
        self.assertAlmostEqual(mean(atoms.C2.bonded.species('H').ms.iso), 28.7851386496)
        self.assertAlmostEqual(mean(atoms.O1.bonded.species('H').ms.iso), 31.9849757497)

        self.assertEqual(atoms.C1.ms.sigma.shape, (3, 3))
        self.assertEqual(len(atoms.C1.ms.evecs), 3)
        self.assertEqual(len(atoms.C1.ms.evals), 3)


if __name__ == "__main__":
    unittest.main()
