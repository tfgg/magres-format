import math
import numpy
import unittest
import os
from magres.format import MagresFile
from magres.atoms import MagresAtoms

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data")


class MagresEfgTest(unittest.TestCase):
    def test_efg(self):
        atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "ethanol/ethanol-nmr.magres"))

        self.assertTrue(hasattr(atoms, "efg"))

        self.assertEqual(len(atoms.efg), 9)

        # Check that the EFGs are traceless
        for efg in atoms.efg:
            self.assertTrue(numpy.trace(efg.V) / 3.0 < 1e-12)

        # Check that the principal components are ordered by the Haeberlen convention
        # |Vzz| >= |Vxx| >= |Vyy|
        for efg in atoms.efg:
            self.assertTrue(abs(efg.evals[2]) >= abs(efg.evals[0]))
            self.assertTrue(abs(efg.evals[0]) >= abs(efg.evals[1]))

        # Check that eta is positive and <=1 (Cq can be positive or negative)
        # eta = (Vxx - Vyy)/Vzz, 0 <= eta <= 1
        for efg in atoms.efg:
            self.assertTrue((efg.evals[1] - efg.evals[0]) / efg.evals[2] >= 0)
            self.assertTrue((efg.evals[1] - efg.evals[0]) / efg.evals[2] <= 1)
            self.assertTrue(efg.eta >= 0)
            self.assertTrue(efg.eta <= 1)

    def test_replace_efg(self):
        """
          Check that we can replace the tensor on an atom.
        """

        atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "ethanol/ethanol-nmr.magres"))

        with self.assertRaises(Exception):
            atoms.C1.efg.V = numpy.array([1, 2, 3])

        orig_Cq = atoms.C1.efg.Cq
        new_V = atoms.C1.efg.V * 2.0

        atoms.C1.efg.V = new_V

        self.assertEqual(atoms.C1.efg.Cq, orig_Cq * 2.0)


if __name__ == "__main__":
    unittest.main()
