import math
import numpy
import unittest
import os
from magres.format import MagresFile
from magres.atoms import MagresAtoms

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data")

class MagresEfgTest(unittest.TestCase):
  def test_nmr(self):
    atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "ethanol/ethanol-nmr.magres"))
    
    self.assertTrue(hasattr(atoms, "efg"))

    self.assertEqual(len(atoms.efg), 9)

    # Check that the EFGs are traceless
    for efg in atoms.efg:
      self.assertTrue(numpy.trace(efg.V)/3.0 < 1e-12)

    # Check that the principal components are ordered by the Haeberlen convention
    # |Vzz| >= |Vxx| >= |Vyy|
    for efg in atoms.efg:
      self.assertTrue(abs(efg.evals[2]) >= abs(efg.evals[0]))
      self.assertTrue(abs(efg.evals[0]) >= abs(efg.evals[1]))

if __name__ == "__main__":
  unittest.main()
