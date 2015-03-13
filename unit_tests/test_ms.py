import math
import numpy
import unittest
import os
from magres.format import MagresFile
from magres.atoms import MagresAtoms

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data")

class MagresMsTest(unittest.TestCase):
  def test_nmr(self):
    """
      Check that the magnetic shielding tensors are being loaded.
    """
    atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "ethanol/ethanol-nmr.magres"))
    
    self.assertTrue(hasattr(atoms, "ms"))

    self.assertEqual(len(atoms.ms), 9)

  def test_replace_nmr(self):
    """
      Check that we can replace the tensor on an atom.
    """

    atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "ethanol/ethanol-nmr.magres"))

    with self.assertRaises(Exception):
      atoms.C1.ms.sigma = numpy.array([1,2,3])

    orig_iso = atoms.C1.ms.iso
    new_sigma = atoms.C1.ms.sigma * 2.0

    atoms.C1.ms.sigma = new_sigma
    
    self.assertEqual(atoms.C1.ms.iso, orig_iso * 2.0)

if __name__ == "__main__":
  unittest.main()
