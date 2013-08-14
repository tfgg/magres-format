import unittest

from magres.format import MagresFile, BadVersion

import math
import numpy
import unittest
import os
from magres.format import MagresFile
from magres.atoms import MagresAtoms

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data")

class FormatTest(unittest.TestCase):
  def test_nmr(self):
    f = MagresFile(open(os.path.join(DATA_DIR, "ethanol", "ethanol-nmr.magres")))

    self.assertTrue("calculation" in f.data_dict)
    self.assertTrue("atoms" in f.data_dict)
    self.assertTrue("magres" in f.data_dict)
    
    self.assertTrue("lattice" in f.data_dict["atoms"])
    self.assertTrue("units" in f.data_dict["atoms"])
    self.assertTrue("atom" in f.data_dict["atoms"])
    
    self.assertTrue("units" in f.data_dict["magres"])
    self.assertTrue("efg" in f.data_dict["magres"])
    self.assertTrue("ms" in f.data_dict["magres"])

    num_atoms = len(f.data_dict['atoms']['atom'])

    self.assertTrue(len(f.data_dict['magres']['efg']), num_atoms)
    self.assertTrue(len(f.data_dict['magres']['ms']), num_atoms)

  def test_nmr_new(self):
    f = MagresFile(open(os.path.join(DATA_DIR, "ethanol", "ethanol-nmr-new.magres")))

    self.assertTrue("calculation" in f.data_dict)
    self.assertTrue("atoms" in f.data_dict)
    self.assertTrue("magres" in f.data_dict)
    
    self.assertTrue("lattice" in f.data_dict["atoms"])
    self.assertTrue("units" in f.data_dict["atoms"])
    self.assertTrue("atom" in f.data_dict["atoms"])
    
    self.assertTrue("units" in f.data_dict["magres"])
    self.assertTrue("efg" in f.data_dict["magres"])
    self.assertTrue("ms" in f.data_dict["magres"])

    num_atoms = len(f.data_dict['atoms']['atom'])

    self.assertTrue(len(f.data_dict['magres']['efg']), num_atoms)
    self.assertTrue(len(f.data_dict['magres']['ms']), num_atoms)

  def test_badversion(self):
    with self.assertRaises(BadVersion):
      f = MagresFile(open(os.path.join(DATA_DIR, "noversion.magres")))
    
    with self.assertRaises(BadVersion):
      f = MagresFile(open(os.path.join(DATA_DIR, "badversion.magres")))

  def test_json(self):
    f1 = MagresFile(open(os.path.join(DATA_DIR, "ethanol", "ethanol-nmr.magres")))

    json = f1.as_json()

    f2 = MagresFile.load_json(json)

    self.assertEqual(f1.data_dict, f2.data_dict)

if __name__ == "__main__":
  unittest.main()

