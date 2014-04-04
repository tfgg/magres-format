import math
import numpy
import unittest
import os
from magres.format import MagresFile
from magres.atoms import MagresAtoms

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data")

class MagresTest(unittest.TestCase):
  def test_isc(self):
    atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "ethanol-isc.magres"))

    self.assertTrue(hasattr(atoms.C2, "isc"))
    self.assertTrue(hasattr(atoms.C2, "isc_orbital_p"))
    self.assertTrue(hasattr(atoms.C2, "isc_orbital_d"))
    self.assertTrue(hasattr(atoms.C2, "isc_spin"))
    self.assertTrue(hasattr(atoms.C2, "isc_fc"))

    self.assertEqual(len(atoms.C2.isc), len(atoms))

    perturb_atom = atoms.get('C', 2)

    for isc in atoms.C2.isc.values():
      self.assertEqual(isc.atom1, perturb_atom)
      self.assertEqual(isc.K_iso, numpy.trace(isc.K)/3.0)
    
    # Check that the principal components are ordered by the Haeberlen convention
    for isc in atoms.C2.isc.values():
      self.assertTrue(abs(isc.K_evals[2] - isc.K_iso) >= abs(isc.K_evals[0] - isc.K_iso))
      self.assertTrue(abs(isc.K_evals[0] - isc.K_iso) >= abs(isc.K_evals[1] - isc.K_iso))

  def test_full_isc(self):
    magres_files = [MagresFile(open(os.path.join(DATA_DIR, "ethanol", f))) for f in os.listdir(os.path.join(DATA_DIR, "ethanol")) if f.startswith('ethanol-jc')]

    self.assertEqual(len(magres_files), 9)

    atoms = MagresAtoms.load_magres(magres_files)

    for atom in atoms:
      self.assertTrue(hasattr(atom, "isc"))

    self.assertEqual(len(atoms.isc), len(atoms))

    # Check every atom has couplings to every other atom
    for atom in atoms:
      self.assertEqual(len(atom.isc), len(atoms))

  def test_nmr(self):
    atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "ethanol/ethanol-nmr.magres"))
    
    self.assertTrue(hasattr(atoms, "ms"))
    self.assertTrue(hasattr(atoms, "efg"))

    self.assertEqual(len(atoms.efg), 9)
    self.assertEqual(len(atoms.ms), 9)

    # Check that the EFGs are traceless
    for efg in atoms.efg:
      self.assertTrue(numpy.trace(efg.V)/3.0 < 1e-12)

    # Check that the principal components are ordered by the Haeberlen convention
    # |Vzz| >= |Vxx| >= |Vyy|
    for efg in atoms.efg:
      self.assertTrue(abs(efg.evals[2]) >= abs(efg.evals[0]))
      self.assertTrue(abs(efg.evals[0]) >= abs(efg.evals[1]))

class AtomsTest(unittest.TestCase):
  cubic = MagresFile.load_json('{"atoms": {"units": [["atom", "Angstrom"], ["lattice", "Angstrom"]], "lattice": [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]], "atom": [{"index": 1, "position": [0.0, 0.0, 0.0], "species": "H", "label": "H"}]}}')

  orthorhombic = MagresFile.load_json('{"atoms": {"units": [["atom", "Angstrom"], ["lattice", "Angstrom"]], "lattice": [[[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]], "atom": [{"index": 1, "position": [0.0, 0.0, 0.0], "species": "H", "label": "H"}]}}')
  
  species = MagresFile.load_json('{"atoms": {"units": [["atom", "Angstrom"], ["lattice", "Angstrom"]], "lattice": [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]], "atom": [{"index": 1, "position": [0.0, 0.0, 0.0], "species": "H", "label": "H"},{"index": 1, "position": [0.0, 0.0, 1.0], "species": "C", "label": "C1"},{"index": 1, "position": [0.0, 0.0, 2.0], "species": "C", "label": "C2"}]}}')

  def test_within_cubic(self):
    atoms = MagresAtoms.load_magres(self.cubic)

    self.assertEquals(len(atoms), 1)

    p = atoms[0].position

    self.assertEquals(len(atoms.within(p, 1.0)), 7)
   
    num_within_3 = 0
    for atom in atoms.within(atoms[0], 5.0):
      if atom.dist(atoms[0]) <= 3.0:
        num_within_3 += 1

    self.assertEquals(len(atoms.within(p, 3.0)), num_within_3)

    self.assertEquals(len(atoms.within(p, 5.0).within(p, 2.0)), len(atoms.within(p, 2.0).within(p, 5.0)))

  def test_within_orthorhombic(self):
    atoms = MagresAtoms.load_magres(self.orthorhombic)

    self.assertEquals(len(atoms), 1)
    
    p = atoms[0].position

    self.assertEquals(len(atoms.within(p, 1.0)), 3)
   
    num_within_3 = 0
    for atom in atoms.within(p, 10.0):
      if atom.dist(atoms[0]) <= 3.0:
        num_within_3 += 1

    self.assertEquals(len(atoms.within(p, 3.0)), num_within_3)
    self.assertEquals(len(atoms.within(p, 5.0).within(p, 2.0)), len(atoms.within(p, 2.0).within(p, 5.0)))

  def test_species(self):
    atoms = MagresAtoms.load_magres(self.species)

    self.assertEquals(len(atoms), 3)
    self.assertEquals(len(atoms.species('H')), 1)
    self.assertEquals(len(atoms.species('C')), 2)

  def test_labels(self):
    atoms = MagresAtoms.load_magres(self.species)

    self.assertEquals(len(atoms), 3)
    self.assertEquals(len(atoms.label('H')), 1)
    self.assertEquals(len(atoms.label('C')), 0)
    self.assertEquals(len(atoms.label('C1')), 1)
    self.assertEquals(len(atoms.label('C2')), 1)

if __name__ == "__main__":
  unittest.main()
