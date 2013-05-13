import unittest
from magres.format import MagresFile
from magres.atoms import MagresAtoms

class AtomsTest(unittest.TestCase):
  cubic = MagresFile.load_json('{"atoms": {"units": [["atom", "Angstrom"], ["lattice", "Angstrom"]], "lattice": [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]], "atom": [{"index": 1, "position": [0.0, 0.0, 0.0], "species": "H", "label": "H"}]}}')
  orthorhombic = MagresFile.load_json('{"atoms": {"units": [["atom", "Angstrom"], ["lattice", "Angstrom"]], "lattice": [[[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]], "atom": [{"index": 1, "position": [0.0, 0.0, 0.0], "species": "H", "label": "H"}]}}')

  def test_within_cubic(self):
    atoms = MagresAtoms.load_magres(self.cubic)

    self.assertEquals(len(atoms), 1)

    self.assertEquals(len(atoms.within(atoms[0], 1.0)), 7)
   
    num_within_3 = 0
    for atom in atoms.within(atoms[0], 5.0):
      if atom.dist(atoms[0]) <= 3.0:
        num_within_3 += 1

    self.assertEquals(len(atoms.within(atoms[0], 3.0)), num_within_3)

  def test_within_orthorhombic(self):
    atoms = MagresAtoms.load_magres(self.orthorhombic)

    self.assertEquals(len(atoms), 1)

    self.assertEquals(len(atoms.within(atoms[0], 1.0)), 3)
   
    num_within_3 = 0
    for atom in atoms.within(atoms[0], 10.0):
      if atom.dist(atoms[0]) <= 3.0:
        num_within_3 += 1

    self.assertEquals(len(atoms.within(atoms[0], 3.0)), num_within_3)

if __name__ == "__main__":
  unittest.main()
