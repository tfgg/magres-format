import math
import numpy
import unittest
import os

from numpy import allclose

from magres.format import MagresFile
from magres.atoms import MagresAtoms

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data")

class MagresAngles(unittest.TestCase):
  def test_angles(self):
    atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "angles/angles1.magres"))
    
    self.assertTrue(allclose(atoms.angle(atoms.C1, atoms.C2, atoms.C4, degrees=True), 45.0))
    self.assertTrue(allclose(atoms.angle(atoms.C1, atoms.C2, atoms.C6, degrees=True), 90.0))
    
    self.assertTrue(allclose(atoms.angle(atoms.C2, atoms.C1, atoms.C4), math.pi/2.0))
    self.assertTrue(allclose(atoms.angle(atoms.C2, atoms.C1, atoms.C6), math.pi/4.0))

  def test_angles_boundary(self):
    atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "angles/angles1.magres"))
    
    self.assertTrue(allclose(atoms.angle(atoms.C2, atoms.C1, atoms.C3, degrees=True), 180.0))
    self.assertTrue(allclose(atoms.angle(atoms.C6, atoms.C1, atoms.C9, degrees=True), 90.0))

class MagresDihedral(unittest.TestCase):
  def test_dihedral1(self):
    atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "angles/dihedral1.magres"))
   
    self.assertTrue(allclose(atoms.dihedral(atoms.C1, atoms.C2, atoms.C3, atoms.C4, degrees=True), 0.0))
   
  def test_dihedral2(self):
    atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "angles/dihedral2.magres"))
   
    self.assertTrue(allclose(atoms.dihedral(atoms.C1, atoms.C2, atoms.C3, atoms.C4, degrees=True), 180.0))

  def test_dihedral3(self):
    atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "angles/dihedral3.magres"))
   
    self.assertTrue(allclose(atoms.dihedral(atoms.C1, atoms.C2, atoms.C3, atoms.C4, degrees=True), 0.0))
  
  def test_dihedral4(self):
    atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "angles/dihedral4.magres"))
  
    self.assertTrue(allclose(atoms.dihedral(atoms.C1, atoms.C2, atoms.C3, atoms.C4, degrees=True), 45.0))

if __name__ == "__main__":
  unittest.main()

