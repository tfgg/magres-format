import sys
from magres.atoms import MagresAtoms

atoms = MagresAtoms.load_magres("samples/T1Si0.magres")

for atom in atoms.within(atoms.get_species('Al', 15), 8.0):
  print atom.efg_isotope, atom, atom.efg.Cq

