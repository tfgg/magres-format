import sys
from magres.atoms import MagresAtoms

# Load T1Si0.magres sample into an atoms structure
atoms = MagresAtoms.load_magres("samples/T1Si0.magres")

# Get all atoms within 8 angstrom of Al 15 and print out their isotope, name and EFG Cq
for atom in atoms.within(atoms.get_species('Al', 15), 8.0):
  print atom.efg_isotope, atom, atom.efg.Cq

