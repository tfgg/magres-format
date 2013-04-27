import sys
from magres.atoms import MagresAtoms

# Load T1Si0.magres sample into an atoms structure
atoms = MagresAtoms.load_magres("samples/T1Si0.magres")

# Get all atoms within 8 angstrom of Al 15 and print out their isotope, name and EFG Cq
atomAl = atoms.get_species('Al', 1)
for atom in atoms.within(atomAl, 2.0):
  print atom, atom.dist(atomAl), atomAl.dist(atom), atom.efg.Cq

