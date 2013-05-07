import sys
from magres.atoms import MagresAtoms

# Load T1Si0.magres sample into an atoms structure
atoms = MagresAtoms.load_magres("samples/T1Si0.magres")

# Get all atoms within 8 angstrom of Al 15 and print out their isotope, name and EFG Cq
atomAl = atoms.get_species('Al', 15)
for image in atoms.within(atomAl, 8.0):
  print image.atom, atomAl.dist(image), image.atom.efg.Cq, image.atom.efg.evecs[2], image.atom.efg.evals

