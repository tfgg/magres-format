import sys
from magres.atoms import MagresAtoms

# Load T1Si0.magres sample into an atoms structure
atoms = MagresAtoms.load_magres("samples/T1Si0.magres")

# Get all atoms within 8 angstrom of Al 15 and print out their isotope, name and EFG Cq
atomAl = atoms.get_species('Al', 5)

#atoms.species(['Si', 'Al']).within(atomAl, 8.0)


#print len(atoms.within(atomAl, 1.0).within(atomAl, 8.0))

for image in atoms.within(atomAl, 8.0).species(['Al', 'Si']):
  print image, atomAl.dist(image), image.efg.Cq, image.efg.evecs[2], image.efg.evals

print "inverted"

for image in atoms.species(['Al', 'Si']).within(atomAl, 8.0):
  print image, atomAl.dist(image), image.efg.Cq, image.efg.evecs[2], image.efg.evals

