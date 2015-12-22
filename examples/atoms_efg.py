import sys
import math
import numpy
from magres.atoms import MagresAtoms

# Load T1Si0.magres sample into an atoms structure
atoms = MagresAtoms.load_magres("../samples/T1Si0.magres")

atomAl15 = atoms.get('Al', 15)

for atom in atoms.species('Al'):
  if atom == atomAl15:
    continue

  Vzz = atom.efg.evecs[2]
  dr = atom.position - atomAl15.position

  Vzz_dr_ang = math.acos(numpy.dot(Vzz, dr) / math.sqrt(numpy.dot(dr,dr) * numpy.dot(Vzz,Vzz)))

  print(atom, atom.efg.Cq, Vzz_dr_ang)

