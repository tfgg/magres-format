import numpy
import sys, os

from matplotlib import pyplot as plt
from numpy import tan, dot, arccos, sqrt
from magres.utils import find_all
from magres.atoms import MagresAtoms

# Find all our calculation output files
files = find_all(sys.argv[1], '.magres')

# Build a list of all our loaded calculations
all_atoms = []

for file in files:
  all_atoms.append(MagresAtoms.load_magres(file))

cqs = []
ghose = []

# Loop over all calculations
for atoms in all_atoms:
  print >>sys.stderr,atoms.magres_file.path

  for Al_atom in atoms.species('Al'):
    # Find all Al and Si neighbours within 3.5 angstroms of this atom
    neighbours = atoms.species(['Al', 'Si']).within(Al_atom, 3.5)

    # Classify the site
    if len(neighbours) == 5: # T1
      cqs.append(abs(Al_atom.efg.Cq))

      # Get the oxygen neighbours
      Os = atoms.species('O').within(Al_atom, 3.5)

      # Get the displacement vectors to the oxygens
      r1 = Os[0].position - Al_atom.position
      r2 = Os[1].position - Al_atom.position
      r3 = Os[2].position - Al_atom.position
      r4 = Os[3].position - Al_atom.position

      # Normalize
      r1 = r1 / sqrt(dot(r1,r1))
      r2 = r2 / sqrt(dot(r2,r2))
      r3 = r3 / sqrt(dot(r3,r3))
      r4 = r4 / sqrt(dot(r4,r4))

      tetra = arccos(-1.0/3)
      distort = tan(abs(arccos(dot(r1,r2)) - tetra))
      distort += tan(abs(arccos(dot(r1,r3)) - tetra))
      distort += tan(abs(arccos(dot(r1,r4)) - tetra))
      distort += tan(abs(arccos(dot(r2,r3)) - tetra))
      distort += tan(abs(arccos(dot(r2,r4)) - tetra))
      distort += tan(abs(arccos(dot(r3,r4)) - tetra))
      distort /= 6

      ghose.append(distort)

      print abs(Al_atom.efg.Cq), distort

plt.plot(ghose, cqs, 'x')
plt.show()

