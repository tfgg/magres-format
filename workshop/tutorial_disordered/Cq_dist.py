import numpy
import sys, os
from matplotlib import pyplot as plt

from magres.utils import find_all
from magres.atoms import MagresAtoms

# Find all our calculation output files
files = find_all(sys.argv[1], '.magres')

# Build a list of all our loaded calculations
all_atoms = []

for file in files:
  all_atoms.append(MagresAtoms.load_magres(file))

# Structure for all our Cqs to go into, indexed by site and number of Si neighbours
Al_Cqs = {'T1': {0: [], 1: [], 2: [], 3: [], 4: []},
          'T2': {0: [], 1: []},}

print "Building histograms of Cq on Al and Si sites"

# Loop over all calculations
for atoms in all_atoms:
  print atoms.magres_file.path

  for Al_atom in atoms.species('Al'):
    # Find all Al and Si neighbours within 3.5 angstroms of this atom
    neighbours = atoms.species(['Al', 'Si']).within(Al_atom, 3.5)

    # Classify the site
    if len(neighbours) == 5:
      site = 'T1'
    else:
      site = 'T2'

    # Count the number of Si neighbours
    num_Si = len(neighbours.species('Si'))

    Al_Cqs[site][num_Si].append(abs(Al_atom.efg.Cq))

# Calculate the mean and standard deviations of Cq for each Al site and number of Si neighbours
for site in Al_Cqs:
  for num_si in Al_Cqs[site]:
    print "Al%s, n_si=%d, Cq=%.2f +- %.2f" % (site, num_si, numpy.mean(Al_Cqs[site][num_si]), numpy.std(Al_Cqs[site][num_si]))

sys.exit(1)

plt.hist(Al_Cqs['T1'][0], bins=20)
plt.hist(Al_Cqs['T1'][1], bins=20)
plt.hist(Al_Cqs['T1'][2], bins=20)
plt.hist(Al_Cqs['T1'][3], bins=20)
plt.hist(Al_Cqs['T1'][4], bins=20)
plt.show()


