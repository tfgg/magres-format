# This is an example script for doing a rough parse of a magres format file without using the proper library

import sys

f = open(sys.argv[1])

lattice = []
atoms = []
isc_tensors = []

for line in f:
  cols = line.strip().split()

  if len(cols) == 0:
    continue

  if cols[0] == "lattice":
    lattice = map(float, cols[1:])
  elif cols[0] == "atom":
    atoms.append((cols[1], cols[2], int(cols[3]), map(float, cols[4:7])))
  elif cols[0] == "isc":
    isc_tensors.append((cols[1], int(cols[2]), cols[3], int(cols[4]), map(float, cols[5:14])))

print lattice
print atoms
print isc_tensors

