#!/home/green/bin/python

"""
  Only show couplings where the coupling in the opposite direction exists in another calculation.
"""

import os
import sys

from magres.format import BadVersion
from magres.atoms import MagresAtoms
from magres.utils import load_all_magres

cwd = "."

if len(sys.argv)>1:
    cwd = str(sys.argv[1])

magres_atoms = load_all_magres(cwd)

find_s1 = str(sys.argv[2])
find_i1 = int(sys.argv[3])

if len(sys.argv) >= 5:
  find_s2 = str(sys.argv[4])
else:
  find_s2 = None

if len(sys.argv) >= 6:
  find_i2 = int(sys.argv[5])
else:
  find_i2 = None

all_Js = {}

tensors = ['isc', 'isc_fc', 'isc_spin', 'isc_orbital_p', 'isc_orbital_d']

print "# Number\tPath\tAtom1\tAtom2\t" + "\t".join(tensors)

lines = []

def get_numeric(s):
  return "".join([c for c in s if ord("0") <= ord(c) <= ord("9")])

for atoms in magres_atoms:
  num = get_numeric(atoms.magres_file.path)

  if num != '':
    idx = int(num)
  else:
    idx = 0

  for atom1 in atoms: 
    if atom1.species == find_s1 and atom1.index == find_i1 and hasattr(atom1, 'isc'):
      for atom2 in atom1.isc:
        if (find_s2 is None or atom2.species == find_s2) and (find_i2 is None or atom2.index == find_i2) and (atom1 != atom2):
          lines.append((idx, atoms.magres_file.path, str(atom1) + "\t" + str(atom2) + "\t" + "\t".join(["%.2f" % getattr(atom1, tensor)[atom2].K_iso for tensor in tensors]) + "\t%.2F" % atom1.dist(atom2)))

lines = sorted(lines, key=lambda (x,y,z): x)

for idx, path, line in lines:
  print idx, path, line


