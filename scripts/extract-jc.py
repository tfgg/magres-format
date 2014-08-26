#!python

"""
  Only show couplings where the coupling in the opposite direction exists in another calculation.
"""

import os
import sys
import argparse

from magres.utils import load_all_magres, get_numeric

parser = argparse.ArgumentParser(description='Extract J-coupling parameters in both directions')
parser.add_argument('-J', '--J_tensor', action="store_const", help="Display J tensor", default=False, const=True)
parser.add_argument('source_dir', help='Directory to look for calculations in')
parser.add_argument('species', nargs=argparse.REMAINDER, help='Species to look for')

a = parser.parse_args(sys.argv[1:])

find_s1 = str(a.species[0])
find_i1 = int(a.species[1])

if len(a.species) >= 3:
  find_s2 = str(a.species[2])
else:
  find_s2 = None

if len(a.species) >= 4:
  find_i2 = int(a.species[3])
else:
  find_i2 = None

all_Js = {}

tensors = ['isc', 'isc_fc', 'isc_spin', 'isc_orbital_p', 'isc_orbital_d']


if a.J_tensor:
  print "# Showing in Hz (J)"
else:
  print "# Showing in 10e19.T^2.J^-1 (K)"

print "# Number\tAtom1\tAtom2\t{}\tDist\tPath".format("\t".join(tensors))

lines = []

magres_atoms = load_all_magres(a.source_dir)

for i, atoms in enumerate(magres_atoms):
  num = get_numeric(atoms.magres_file.path)

  if num:
    idx = num
  else:
    idx = [i]

  for atom1 in atoms: 
    if atom1.species == find_s1 and atom1.index == find_i1 and hasattr(atom1, 'isc'):
      for atom2 in atom1.isc:

        if (find_s2 is None or atom2.species == find_s2) and \
           (find_i2 is None or atom2.index == find_i2) and \
           (atom1 != atom2):

          if a.J_tensor:
            tensor_strs = ["{:.3f}".format(getattr(atom1, tensor)[atom2].J_iso) for tensor in tensors]
          else:
            tensor_strs = ["{:.3f}".format(getattr(atom1, tensor)[atom2].K_iso) for tensor in tensors]

          lines.append((idx,
                        atoms.magres_file.path,
                        str(atom1),
                        str(atom2),
                        tensor_strs,
                        "\t%.3F" % atom1.dist(atom2)))

lines = sorted(lines, key=lambda xs: xs[0])

for idx, path, atom1, atom2, data, dist in lines:
  print " ".join(map(str,idx)), atom1, atom2, "\t".join(data), dist, path


