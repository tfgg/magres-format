#!python

"""
  Only show couplings where the coupling in the opposite direction exists in another calculation.
"""

import os, os.path
import sys
import argparse

from magres.atoms import MagresAtoms
from magres.utils import load_all_magres, get_numeric

parser = argparse.ArgumentParser(description='Extract J-coupling parameters in both directions')
parser.add_argument('-J', '--J_tensor', action="store_const", help="Display J tensor", default=False, const=True)
parser.add_argument('-N', '--numbers', action="store_const", help="Parse numbers from path and print. This is useful for e.g. convergence calculations.", default=False, const=True)
parser.add_argument('source_dir', help='Directory to look for calculations in')
parser.add_argument('atom_species1', nargs='?', type=str, default=None, help='Only print couplings from this atomic species.')
parser.add_argument('atom_index1', nargs='?', type=int, default=None, help='Only print couplings from this atom.')
parser.add_argument('atom_species2', nargs='?', type=str, default=None, help='Only print couplings to this atomic species.')
parser.add_argument('atom_index2', nargs='?', type=int, default=None, help='Only print couplings to this atom.')

a = parser.parse_args(sys.argv[1:])

find_s1 = a.atom_species1
find_i1 = a.atom_index1
find_s2 = a.atom_species2
find_i2 = a.atom_index2

all_Js = {}

tensors = ['isc', 'isc_fc', 'isc_spin', 'isc_orbital_p', 'isc_orbital_d']

if a.J_tensor:
  print "# Showing in Hz (J)"
  property = "J_iso"
else:
  print "# Showing in 10e19.T^2.J^-1 (K)"
  property = "K_iso"

print "# Number\tAtom1\tAtom2\t{}\tDist\tPath".format("\t".join(tensors))

lines = []

if os.path.isfile(a.source_dir):
  magres_atoms = [MagresAtoms.load_magres(a.source_dir)]
else:
  magres_atoms = load_all_magres(a.source_dir)

for i, atoms in enumerate(magres_atoms):
  num = get_numeric(atoms.magres_file.path)

  if num:
    idx = num
  else:
    idx = [i]

  for isc in atoms.isc.perturbing(find_s1, find_i1).receiving(find_s2, find_i2):
    atom1 = isc.atom1
    atom2 = isc.atom2

    all_tensors = [getattr(atoms, tensor).perturbing(atom1).receiving(atom2)[0] for tensor in tensors]

    tensor_strs = ["{:.3f}".format(getattr(isc_, property)) for isc_ in all_tensors]

    lines.append((idx,
                  atoms.magres_file.path,
                  str(atom1),
                  str(atom2),
                  tensor_strs,
                  "\t%.3F" % atom1.dist(atom2)))

lines = sorted(lines, key=lambda xs: xs[0])

for idx, path, atom1, atom2, data, dist in lines:
  if a.numbers:
    print " ".join(map(str,idx)), atom1, atom2, "\t".join(data), dist, path
  else:
    print atom1, atom2, "\t".join(data), dist, path


