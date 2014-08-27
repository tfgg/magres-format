#!python

import os
import sys
import argparse

from magres.utils import load_all_magres, get_numeric

parser = argparse.ArgumentParser(description='Extract magnetic shielding parameters from a set of calculations.')
parser.add_argument('-N', '--numbers', action="store_const", help="Parse numbers from path and print. This is useful for e.g. convergence calculations.", default=False, const=True)
parser.add_argument('source_dir', help='Directory to look for calculations below.')
parser.add_argument('species', nargs=argparse.REMAINDER, help='Only print this species.')

a = parser.parse_args(sys.argv[1:])

magres_atoms = load_all_magres(a.source_dir)

find_s = str(a.species[0])

if len(a.species) >= 2:
  find_i = int(a.species[1])
else:
  find_i = None

tensors = ['ms']

lines = []

print "# Number\tAtom\tIso\tAniso\tAsym\tPath"

for i, atoms in enumerate(magres_atoms):
  num = get_numeric(atoms.magres_file.path)

  if num:
    idx = num
  else:
    idx = [i]

  for atom in atoms: 
    if atom.species == find_s and \
       (find_i is None or atom.index == find_i) and \
       hasattr(atom, 'efg'):

      lines.append((idx,
                    atoms.magres_file.path,
                    str(atom),
                    [atom.ms.iso, atom.ms.aniso, atom.ms.eta]))

lines = sorted(lines, key=lambda xs: xs[0])

for idx, path, atom, data in lines:
  if a.numbers:
    print " ".join(map(str, idx)) + "\t" + atom + "\t" + "\t".join("{:.3f}".format(x) for x in data) + "\t" + path
  else:
    print atom + "\t" + "\t".join("{:.3f}".format(x) for x in data) + "\t" + path

