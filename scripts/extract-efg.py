#!python

import os
import sys
import argparse

from magres.utils import load_all_magres, get_numeric

parser = argparse.ArgumentParser(description='Extract J-coupling parameters in both directions')
parser.add_argument('-J', '--J_tensor', action="store_const", help="Display J tensor", default=False, const=True)
parser.add_argument('source_dir', help='Directory to look for calculations in')
parser.add_argument('species', nargs=argparse.REMAINDER, help='Species to look for')

a = parser.parse_args(sys.argv[1:])

magres_atoms = load_all_magres(a.source_dir)

find_s = str(a.species[0])

if len(a.species) >= 2:
  find_i = int(a.species[1])
else:
  find_i = None

tensors = ['efg', 'efg_local', 'efg_nonlocal']

lines = []

print "# Number\tAtom\tCq\tCq_local\tCq_nonlocal\tEta\tEta_local\tEta_nonlocal\tPath"

for atoms in magres_atoms:
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
                    ["%.2f" % getattr(atom, tensor).Cq for tensor in tensors],
                    ["%.2f" % getattr(atom, tensor).eta for tensor in tensors]))

lines = sorted(lines, key=lambda xs: xs[0])

for idx, path, atom, data1, data2 in lines:
  print " ".join(map(str,idx)) + "\t" + atom + "\t" + "\t".join(data1) + "\t" + "\t".join(data2), path


