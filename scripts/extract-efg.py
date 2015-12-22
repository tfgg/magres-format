#!python

import os
import sys
import argparse

from magres.atoms import MagresAtoms
from magres.utils import load_all_magres, get_numeric, parse_atom_list

parser = argparse.ArgumentParser(description='Extract and print EFG values from a set of calculations.')
parser.add_argument('-N', '--numbers', action="store_const", help="Parse numbers from path and print. This is useful for e.g. convergence calculations.", default=False, const=True)
parser.add_argument('source', help='Directory to look for calculations below or specific file.')

parser.add_argument('atoms', nargs='?', type=str, default=None, help='Which atoms to print shieldings of. Specify with atom list notation, e.g. "H1" or "H1,H2,H3" or "H,C" or "H1-3".')

a = parser.parse_args(sys.argv[1:])

atoms_filter_str = a.atoms

if atoms_filter_str:
  atoms_filter = parse_atom_list(atoms_filter_str)
else:
  atoms_filter = lambda x: True

tensors = ['efg',]# 'efg_local', 'efg_nonlocal']

lines = []

if a.numbers:
  print("# Number\tAtom\tCq\tEta\tPath")
else:
  print("# Atom\tCq\tEta\tPath")

if os.path.isfile(a.source):
  magres_atoms = [MagresAtoms.load_magres(a.source)]
else:
  magres_atoms = load_all_magres(a.source)

for i, atoms in enumerate(magres_atoms):
  num = get_numeric(atoms.magres_file.path)

  if num:
    idx = num
  else:
    idx = [i]

  for atom in atoms: 
    if atoms_filter(atom) and \
      hasattr(atom, 'efg'):

      lines.append((idx,
                    atoms.magres_file.path,
                    str(atom),
                    ["%.2f" % getattr(atom, tensor).Cq for tensor in tensors],
                    ["%.2f" % getattr(atom, tensor).eta for tensor in tensors]))

lines = sorted(lines, key=lambda xs: xs[0])

for idx, path, atom, data1, data2 in lines:
  if a.numbers:
    print(" ".join(map(str,idx)) + "\t" + atom + "\t" + "\t".join(data1) + "\t" + "\t".join(data2), path)
  else:
    print(atom + "\t" + "\t".join(data1) + "\t" + "\t".join(data2), path)


