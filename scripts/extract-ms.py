#!python

import os
import sys

from magres.format import BadVersion
from magres.atoms import MagresAtoms
from magres.utils import load_all_magres

cwd = "."

if len(sys.argv)>1:
    cwd = str(sys.argv[1])

magres_atoms = load_all_magres(cwd)

find_s = str(sys.argv[2])

if len(sys.argv) >= 4:
  find_i = int(sys.argv[3])
else:
  find_i = None

tensors = ['efg', 'efg_local', 'efg_nonlocal']

lines = []

def get_numeric(s):
  return "".join([c for c in s if ord("0") <= ord(c) <= ord("9")])

print "# Number\tPath\tAtom\tCq\tCq_local\tCq_nonlocal\tEta\tEta_local\tEta_nonlocal"

for atoms in magres_atoms:
  num = get_numeric(atoms.magres_file.path)

  if num != '':
    idx = int(num)
  else:
    idx = 0

  for atom in atoms: 
    if atom.species == find_s and (find_i is None or atom.index == find_i) and hasattr(atom, 'efg'):
      lines.append((idx, atoms.magres_file.path, str(atom) + "\t" + "\t".join(["%.2f" % getattr(atom, tensor).Cq for tensor in tensors]) + "\t" + "\t".join(["%.2f" % getattr(atom, tensor).eta for tensor in tensors])))

lines = sorted(lines, key=lambda (x,y,z): x)

for idx, path, line in lines:
  print idx, path, line


