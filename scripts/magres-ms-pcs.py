#!python
import sys
from magres.atoms import MagresAtoms

file = sys.argv[1]
atoms = MagresAtoms.load_magres(file)

print("# Principle components of sigma (magnetic shielding, unreferenced) tensor")

for atom in atoms:
  pcs = atom.ms.evals
  print(", ".join(map(str, [atom.species, atom.index, pcs[0], pcs[1], pcs[2]])))

