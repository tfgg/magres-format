from magres.format import MagresFile
import sys

magres_file = MagresFile(sys.argv[1])
atoms = magres_file.as_ase()

print atoms

for atom in atoms:
  print atom.symbol, atom.tag, atom.position
