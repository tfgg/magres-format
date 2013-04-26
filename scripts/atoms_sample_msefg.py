import sys
from magres.atoms import MagresAtoms

atoms = MagresAtoms.load_magres(open(sys.argv[1]))

for atom in atoms:
  print atom.efg_isotope, atom.label, atom.index, atom.ms.iso, atom.efg.Cq
