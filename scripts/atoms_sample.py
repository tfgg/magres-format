import sys
from magres.atoms import MagresAtoms

atoms = MagresAtoms.load_magres(open(sys.argv[1]))

for isc in atoms.isc:
  print isc.atom1.label, isc.atom1.index, isc.atom2.label, isc.atom2.index, isc.K_iso

