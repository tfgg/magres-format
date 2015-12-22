from __future__ import print_function

import sys
from magres.atoms import MagresAtoms

# Load ethanol-jc.magres sample into an atoms structure
atoms = MagresAtoms.load_magres("../samples/ethanol-jc.magres")

# For fun, set all the hydrogens to be tritium
for atom in atoms.species('H'):
  atom.isotope = 3

# Set all the carbons to be 12C (won't work!)
try:
  for atom in atoms.species('C'):
    atom.isotope = 12
except ValueError as e:
  print("Error changing C isotopes:", e, file=sys.stderr)

# Loop over and print out coupling symbol, distance between atoms, isotropic coupling, anisotropic coupling and asymmetry
for isc in atoms.isc:
  print(isc.atom1, isc.atom2, isc.dist, isc.J_iso, isc.J_aniso, isc.J_eta)

