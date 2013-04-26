import sys
from magres.atoms import MagresAtoms

atoms = MagresAtoms.load_magres("samples/ethanol-jc.magres")

# Set all the hydrogens to be tritium
for atom in atoms.get_species('H'):
  atom.isc_isotope = 3

# Set all the carbons to be 12C (won't work!)
try:
  for atom in atoms.get_species('C'):
    atom.isc_isotope = 12
except ValueError, e:
  print >>sys.stderr, "Error changing C isotopes:", e

# Loop over and print out all couplings
for isc in atoms.isc:
  print isc.symbol, isc.dist, isc.K_iso, isc.K_aniso, isc.K_eta

