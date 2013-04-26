import sys
from magres.atoms import MagresAtoms

# Load ethanol-jc.magres sample into an atoms structure
atoms = MagresAtoms.load_magres("samples/ethanol-jc.magres")

# For fun, set all the hydrogens to be tritium
for atom in atoms.get_species('H'):
  atom.isc_isotope = 3

# Set all the carbons to be 12C (won't work!)
try:
  for atom in atoms.get_species('C'):
    atom.isc_isotope = 12
except ValueError, e:
  print >>sys.stderr, "Error changing C isotopes:", e

# Loop over and print out coupling symbol, distance between atoms, isotropic coupling, anisotropic coupling and asymmetry
for isc in atoms.isc:
  print isc.symbol, isc.dist, isc.J_iso, isc.J_aniso, isc.J_eta

