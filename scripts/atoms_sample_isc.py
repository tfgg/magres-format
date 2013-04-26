import sys
from magres.atoms import MagresAtoms

atoms = MagresAtoms.load_magres(open(sys.argv[1]))

# Set all the hydrogens to be tritium
for atom in atoms.get_species('H'):
  atom.isc_isotope = 3

# Loop over and print out all couplings
for isc in atoms.isc:
  print "%d%s" % (isc.atom1.isc_isotope, isc.atom1), "%d%s" % (isc.atom2.isc_isotope, isc.atom2), isc.dist, isc.K_iso, isc.K_aniso, isc.K_eta

