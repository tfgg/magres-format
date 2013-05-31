#!python
import sys
from magres.format import MagresFile, BadVersion
from magres.atoms import MagresAtoms
from magres.utils import find_all_magres

magres_files = []
for f in find_all_magres(sys.argv[1]):
  try:
    magres_files.append(MagresFile(f))
  except BadVersion:
    continue

s1 = sys.argv[2]
i1 = int(sys.argv[3])
s2 = sys.argv[4]
i2 = int(sys.argv[5])

tensor = 'isc'
if len(sys.argv) >= 7:
  tensor = sys.argv[6]

atoms = MagresAtoms.load_magres(magres_files)

atom1 = atoms.get_species(s1, i1)
atom2 = atoms.get_species(s2, i2)

if tensor == "isc_spin_total":
  K1 = atom1.isc_fc[atom2].K + atom1.isc_spin[atom2].K
  K2 = atom2.isc_fc[atom1].K + atom2.isc_spin[atom1].K
else:
  K1 = getattr(atom1, tensor)[atom2].K
  K2 = getattr(atom2, tensor)[atom1].K

print "K tensors"
print " ".join(["%10s" % ("%.3f" % x) for x in K1[0]]) + "   " + " ".join(["%10s" % ("%.3f" % x) for x in K2[0]])
print " ".join(["%10s" % ("%.3f" % x) for x in K1[1]]) + "   " + " ".join(["%10s" % ("%.3f" % x) for x in K2[1]])
print " ".join(["%10s" % ("%.3f" % x) for x in K1[2]]) + "   " + " ".join(["%10s" % ("%.3f" % x) for x in K2[2]])

#print "K_sym difference, K_asym difference"
#delta_K_sym = isc1.K_sym - isc2.K_sym
#delta_K_asym = isc1.K_asym + isc2.K_asym

#print " ".join(["%10s" % ("%.3f" % x) for x in delta_K_sym[0]]) + "   " + " ".join(["%10s" % ("%.3f" % x) for x in delta_K_asym[0]])
#print " ".join(["%10s" % ("%.3f" % x) for x in delta_K_sym[1]]) + "   " + " ".join(["%10s" % ("%.3f" % x) for x in delta_K_asym[1]])
#print " ".join(["%10s" % ("%.3f" % x) for x in delta_K_sym[2]]) + "   " + " ".join(["%10s" % ("%.3f" % x) for x in delta_K_asym[2]])


