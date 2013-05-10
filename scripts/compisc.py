#!python
import sys
from magres.format import MagresFile
from magres.atoms import MagresAtoms
from magres.utils import find_all_magres

magres_files = [MagresFile(f) for f in find_all_magres(sys.argv[1])]

s1 = sys.argv[2]
i1 = int(sys.argv[3])
s2 = sys.argv[4]
i2 = int(sys.argv[5])

atoms = MagresAtoms.load_magres(magres_files)

atom1 = atoms.get_species(s1, i1)
atom2 = atoms.get_species(s2, i2)

isc1 = atom1.isc[atom2]
isc2 = atom2.isc[atom1]

print "K tensors"
print isc1.K[0], isc2.K[0]
print isc1.K[1], isc2.K[1]
print isc1.K[2], isc2.K[2]

print "K_sym difference, K_asym difference"
delta_K_sym = isc1.K_sym - isc2.K_sym
delta_K_asym = isc1.K_asym - isc2.K_asym

print delta_K_sym[0], delta_K_asym[0]
print delta_K_sym[1], delta_K_asym[1]
print delta_K_sym[2], delta_K_asym[2]


