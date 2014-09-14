#!python

"""
  Only show couplings where the coupling in the opposite direction exists in another calculation.
"""

import os
import sys
import argparse

from magres.format import BadVersion
from magres.atoms import MagresAtoms
from magres.utils import load_all_magres

cwd = "."

parser = argparse.ArgumentParser(description='Extract J-coupling parameters in both directions')
parser.add_argument('-J', '--J_tensor', action="store_const", help="Display J tensor", default=False, const=True)
parser.add_argument('source_dir', help='Directory to look for calculations in')
parser.add_argument('atom_species1', nargs='?', type=str, default=None, help='Only print couplings from this atomic species.')
parser.add_argument('atom_index1', nargs='?', type=int, default=None, help='Only print couplings from this atom.')
parser.add_argument('atom_species2', nargs='?', type=str, default=None, help='Only print couplings to this atomic species.')
parser.add_argument('atom_index2', nargs='?', type=int, default=None, help='Only print couplings to this atom.')

a = parser.parse_args(sys.argv[1:])

find_s1 = a.atom_species1
find_i1 = a.atom_index1

find_s2 = a.atom_species2
find_i2 = a.atom_index2

all_Js = {}

tensors = ['isc', 'isc_fc', 'isc_spin', 'isc_orbital_p', 'isc_orbital_d']

magres_atoms = load_all_magres(a.source_dir)

for atoms in magres_atoms:
  have_all_tensors = True
  for tensor in tensors:
    if not hasattr(atoms, tensor):
      print "# J-coupling %s not found" % (tensor,)
      have_all_tensors = False

  if not have_all_tensors:
    continue

  for tensor in tensors:
    for atom in atoms:
      if hasattr(atom, tensor):
        for isc in getattr(atom, tensor).values():
          if (find_s1 is None and find_i1 is None) or \
             (isc.atom2.species == find_s1 and isc.atom2.index == find_i2) or \
             (isc.atom1.species == find_s1 and isc.atom1.index == find_i2) and \
             (find_s2 is None or isc.atom2.species == find_s2):

            idx = (isc.atom1.species,isc.atom1.index,isc.atom2.species,isc.atom2.index)
            if idx not in all_Js:
              all_Js[idx] = {}

            all_Js[idx][tensor] = isc

matching_Js = []

for (s1,i1,s2,i2),iscs in all_Js.items():
  if not (s2 == s1 and i2 == i1) and (s2,i2,s1,i1) in all_Js:
    matching_Js.append((s1,i1,s2,i2,iscs))

matching_Js = sorted(matching_Js, key=lambda (s1,i1,s2,i2,isc): sorted(((s1,i1),(s2,i2))))

print "#atm1\tatm2\t", "\t".join(tensors), "\tr"

for s1,i1,s2,i2,iscs in matching_Js:
  atom1 = iscs['isc'].atom1
  atom2 = iscs['isc'].atom2

  d = atom1.dist(atom2)

  J_isos = []
  K_isos = []
  for tensor in tensors:
    J_isos.append(iscs[tensor].J_iso)
    K_isos.append(iscs[tensor].K_iso)

  if a.J_tensor:
    print "J %s\t" % atom1 + "%s\t" % atom2 + "\t".join(["%.2f" % J_iso for J_iso in J_isos]) + "\t" + str(d)
  else:
    print "K %s\t" % atom1 + "%s\t" % atom2 + "\t".join(["%.2f" % K_iso for K_iso in K_isos]) + "\t" + str(d)

