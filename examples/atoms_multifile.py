import sys
from magres.format import MagresFile
from magres.utils import find_all_magres
from magres.atoms import MagresAtoms

magres_files = [MagresFile(f) for f in find_all_magres(sys.argv[1])]

atoms = MagresAtoms.load_magres(magres_files)

coupling_atom = atoms.get_species('C', 1)

# Set a bunch of arbitrary references for the magnetic shielding
atoms.set_reference('H', 100.0)
atoms.set_reference('C', 55.0)
atoms.set_reference('O', -2.0)

for atom in atoms:
  print(atom, atom.efg.Cq, atom.ms.iso, atom.isc[coupling_atom].J_iso)

