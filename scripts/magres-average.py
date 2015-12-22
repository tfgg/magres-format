#!python
import sys
from numpy import mean
from magres.atoms import MagresAtoms

# Load the magres files
out_atoms = MagresAtoms.load_magres(sys.argv[1])
atomss = [MagresAtoms.load_magres(f) for f in sys.argv[1:]]

for atom in out_atoms:
  other_atoms = [atoms.species(atom.species)[atom.index-1] for atoms in atomss]

  atom.magres_atom['position'] = mean([atom_.position for atom_ in other_atoms],0).tolist()

  if hasattr(atom, 'ms'):
    atom.ms.magres_ms['sigma'] = mean([atom_.ms.sigma for atom_ in other_atoms],0).tolist()

  if hasattr(atom, 'efg'):
    atom.efg.magres_efg['V'] = mean([atom_.efg.V for atom_ in other_atoms],0).tolist()
  
print(out_atoms.magres_file)

