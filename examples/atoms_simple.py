from magres.atoms import MagresAtoms

atoms = MagresAtoms.load_magres('samples/simple.magres')

print "5-1"
for atom in atoms.within(atoms[0], 5.0).within(atoms[0], 1.0):
  print atom, atom.position

print "1-5"
for atom in atoms.within(atoms[0], 1.0).within(atoms[0], 5.0):
  print atom, atom.position


