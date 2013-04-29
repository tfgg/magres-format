from magres.atoms import MagresAtoms

atoms = MagresAtoms.load_magres('samples/simple.magres')

for atom_image in atoms.within(atoms[0], 3.0):
  print atoms[0].dist(atom_image), atom_image.position
