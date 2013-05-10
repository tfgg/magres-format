from magres.atoms import MagresAtoms

atoms = MagresAtoms.load_magres('samples/simple.magres')

images = atoms.within(atoms[0], 5.0)

for image in images:
  print image, atoms[0].dist(image), " ".join(map(str,image.position))

