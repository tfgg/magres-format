from magres.atoms import MagresAtoms

atoms = MagresAtoms.load_magres('samples/simple.magres')

images = atoms.within(atoms[0], 3.0)

for image in images:
  print image, atoms[0].dist(image), image.position

images[0].species = images[0].label = 'K'

for image in images:
  print image, atoms[0].dist(image), image.position
