from __future__ import print_function
from magres.atoms import MagresAtoms

atoms = MagresAtoms.load_magres('../samples/simple.magres')

print(atoms[0].position)

Catoms = atoms.species('C')

print("1")
for atom in Catoms.within(atoms[0], 1.0):
    print(atom, atom.position)

print("5")
for atom in Catoms.within(atoms[0], 3.0):
    print(atom, atom.position)

print("5-1")
for atom in Catoms.within(atoms[0], 3.0).within(atoms[0], 1.0):
    print(atom, atom.position)

print("1-5")
for atom in Catoms.within(atoms[0], 1.0).within(atoms[0], 3.0):
    print(atom, atom.position)
