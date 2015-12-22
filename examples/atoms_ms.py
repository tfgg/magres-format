from __future__ import print_function

import sys
from magres.atoms import MagresAtoms

atoms = MagresAtoms.load_magres("../samples/ethanol-nmr.magres")

for atom in atoms:
  print(atom, atom.ms.iso, atom.ms.aniso, atom.ms.zeta, atom.ms.eta, atom.ms.span, atom.ms.skew, atom.ms.evals[0], atom.ms.evals[1], atom.ms.evals[2])

