import os
from format import BadVersion
import atoms

def find_all_magres(dir):
  calcs = []
  for f in os.listdir(dir):
    path = os.path.join(dir, f)
    if ".magres" in f:
      calcs.append(path)
    elif os.path.isdir(path):
      calcs += find_all_magres(path)
  
  return calcs

def load_all_magres(dir):
  atoms = []
  for magres_file in find_all_magres(dir):
    try:
      atoms.append(MagresAtoms.load_magres(magres_file))
    except BadVersion:
      pass

  return atoms
