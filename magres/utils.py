import os
from format import BadVersion
from atoms import MagresAtoms

def get_numeric(s):
  """
    Turn a path, e.g. "calcs_gs2fg4/damp_scale=0.5/x=0.02/jc_site=C_1/" into
    a whitespace deliminated sequence of numbers e.g. ["2","4","0.5","0.02,"1"].

    This is useful for plotting things like convergence.
  """
  return [float(x) for x in re.split("[^0-9\.]+", s) if len(x) != 0]

def find_all(dir, suffix=".cell"):
  """
    Recursively find all files with a particular suffix starting in directory dir. Returns a list of the relative file paths.

    >>> from magres.utils import find_all
    >>> find_all('.', '.magres')
    ['./TlCl/TlCl.magres', './Tl(CN)Cl2/Tl(CN)Cl2.magres', './Tl4(OCH3)4/Tl4(OCH3)4.magres', './TlF/TlF.magres', './Tl(CN)2Cl/Tl(CN)2Cl.magres', './TlBr/TlBr.magres', './TlI/TlI.magres', './Tl(CN)3/Tl(CN)3.magres']

  """

  calcs = []
  for f in os.listdir(dir):
    path = os.path.join(dir, f)
    if f.endswith(suffix):
      calcs.append(path)
    elif os.path.isdir(path):
      calcs += find_all(path, suffix)
  
  return calcs

def find_all_magres(dir):
  """
    Find all magres files starting in directory dir.
  """
  calcs = []
  for f in os.listdir(dir):
    path = os.path.join(dir, f)
    if ".magres" in f:
      calcs.append(path)
    elif os.path.isdir(path):
      calcs += find_all_magres(path)
  
  return calcs

def load_all_magres(dir):
  """
    Find all magres files starting in directory dir and load them into a :py:class:`magres.atoms.MagresAtoms` structure. Returns a list.
  """

  atoms = []
  for magres_file in find_all_magres(dir):
    try:
      atoms.append(MagresAtoms.load_magres(magres_file))
    except BadVersion:
      pass

  return atoms


