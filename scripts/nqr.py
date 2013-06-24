#!python
"""
  Script to calculate NQR transition frequencies from .magres ab-initio EFG calculations.

  Equation from "Nuclear Quadrupole Resonance Spectroscopy" by Hand and Das, Solid State Physics Supplement 1.

    A = eVzzQ / 4S(2S-1)
    omega = 3A(2|m| + 1)/hbar

  for all available transitions m->m+1

"""
import sys
import math
from magres.atoms import MagresAtoms
from magres.utils import load_all_magres
import magres.units as units

s = sys.argv[1]

atomss = []
for f in sys.argv[2:]:
  atomss += load_all_magres(f)

def frange(a, b, x):
  while a < b:
    yield a
    a += x

for atoms in atomss:
  print atoms.magres_file.path

  vals = []

  first = True

  for atom in atoms:
    if atom.species == s:
      if first:
        print "%d%s Q=%f mb" % (atom.isotope, atom.species, atom.Q)
        first = False

      ms = [m for m in frange(-atom.spin, atom.spin+1, 1) if m >= 0.0]

      for m in ms[:-1]:
        Vzz = atom.efg.evals[2]

        A = Vzz * (atom.Q * units.millibarn) / (4.0 * atom.spin * (2.0*atom.spin - 1.0))
        fq = 3*A * (2.0*abs(m) + 1.0)

        vals.append(fq)
        print atom, "m=%.1f-->%.1f" % (m,m+1.0), "f=%f MHz" % (fq / units.megahertz)
     
  print "mean=%f MHz" % (sum(vals)/len(vals) / units.megahertz)
 
