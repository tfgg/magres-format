import numpy
import math

import constants

from decorators import lazyproperty

class MagresAtom(object):
  __slots__ = ["magres_atom",
               "reference",
               "efg",
               "efg_nonlocal",
               "efg_local",
               "isc",
               "isc_fc",
               "isc_spin",
               "isc_orbital_p",
               "isc_orbital_d",
               "ms",
               "bonded",
               "_isotope",]

  def __init__(self, magres_atom):
    self.magres_atom = magres_atom
    self.reference = 0.0

  def __str__(self):
    if self.species != self.label:
      return "%d%s(%s)%d" % (self.isotope, self.species, self.label, self.index)
    else:
      return "%d%s%d" % (self.isotope, self.species, self.index)

  def dist(self, r):
    """
      Calculate distance from this atom to another position or atom.
    """

    if hasattr(r, 'position'):
      r = r.position
    dr = self.position - r
    return math.sqrt(numpy.dot(dr, dr))

  @property
  def label(self):
    """
      This atom's label.
    """
    return self.magres_atom['label']
  
  @label.setter
  def label(self, value):
    self.magres_atom['label'] = value

  @property
  def species(self):
    """
      This atom's species.
    """
    return self.magres_atom['species']
 
  @species.setter
  def species(self, value):
    self.magres_atom['species'] = value

  @property
  def index(self):
    """
      This atom's label index.
    """
    return self.magres_atom['index']
  
  @property
  def position(self):
    """
      This atom's position in cartesian coordinates. Units are Angstroms.
    """
    return numpy.array(self.magres_atom['position'])

  @property
  def isotope(self):
    """
      The isotope of this atom. Assumed to be most common NMR-active nucleus unless specified otherwise.
    """
    if hasattr(self, '_isotope'):
      return self._isotope
    else:
      if self.species in constants.gamma_iso:
        return constants.gamma_iso[self.species]
      else:
        return None

  @isotope.setter
  def isotope(self, value):
    if (self.species, value) in constants.gamma:
      self._isotope = value
    else:
      raise ValueError("Unknown NMR isotope %d%s" % (value, self.species))

  @property
  def gamma(self):
    """
      The gyromagnetic ratio constant of this atom's species and isotope.
    """
    if (self.species, self.isotope) in constants.gamma:
      return constants.gamma[(self.species, self.isotope)]
    else:
      return 0.0
  
  @property
  def Q(self):
    """
      The quadrupole moment of this atom's species and isotope.
    """
    if (self.species, self.isotope) in constants.Q:
      return constants.Q[(self.species, self.isotope)]
    else:
      return 0.0

  @property
  def spin(self):
    if (self.species, self.isotope) in constants.iso_spin:
      return constants.iso_spin[(self.species, self.isotope)]
    else:
      return None # We don't even know what spin this isotope is.

  def __eq__(self, other):
    idx1 = (self.species, self.index)

    if not (hasattr(other, 'species') and hasattr(other, 'index')):
      return False
    else:
      idx2 = (other.species, other.index)

      if idx1 == idx2:
        return True
      else:
        return False

  def __ne__(self, other):
    if not (self == other):
      return True
    else:
      return False


class MagresAtomImage(object):
  """
    A periodic image of a particular atom. Exactly like the underlying atom except for its position.
  """

  __slots__ = ["position", "atom"]

  def __init__(self, atom, position):
    object.__setattr__(self, "position", position)
    object.__setattr__(self, "atom", atom)

  def dist(self, r):
    if hasattr(r, 'position'):
      r = r.position

    dr = self.position - r
    return math.sqrt(numpy.dot(dr, dr))

  def __str__(self):
    return str(self.atom)
  
  def __unicode__(self):
    return unicode(self.atom)

  def __getattribute__(self, name):
    if name == "atom":
      return object.__getattribute__(self, "atom")
    elif name == "position":
      return object.__getattribute__(self, "position")
    elif name == "dist":
      return object.__getattribute__(self, "dist")
    else:
      return getattr(object.__getattribute__(self, "atom"), name)

  def __setattr__(self, name, value):
    if name == "position":
      object.__setattr__(self, "position", value)
    elif name == "atom":
      object.__setattr__(self, "atom", value)
    else:
      setattr(object.__getattribute__(self, "atom"), name, value)

