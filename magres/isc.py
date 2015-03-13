import constants
import html_repr
import numpy

from view import ListPropertyView

class MagresAtomIsc(object):
  """
    Representation of the indirect spin coupling between two atoms.
  """
  __slots__ = ["atom1", "atom2", "magres_isc"]

  def __init__(self, atom1, atom2, magres_isc):
    self.atom1 = atom1
    self.atom2 = atom2
    self.magres_isc = magres_isc

  def perturbing(self, species=None, index=None):
    if hasattr(species, "__call__"):
      if species(self.atom1):
        return ListPropertyView([self])
      else:
        return ListPropertyView([])

    elif hasattr(species, 'species'):
      index = species.index
      species = species.species

    if species is None or self.atom1.species == species:
      if index is None or self.atom1.index == index:
        return ListPropertyView([self])
      else:
        return ListPropertyView([])
    else:
      return ListPropertyView([])

  def receiving(self, species=None, index=None):
    if hasattr(species, "__call__"):
      if species(self.atom2):
        return ListPropertyView([self])
      else:
        return ListPropertyView([])

    elif hasattr(species, 'species'):
      index = species.index
      species = species.species

    if species is None or self.atom2.species == species:
      if index is None or self.atom2.index == index:
        return ListPropertyView([self])
      else:
        return ListPropertyView([])
    else:
      return ListPropertyView([])

  @property
  def symbol(self):
    """
      A textual symbol representing this coupling.
    """
    return "%s -> %s" % (self.atom1, self.atom2)

  @property
  def dist(self):
    """
      The distance between the two atoms involved.
    """
    return self.atom1.dist(self.atom2.position)

  @property
  def K(self):
    """
      The reduced indirect spin-spin coupling K tensor.
    """
    return numpy.array(self.magres_isc['K'])
 
  @K.setter
  def K(self, value):
    sh = numpy.shape(value)
    if sh != (3,3):
      raise Exception("Wrong shape for new indirect spin-spin coupling tensor, {}. Should be (3,3)".format(sh))

    self.magres_isc['K'] = value 

  @property
  def K_iso(self):
    """
      The isotropic component of the reduced indirect spin-spin coupling tensor.
    """
    return numpy.trace(self.K)/3.0

  @property
  def J(self):
    """
      The spin-spin coupling J tensor.
    """
    return constants.K_to_J_iso(self.K, self.atom1.species, self.atom1.isotope, self.atom2.species, self.atom2.isotope)

  @property
  def J_iso(self):
    """
      The isotropic component of the indirect spin-spin coupling J tensor.
    """
    return numpy.trace(self.J)/3.0

  @property
  def K_sym(self):
    """
      The symmetric component of the reduced indirect spin-spin coupling tensor K.
    """
    return (self.K + self.K.T)/2.0
  
  @property
  def K_asym(self):
    """
      The asymmetric component of the reduced indirect spin-spin coupling tensor K.
    """
    return (self.K - self.K.T)/2.0
  
  @property
  def J_sym(self):
    """
      The symmetric component of the indirect spin-spin coupling tensor J.
    """
    return (self.J + self.J.T)/2.0
  
  @property
  def J_asym(self):
    """
      The asymmetric component of the indirect spin-spin coupling tensor J.
    """
    return (self.J - self.J.T)/2.0

  @property
  def K_evalsvecs(self):
    evals, evecs = numpy.linalg.eig(self.K_sym)

    se = zip(*sorted(zip(evals, evecs), key=lambda (x,y): abs(x-self.K_iso)))

    return ([se[0][1], se[0][0], se[0][2]], [se[1][1], se[1][0], se[1][2]])

  @property
  def K_evecs(self):
    return self.K_evalsvecs[1]
  
  @property
  def K_evals(self):
    return self.K_evalsvecs[0]

  @property
  def J_evalsvecs(self):
    evals, evecs = numpy.linalg.eig(self.J_sym)

    se = zip(*sorted(zip(evals, evecs), key=lambda (x,y): abs(x-self.J_iso)))

    return ([se[0][1], se[0][0], se[0][2]], [se[1][1], se[1][0], se[1][2]])

  @property
  def J_evecs(self):
    return self.J_evalsvecs[1]
  
  @property
  def J_evals(self):
    return self.J_evalsvecs[0]

  @property
  def K_aniso(self):
    """
      The K anisotropy.
    """
    jx = self.K_evals
    return jx[2] - (jx[0] + jx[1])/2.0
  
  @property
  def J_aniso(self):
    """
      The J anisotropy.
    """
    jx = self.J_evals
    return jx[2] - (jx[0] + jx[1])/2.0

  @property
  def K_eta(self):
    """
      The K principal component asymmetry.
    """
    jx = self.K_evals
    return (jx[1] - jx[0]) / (jx[2] - sum(jx)/3.0)
  
  @property
  def J_eta(self):
    """
      The K principal component asymmetry.
    """
    jx = self.J_evals
    return (jx[1] - jx[0]) / (jx[2] - sum(jx)/3.0)

  def _repr_html_(self):
    return html_repr.isc(self)


