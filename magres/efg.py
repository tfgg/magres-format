import numpy
from . import constants
from . import html_repr
#from decorators import property

class MagresAtomEfg(object):
  """
    Representation of the electric field gradient on a particular atom.
  """

  __slots__ = ["atom", "magres_efg"]

  def __init__(self, atom, magres_efg):
    self.atom = atom
    self.magres_efg = magres_efg

  @property
  def V(self):
    """
      The EFG V tensor in atomic units.
    """
    return numpy.array(self.magres_efg['V'])

  @V.setter
  def V(self, value):
    sh = numpy.shape(value)
    if sh != (3,3):
      raise Exception("Wrong shape for new electric field gradient tensor, {}. Should be (3,3)".format(sh))

    self.magres_efg['V'] = value

  @property
  def Cq(self):
    """
      The Cq of the V tensor.

      :math:`C_q = eV_{ZZ} Q / h`

      Where Q is the quadrupole moment of this particular species, e is the electron charge and h is Planck's constant.
    """
    try:
      return constants.efg_to_Cq_isotope(self.V, self.atom.species, self.atom.isotope)
    except KeyError:
      return 0.0

  @property
  def eta(self):
    """
      The asymmetry of the V tensor

      :math: `(V_{XX}-V_{YY})/V_{ZZ}`
    """
    evals = self.evals

    return (evals[1] - evals[0])/evals[2]

  @property
  def evalsvecs(self):
    """
      The eigenvalues and eigenvectors of V, ordered according to the Haeberlen convention:

      :math:`|V_{ZZ}| \geq |V_{XX}| \geq |V_{YY}|`

      where

        :math:`V_{XX}` = evals[0]

        :math:`V_{YY}` = evals[1]

        :math:`V_{ZZ}` = evals[2]
    """
    evals, evecs = numpy.linalg.eig(self.magres_efg['V'])

    se = list(zip(*sorted(zip(evals, evecs), key=lambda x_y: abs(x_y[0]))))

    return ([se[0][1], se[0][0], se[0][2]], [se[1][1], se[1][0], se[1][2]])

  @property
  def evecs(self):
    """
      The eigenvectors of V, ordered according to the Haeberlen convention:

      :math:`|V_{ZZ}| \geq |V_{XX}| \geq |V_{YY}|`

      where

        :math:`V_{XX}` = evals[0]

        :math:`V_{YY}` = evals[1]

        :math:`V_{ZZ}` = evals[2]
    """
    return self.evalsvecs[1]

  @property
  def evals(self):
    """
      The eigenvalues of V, ordered according to the Haeberlen convention:

      :math:`|V_{ZZ}| \geq |V_{XX}| \geq |V_{YY}|`

      where

        :math:`V_{XX}` = evals[0]

        :math:`V_{YY}` = evals[1]

        :math:`V_{ZZ}` = evals[2]
    """
    return self.evalsvecs[0]


  def _repr_html_(self):
    return html_repr.efg(self)

  def _repr_html_row_(self):
    return html_repr.efg_row(self)

  def _repr_html_row_head_(self):
    return html_repr.efg_row_head(self)
