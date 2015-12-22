import numpy

from . import html_repr

#from decorators import property

class MagresAtomMs(object):
  """
    Representation of the magnetic shielding of a particular atom.
  """

  __slots__ = ["atom", "magres_ms"]

  def __init__(self, atom, magres_ms):
    self.atom = atom
    self.magres_ms = magres_ms

  @property
  def sigma(self):
    """
      The sigma tensor, i.e. the magnetic shielding.
    """
    return numpy.array(self.magres_ms['sigma'])

  @sigma.setter
  def sigma(self, value):
    sh = numpy.shape(value)
    if sh != (3,3):
      raise Exception("Wrong shape for new magnetic shielding tensor, {}. Should be (3,3)".format(sh))

    self.magres_ms['sigma'] = value

  @property
  def sym(self):
    """
      The symmetric part of sigma.
      
      :math:`\sigma_{sym} = (\sigma + \sigma^T)/2`
    """
    return (self.sigma + self.sigma.T)/2.0
  
  @property
  def asym(self):
    """
      The asymmetric part of sigma.

      :math:`\sigma_{asym} = (\sigma - \sigma^T)/2`
    """
    return (self.sigma - self.sigma.T)/2.0

  @property
  def iso(self):
    """
      The isotropic part of sigma. Defined by

      :math:`\sigma_{iso} = (\sigma_{XX} + \sigma_{YY} + \sigma_{ZZ})/3`

    """
    return numpy.trace(self.sigma)/3.0

  @property
  def cs(self):
    """
      The chemical shift, referenced.
    """

    return self.atom.reference - self.iso

  @property
  def aniso(self):
    """
      The shielding anisotropy. Defined by

      :math:`\Delta \sigma = \sigma_{ZZ} - (\sigma_{XX} + \sigma_{YY})/2`
    """
    ev = self.evals
    return ev[2] - (ev[0] + ev[1])/2.0

  @property
  def zeta(self):
    """
      The shielding anisotropy (alternative). Defined by

      :math:`\zeta = \sigma_{ZZ} - \sigma_{iso}`
    """
    return self.evals[2] - self.iso

  @property
  def eta(self):
    """
      The shielding asymmetry. Defined by

      :math:`\eta = (\sigma_{YY} - \sigma_{XX}) / \zeta`
    """
    return (self.evals[1] - self.evals[0]) / self.zeta

  @property
  def evalsvecs(self):
    """
      The eigenvalues and eigenvectors of the symmetric part of sigma, ordered according to the Haeberlen convention:

      :math:`|\sigma_{ZZ} - \sigma_{iso}| \geq |\sigma_{XX} - \sigma_{iso}| \geq |\sigma_{YY} - \sigma_{iso}|`
      
      where

        sigma_XX = evals[0]
        sigma_YY = evals[1]
        sigma_ZZ = evals[2]
    """

    evals, evecs = numpy.linalg.eig(self.sym)

    se = list(zip(*sorted(zip(evals, evecs), key=lambda x_y: abs(x_y[0] - self.iso))))

    return ([se[0][1], se[0][0], se[0][2]], [se[1][1], se[1][0], se[1][2]])

  @property
  def evecs(self):
    """
      The eigenvectors of sigma, ordered according to the Haeberlen convention:

      :math:`|\sigma_{ZZ} - \sigma_{iso}| \geq |\sigma_{XX} - \sigma_{iso}| \geq |\sigma_{YY} - \sigma_{iso}|`

      where

        sigma_XX = evals[0]
        sigma_YY = evals[1]
        sigma_ZZ = evals[2]
    """
    return self.evalsvecs[1]
  
  @property
  def evals(self):
    """
      The eigenvalues of sigma, ordered according to the Haeberlen convention:

      :math:`|\sigma_{ZZ} - \sigma_{iso}| \geq |\sigma_{XX} - \sigma_{iso}| \geq |\sigma_{YY} - \sigma_{iso}|`

      where

        sigma_XX = evals[0]
        sigma_YY = evals[1]
        sigma_ZZ = evals[2]
    """
    return self.evalsvecs[0]

  @property
  def evalsvecs_mehring(self):
    """
      The eigenvalues and eigenvectors of the symmetric part of sigma ordered according to the Mehring notation:

      :math:`\sigma_{11} \leq \sigma_{22} \leq \sigma_{33}`
    """
 
    evals, evecs = numpy.linalg.eig(self.sym)

    se = list(zip(*sorted(zip(evals, evecs), key=lambda x_y1: x_y1[0])))

    return ([se[0][0], se[0][1], se[0][2]], [se[1][0], se[1][1], se[1][2]])
  
  @property
  def evals_mehring(self):
    """
      The eigenvalues of sigma ordered according to the Mehring notation:

      :math:`\sigma_{11} \leq \sigma_{22} \leq \sigma_{33}`
    """
    return self.evalsvecs_mehring[0]

  @property
  def evecs_mehring(self):
    """
      The eigenvectors of sigma ordered according to the Mehring notation:

      :math:`\sigma_{11} \leq \sigma_{22} \leq \sigma_{33}`
    """
    return self.evalsvecs_mehring[1]

  @property
  def span(self):
    """
      The span of sigma. Defined by

      :math:`\Omega = \sigma_{33} - \sigma_{11}`
    """
    return self.evals_mehring[2] - self.evals_mehring[0]
  
  @property
  def skew(self):
    """
      The skew of sigma. Defined by

      :math:`\kappa = 3(\sigma_{iso} - \sigma_{22}) / \Omega`
    """
    return 3.0*(self.iso - self.evals_mehring[1]) / self.span

  def _repr_html_(self):
    return html_repr.ms(self)
  
  def _repr_html_row_(self):
    return html_repr.ms_row(self)

  def _repr_html_row_head_(self):
    return html_repr.ms_row_head(self)

