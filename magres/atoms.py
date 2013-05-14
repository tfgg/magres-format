"""
  magres.atoms is a collection of user-friendly data structures for representing
  groups of atoms and their NMR parameters. The NMR parameter structures can
  calculate a variety of properties from the underlying tensors, according to
  a variety of conventions.
"""

import sys,os
import math
import constants
import numpy

import utils
from utils import lazyproperty

from format import MagresFile

class MagresAtomEfg(object):
  """
    Representation of the electric field gradient on a particular atom.
  """

  def __init__(self, atom, magres_efg):
    self.atom = atom
    self.magres_efg = magres_efg

  @lazyproperty
  def V(self):
    """
      Return the EFG V tensor in atomic units (?)
    """
    return numpy.array(self.magres_efg['V'])

  @lazyproperty
  def Cq(self):
    """
      Calculate the Cq from the V tensor
    """
    try:
      return constants.efg_to_Cq_isotope(self.V, self.atom.species, self.atom.isotope)
    except KeyError:
      return 0.0

  @lazyproperty
  def evalsvecs(self):
    evals, evecs = numpy.linalg.eig(self.magres_efg['V'])

    se = zip(*sorted(zip(evals, evecs), key=lambda (x,y): abs(x)))

    return ([se[0][1], se[0][0], se[0][2]], [se[1][1], se[1][0], se[1][2]])

  @lazyproperty
  def evecs(self):
    return self.evalsvecs[1]
  
  @lazyproperty
  def evals(self):
    return self.evalsvecs[0]

class MagresAtomIsc(object):
  """
    Representation of the indirect spin coupling between two atoms.
  """
  def __init__(self, atom1, atom2, magres_isc):
    self.atom1 = atom1
    self.atom2 = atom2
    self.magres_isc = magres_isc

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

  @lazyproperty
  def K_evalsvecs(self):
    evals, evecs = numpy.linalg.eig(self.K)

    se = zip(*sorted(zip(evals, evecs), key=lambda (x,y): abs(x-self.K_iso)))

    return ([se[0][1], se[0][0], se[0][2]], [se[1][1], se[1][0], se[1][2]])

  @lazyproperty
  def K_evecs(self):
    return self.K_evalsvecs[1]
  
  @lazyproperty
  def K_evals(self):
    return self.K_evalsvecs[0]

  @lazyproperty
  def J_evalsvecs(self):
    evals, evecs = numpy.linalg.eig(self.J)

    se = zip(*sorted(zip(evals, evecs), key=lambda (x,y): abs(x-self.J_iso)))

    return ([se[0][1], se[0][0], se[0][2]], [se[1][1], se[1][0], se[1][2]])

  @lazyproperty
  def J_evecs(self):
    return self.J_evalsvecs[1]
  
  @lazyproperty
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

class MagresAtomMs(object):
  """
    Representation of the magnetic shielding of a particular atom.
  """

  def __init__(self, atom, magres_ms, reference=None):
    self.atom = atom
    self.magres_ms = magres_ms
    self.reference = reference

  @property
  def sigma(self):
    """
      The sigma tensor.
    """
    return self.magres_ms['sigma']

  @property
  def iso(self):
    """
      The referenced isotropic component of the sigma tensor.
    """

    if self.reference is None:
      reference = 0.0
    else:
      reference = self.reference

    return reference - numpy.trace(self.sigma)/3.0

class MagresAtom(object):
  def __init__(self, magres_atom):
    self.magres_atom = magres_atom

  def __str__(self):
    if self.species != self.label:
      return "%d%s(%s)%d" % (self.isotope, self.species, self.label, self.index)
    else:
      return "%d%s%d" % (self.isotope, self.species, self.index)

  def dist(self, r):
    if hasattr(r, 'position'):
      r = r.position
    dr = self.position - r
    return math.sqrt(numpy.dot(dr, dr))

  @property
  def label(self):
    return self.magres_atom['label']
  
  @label.setter
  def label(self, value):
    self.magres_atom['label'] = value

  @property
  def species(self):
    return self.magres_atom['species']
 
  @species.setter
  def species(self, value):
    self.magres_atom['species'] = value

  @property
  def index(self):
    return self.magres_atom['index']
  
  @property
  def position(self):
    return numpy.array(self.magres_atom['position'])

  @property
  def isotope(self):
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
    if self.isotope in constants.gamma:
      return constants.gamma[self.isotope]
    else:
      return 0.0
  
  @property
  def Q(self):
    if self.isotope in constants.Q:
      return constants.Q[self.isotope]
    else:
      return 0.0

class MagresAtomImage(object):
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

class LabelNotFound(Exception):
  pass

class SpeciesNotFound(Exception):
  pass

class AtomNotFound(Exception):
  pass

class MagresAtomsView(object):
  def __init__(self, atoms=None, lattice=None):
    if atoms is not None:
      self.atoms = atoms
    else:
      self.atoms = []

    if lattice is not None:
      self.lattice = lattice
    else:
      self.lattice = None

    self.label_index = {}
    self.species_index = {}

    self.build_index()

  def add(self, atoms):
    """
      Add an extra atom or list of atoms.
    """

    if type(atoms) == Atom:
      self.atoms.append(atoms)
    elif type(atom) == list:
      for atom in atoms:
        self.atoms.append(atom)

    self.build_index()

  def build_index(self):
    self.label_index.clear()
    self.species_index.clear()

    for atom in self.atoms:
      if atom.label in self.label_index:
        self.label_index[atom.label].append(atom)
      else:
        self.label_index[atom.label] = [atom]
      
      if atom.species in self.species_index:
        self.species_index[atom.species].append(atom)
      else:
        self.species_index[atom.species] = [atom]

  def get_label(self, label, index=None):
    """
      Get all atoms of a particular label or a single atom of a particular label and index.

      >>> atoms.get_label("Al2")
      >>> atoms.get_label("Al2", 10)
    """

    if label not in self.label_index:
      return []
    else:
      if index is None:
        return self.label_index[label]
      elif len(self.label_index[label]) >= index:
        return self.label_index[label][index-1]
      else:
        raise AtomNotFound("Atom %s %d does not exist in this system. There are %d atoms at the %s label." % (label, index, len(self.label_index[label]), label))

  def label(self, label):
    """
      test version of get_label that returns the list of atoms as a magresatoms object, to make chaining easy.
    """
    if type(label) != list:
      label = [label]
      
    rtn_atoms = []
    for l in label:
      if l in self.label_index:
        rtn_atoms += self.label_index[l]
    return MagresAtomsView(rtn_atoms, self.lattice)

  def get_species(self, species, index=None):
    """
      Get all atoms of a particular species or a single atom of a particular species and index.
      
      >>> atoms.get_species("Al")
      >>> atoms.get_species("Al", 10)
    """

    if species not in self.species_index:
      return []
    else:
      if index is None:
        return self.species_index[species]
      elif len(self.species_index[species]) >= index:
        return self.species_index[species][index-1]
      else:
        raise AtomNotFound("Atom %s %d does not exist in this system. There are %d atoms in the %s species." % (species, index, len(self.species_index[species]), species))

  def species(self, species):
    """
      test version of get_species that returns the list of atoms as a magresatoms object, to make chaining easy.
    """
    if type(species) != list:
      species = [species]
      
    rtn_atoms = []
    for s in species:
      if s in self.species_index:
        rtn_atoms += self.species_index[s]
    return MagresAtomsView(rtn_atoms, self.lattice)

  def within(self, pos, max_dr):
    """
      Gives you all atoms within max_dr angstroms of pos, including all images.
    """

    if type(pos) is MagresAtom:
      pos = pos.position

    atoms = []

    for atom in self.atoms:
      if type(atom) == MagresAtom:
        images = self.all_images_within(atom.position, pos, max_dr)

        for image_dist, image_pos in images:
          if image_dist <= max_dr:
            atoms.append(MagresAtomImage(atom, image_pos))
      elif type(atom) == MagresAtomImage:
        if atom.dist(pos) <= max_dr:
          atoms.append(atom)

    return MagresAtomsView(atoms, self.lattice)

  #def transform(self, M):
  #  """
  #    Apply matrix M to positions of all atoms and return them.
  #
  #    TODO: Is this useful? What about rotating the tensors?
  #    Active vs. passive rotations? What about the lattice?
  #  """

  #  return [MagresAtomImage(numpy.dot(M, atom.position.T).T, atom) for atom in self.atoms]

  def least_mirror(self, a, b):
    """
      Give the closest periodic image of a to b given the current lattice.
    """
    min = None
    min_p = None

    for i in range(-1,2):
      for j in range(-1,2):
        for k in range(-1,2):
          ap = numpy.add(a, numpy.dot(self.lattice.T, (float(i), float(j), float(k))))
          r = numpy.subtract(ap, b)
          d = numpy.dot(r, r)

          if min is None or d < min:
            min = d
            min_p = ap

    return (math.sqrt(min), min_p)

  def all_images_within(self, a, b, r):
    """
      Give all images of a to b within distance r.
    """

    images = []

    for i in utils.insideout():
      any_j = False

      for j in utils.insideout():
        any_k = False

        for k in utils.insideout():
          R = numpy.dot(self.lattice.T, numpy.array([float(i), float(j), float(k)]))
          
          if numpy.dot(a+R-b,a+R-b) > r*r:
            if k < 0:
              break
            else:
              continue

          any_k = True
          any_j = True
      
          ap = numpy.add(a, R)
          dr = numpy.subtract(ap, b)
          d = numpy.dot(dr, dr)

          images.append((math.sqrt(d), ap))

        if not any_k and j < 0:
          break

      if not any_j and i < 0:
        break

    images = sorted(images, key=lambda (d,p): d)

    return images

  def __getitem__(self, idx):
    if type(idx) == tuple:
      s, i = idx
      return self.get_species(s, i)
    else:
      return self.atoms[idx]

  def __iter__(self):
    return self.atoms.__iter__()

  def __len__(self):
    return len(self.atoms)

class MagresAtoms(MagresAtomsView):
  """
    A collection of atoms, including lattice parameters, and (if available) lists of NMR parameters.
  """

  def __init__(self, atoms=None, lattice=None):
    if atoms is not None:
      # We've been passed a list of MagresAtom or MagresAtomImage, just drop it in.
      if type(atoms) == list and all([type(atom) in [MagresAtom, MagresAtomImage] for atom in atoms]):
        pass
      # We've been passed a MagresFile, build all the atoms off it
      elif type(atoms) == MagresFile:
        atoms, lattice = self._from_magres(atoms) 
      # We've been passed a list of MagresFiles, merge them and build the atoms
      elif type(atoms) == list and all([type(magres_file) is MagresFile for magres_file in atoms]):
        atoms, lattice = self._from_magres(MagresFile.merge(atoms))
    else:
      atoms = []

    super(MagresAtoms, self).__init__(atoms, lattice)

  def _from_magres(self, magres_file):
    """
      Take a MagresFile and create the MagresAtoms structure.
    """

    self.magres_file = magres_file
    atoms = []

    if 'atoms' in magres_file.data_dict and 'lattice' in magres_file.data_dict['atoms'] and len(magres_file.data_dict['atoms']['lattice']) == 1:
      lattice = numpy.array(magres_file.data_dict['atoms']['lattice'][0])

    if 'atoms' in magres_file.data_dict and 'atom' in magres_file.data_dict['atoms']:
      temp_label_index = {}
      for magres_atom in magres_file.data_dict['atoms']['atom']:
        atom = MagresAtom(magres_atom)

        temp_label_index[(atom.label, atom.index)] = atom

        atoms.append(atom)

    if 'magres' in magres_file.data_dict:
      for tag in magres_file.data_dict['magres']:
        if not (tag.startswith("ms_") or tag == "ms"):
          continue

        ms_type = tag
        
        setattr(self, ms_type, [])

        for magres_ms in magres_file.data_dict['magres'][ms_type]:
          atom = temp_label_index[(magres_ms['atom']['label'], magres_ms['atom']['index'])]
          magres_atom_ms = MagresAtomMs(atom, magres_ms)
          getattr(self, ms_type).append(magres_atom_ms)
          setattr(atom, ms_type, magres_atom_ms)

      for tag in magres_file.data_dict['magres']:
        if not (tag.startswith("efg_") or tag == "efg"):
          continue

        efg_type = tag
        
        setattr(self, efg_type, [])

        for magres_efg in magres_file.data_dict['magres'][efg_type]:
          atom = temp_label_index[(magres_efg['atom']['label'], magres_efg['atom']['index'])]
          magres_atom_efg = MagresAtomEfg(atom, magres_efg)
          getattr(self, efg_type).append(magres_atom_efg)
          setattr(atom, efg_type, magres_atom_efg)
      
      for tag in magres_file.data_dict['magres']:
        if not (tag.startswith("isc_") or tag == "isc"):
          continue

        isc_type = tag

        setattr(self, isc_type, [])

        for magres_isc in magres_file.data_dict['magres'][isc_type]:
          atom1 = temp_label_index[(magres_isc['atom1']['label'], magres_isc['atom1']['index'])]
          atom2 = temp_label_index[(magres_isc['atom2']['label'], magres_isc['atom2']['index'])]
          magres_atom_isc = MagresAtomIsc(atom1, atom2, magres_isc)
          getattr(self, isc_type).append(magres_atom_isc) 

          if not hasattr(atom1, isc_type):
            setattr(atom1, isc_type, {})

          getattr(atom1, isc_type)[atom2] = magres_atom_isc

    return (atoms, lattice)

  @classmethod
  def load_magres(self, f):
    """
      A class method to easily load a :py:class:`magres.format.MagresFile` and return the corresponding MagresAtoms.

      >>> MagresAtoms.load_magres("path/to/magres/file.magres")

      or

      >>> MagresAtoms.load_magres(open("path/to/magres/file.magres"))
    """

    if type(f) == str:
      magres_file = MagresFile(open(f))
    elif type(f) == MagresFile:
      magres_file = f
    elif type(f) == list:
      magres_file = f
    else:
      magres_file = MagresFile(f)
    
    return MagresAtoms(magres_file)

