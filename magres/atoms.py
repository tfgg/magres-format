import constants
import numpy
from format import MagresFile

class MagresAtomEfg(object):
  def __init__(self, atom, magres_efg):
    self.atom = atom
    self.magres_efg = magres_efg

  @property
  def V(self):
    return self.magres_efg['V']

  @property
  def Cq(self):
    try:
      return constants.efg_to_Cq(self.V, self.atom.species)
    except KeyError:
      return None

class MagresAtomIsc(object):
  def __init__(self, atom1, atom2, magres_isc):
    self.atom1 = atom1
    self.atom2 = atom2
    self.magres_isc = magres_isc

  @property
  def K(self):
    return self.magres_isc['K']
  
  @property
  def K_iso(self):
    return numpy.trace(self.K)/3.0

  @property
  def J(self):
    return constants.K_to_J(self.K, self.atom1.species, self.atom2.species)

  @property
  def J_iso(self):
    return numpy.trace(self.J)/3.0

class AtomMs(object):
  def __init__(self, atom, ms, reference=None):
    self.atom = atom
    self.ms = ms
    self.reference = reference

  @property
  def iso(self):
    if self.reference is None:
      reference = 0.0
    else:
      reference = self.reference

    return numpy.linalg.trace(self.ms)/3.0

class MagresAtom(object):
  def __init__(self, magres_atom, efg=None, isc=None, ms=None):
    self.magres_atom = magres_atom

    if efg is not None:
      self.efg = AtomEfg(self, efg)
    
    if isc is not None:
      self.isc = AtomIsc(self, isc)
    
    if ms is not None:
      self.ms = AtomMs(self, ms)

  @property
  def label(self):
    return self.magres_atom['label']
  
  @property
  def species(self):
    return self.magres_atom['species']
  
  @property
  def index(self):
    return self.magres_atom['index']
  
  @property
  def position(self):
    return numpy.array(self.magres_atom['position'])

class MagresAtoms(object):
  def __init__(self, atoms=None):
    if atoms is not None:
      if type(atoms) == list:
        self.atoms = atoms
      elif type(atoms) == MagresFile:
        self.from_magres(atoms) 
    else:
      self.atoms = []

    self.label_index = {}
    self.species_index = {}

    self.build_index()

  def from_magres(self, magres_file):
    self.magres_file = magres_file
    atoms = []

    temp_label_index = {}
    for magres_atom in magres_file.data_dict['atoms']['atom']:
      atom = MagresAtom(magres_atom)

      temp_label_index[(atom.label, atom.index)] = atom

      atoms.append(atom)

    efg_types = ['efg', 'efg_local', 'efg_nonlocal']
    for efg_type in efg_types:
      if efg_type not in magres_file.data_dict['magres']:
        continue
      
      setattr(self, efg_type, [])

      for magres_efg in magres_file.data_dict['magres'][efg_type]:
        atom = temp_label_index[(magres_efg['atom']['label'], magres_efg['atom']['index'])]
        magres_atom_efg = MagresAtomEfg(atom, magres_efg)
        getattr(self, efg_type).append(magres_atom_efg)
        setattr(atom, efg_type, magres_atom_efg)
    
    isc_types = ['isc', 'isc_spin', 'isc_fc', 'isc_orbital_p', 'isc_orbital_d']
    for isc_type in isc_types:
      if isc_type not in magres_file.data_dict['magres']:
        continue

      setattr(self, isc_type, [])

      for magres_isc in magres_file.data_dict['magres'][isc_type]:
        atom1 = temp_label_index[(magres_isc['atom1']['label'], magres_isc['atom1']['index'])]
        atom2 = temp_label_index[(magres_isc['atom2']['label'], magres_isc['atom2']['index'])]
        magres_atom_isc = MagresAtomIsc(atom1, atom2, magres_isc)
        getattr(self, isc_type).append(magres_atom_isc) 

    self.atoms = atoms

  @classmethod
  def load_magres(self, f):
    magres_file = MagresFile(f)
    return MagresAtoms(magres_file)

  def add(self, atoms):
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
    if index is None:
      return self.label_index[label]
    else:
      return self.label_index[label][index]

  def get_species(self, label, index=None):
    if index is None:
      return self.species_index[label]
    else:
      return self.species_index[label][index]

  def get_efg_isotope(self):
    if hasattr(self, 'efg_isotope'):
      return self.efg_isotope
    else:
      return constants.Q_common[self.species]
  
  def get_isc_isotope(self):
    if hasattr(self, 'isc_isotope'):
      return self.isc_isotope
    else:
      return constants.gamma_common[self.species]

  def __getitem__(self, idx):
    if type(idx) == tuple:
      s, i = idx
      return self.get_species(s, i)
    else:
      return self.atoms[idx]

  def __iter__(self):
    return self.atoms.__iter__()
