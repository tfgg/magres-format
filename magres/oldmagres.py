# -*- coding: utf-8 -*-
import sys
import re
import math
import numpy
import format

def castep_get_lattice(castep_file):
#        Real Lattice(A)                      Reciprocal Lattice(1/A)
#   15.3790010   0.0000000   0.0000000        0.4085561   0.0000000   0.0000000
#   0.0000000  15.3790010   0.0000000        0.0000000   0.4085561   0.0000000
#   0.0000000   0.0000000  15.1986010        0.0000000   0.0000000   0.4134055

   lattice_regex = re.compile(r"\s+Real Lattice\(A\)\s+Reciprocal Lattice\(1\/A\)\n(.*?)\n(.*?)\n(.*?)\n")

   lattice = lattice_regex.findall(castep_file)[-1] # Get the last lattice in the output

   xx, xy, xz, _, _, _ = map(float, lattice[0].split())
   yx, yy, yz, _, _, _ = map(float, lattice[1].split())
   zx, zy, zz, _, _, _ = map(float, lattice[2].split())

   return [[xx, xy, xz],
           [yx, yy, yz],
           [zx, zy, zz],]

class OldMagres(object):
  def __init__(self, magres_file=None, castep_file=None):
    if magres_file is not None:
      self.parse(magres_file)

    if castep_file is not None:
      self.parse_castep(castep_file)

  def parse_castep(self, castep_file):
    lattice = castep_get_lattice(castep_file)

    if 'lattice' not in self.data['atoms']:
      self.data['atoms']['lattice'] = []
    
    if 'units' not in self.data['atoms']:
      self.data['atoms']['units'] = []

    self.data['atoms']['units'].append(('lattice', 'Angstrom'))
    self.data['atoms']['lattice'].append(lattice)

  def parse(self, magres_file):
    """
      Parse an CASTEP old-style .magres file for total tensors.
    """

    atom_regex = re.compile("[=]+[\r\n]+( Perturbing Atom|Atom): ([A-Za-z\:0-9]+)\s+([0-9]+)[\r\n]+[=]+[\r\n]+([^=]+)[\r\n]+", re.M | re.S)
    shielding_tensor_regex = re.compile("\s{0,}(.*?) Shielding Tensor[\r\n]+\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s+([0-9\.\-]+)[\n\r]+\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s+([0-9\.\-]+)[\n\r]+\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s+")

    jc_tensor_regex = re.compile("\s{0,}J-coupling (.*?)[\r\n]+\s+([0-9eE\.\-]+)\s+([0-9eE\.\-]+)\s+([0-9eE\.\-]+)[\n\r]+\s+([0-9eE\.\-]+)\s+([0-9eE\.\-]+)\s+([0-9eE\.\-]+)[\r\n]+\s+([0-9eE\.\-]+)\s+([0-9eE\.\-]+)\s+([0-9eE\.\-]+)\s+")

    efg_tensor_regex = re.compile("\s{0,}(.*?) tensor[\r\n]+\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s+([0-9\.\-]+)[\r\n]+\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s+([0-9\.\-]+)[\r\n]+\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s+")

    coords_regex = re.compile("([A-Za-z\:0-9]+)\s+([0-9]+)\s+Coordinates\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s+A[\r\n]+")

    atoms = atom_regex.findall(magres_file)

    coords = coords_regex.findall(magres_file)

    self.data = {'atoms': {'atom': [], 'units': []}, 'magres': {'units': []}}

    self.data['atoms']['units'].append(('atom', 'Angstrom'))

    found_atoms = set()

    for s, i, x, y, z in coords:
      i = int(i)
      index = (s,i)

      if index not in found_atoms:
        self.data['atoms']['atom'].append((s, s, i, [[x,y,z]]))
        found_atoms.add(index)

    perturbing_index = ('Al', 15)

    def shape_tensor(els):
      return numpy.reshape(numpy.array(map(float,els)), (3,3)).tolist()

    for atom in atoms:
      index = atom[1].split(":")[0], int(atom[2])
      if atom[0] == " Perturbing Atom":
        perturbing_index = index

    ms_units = False
    efg_units = False
    jc_units = False


    for atom in atoms:
      index = atom[1], int(atom[2])

      shielding_tensors = shielding_tensor_regex.findall(atom[3])
      if len(shielding_tensors) != 0:
        if not ms_units:
          self.data['magres']['units'].append(('ms', 'ppm'))
          ms_units = True

        if 'ms' not in self.data['magres']:
          self.data['magres']['ms'] = []

        for tensor in shielding_tensors:
          self.data['magres']['ms'].append((index[0], index[1], shape_tensor(tensor[1:])))

      efg_tensors = efg_tensor_regex.findall(atom[3])

      if len(efg_tensors) != 0:
        if not efg_units:
          self.data['magres']['units'].append(('efg', 'au'))
          efg_units = True

        if 'efg' not in self.data['magres']:
          self.data['magres']['efg'] = []

        for tensor in efg_tensors:
          self.data['magres']['efg'].append((index[0], index[1], shape_tensor(tensor[1:])))

      jc_tensors = jc_tensor_regex.findall(atom[3])
      if len(jc_tensors) != 0:
        if not jc_units:
          self.data['magres']['units'].append(('isc', '10^19.T^2.J^-1'))
          self.data['magres']['units'].append(('isc_fc', '10^19.T^2.J^-1'))
          self.data['magres']['units'].append(('isc_spin', '10^19.T^2.J^-1'))
          self.data['magres']['units'].append(('isc_orbital_p', '10^19.T^2.J^-1'))
          self.data['magres']['units'].append(('isc_orbital_d', '10^19.T^2.J^-1'))
          jc_units = True

        for tensor in jc_tensors:
          if tensor[0] == "Fermi Contact": tag = "isc_fc"
          if tensor[0] == "Spin Dipole": tag = "isc_spin"
          if tensor[0] == "Diamagnetic": tag = "isc_orbital_d"
          if tensor[0] == "Paramagnetic": tag = "isc_orbital_p"
          if tensor[0] == "Total": tag = "isc"

          if tag not in self.data['magres']:
            self.data['magres'][tag] = []

          self.data['magres'][tag].append((perturbing_index[0], perturbing_index[1], index[0], index[1],  shape_tensor(tensor[1:])))
 
  def as_new_format(self):
    magres_file = format.MagresFile()
    magres_file.data_dict = self.data

    return magres_file

