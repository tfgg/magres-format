import re
import json
import magres.schema.validate

try:
  import numpy
except:
  numpy = None

import sys

blocks_re = re.compile(r"[\[<](?P<block_name>.*?)[>\]](.*?)[<\[]/(?P=block_name)[\]>]", re.M | re.S)

if numpy is not None:
  def tensor33(x):
    return numpy.squeeze(numpy.reshape(x, (3,3))).tolist()
  def tensor31(x):
    return numpy.squeeze(numpy.reshape(x, (3,1))).tolist()
else:
  def tensor33(x):
    return [x[0:3], x[3:6], x[6:]]

  def tensor31(x):
    return x

class BadVersion(Exception):
  pass

class BadMagresFile(Exception):
  pass

class BadUnits(Exception):
  pass

def get_version(file_contents):
  """
    Look for and parse the magres file format version line
  """

  lines = file_contents.split('\n')
  match = re.match("\#\$magres-abinitio-v([0-9]+).([0-9]+)", lines[0])

  if match:
    version = match.groups()
    version = tuple(map(int, version))
  else:
    version = None

  return version

def parse_blocks(file_contents):
  """
    Parse series of XML-like deliminated blocks into a list of (block_name, contents) tuples
  """

  blocks = blocks_re.findall(file_contents)

  return blocks

def parse_block(block):
  """
    Parse block contents into a series of (tag, data) records
  """

  def clean_line(line):
    # Remove comments and whitespace at start and ends of line
    line = re.sub('#(.*?)\n', '', line)
    line = line.strip()

    return line

  name, data = block

  lines = [clean_line(line) for line in data.split('\n')]

  records = []

  for line in lines:
    xs = line.split()

    if not xs:
      continue

    tag = xs[0]
    data = xs[1:]

    records.append((tag, data))

  return (name, records)

def check_units(d):
  """
    Verify that given units for a particular tag are correct.
  """

  allowed_units = {'lattice': 'Angstrom',
                   'atom': 'Angstrom',
                   'ms': 'ppm',
                   'efg': 'au',
                   'efg_local': 'au',
                   'efg_nonlocal': 'au',
                   'isc': '10^19.T^2.J^-1',
                   'isc_fc': '10^19.T^2.J^-1',
                   'isc_orbital_p': '10^19.T^2.J^-1',
                   'isc_orbital_d': '10^19.T^2.J^-1',
                   'isc_spin': '10^19.T^2.J^-1',
                   'isc': '10^19.T^2.J^-1',
                   'sus': '10^-6.cm^3.mol^-1',
                   'calc_cutoffenergy': 'Hartree',}

  if d[0] in d and d[1] == allowed_units[d[0]]:
    pass
  else:
    raise BadUnits("Unrecognized units: %s %s" % (d[0], d[1]))

  return d

def parse_magres_block(block):
  """
    Parse magres block into data dictionary given list of record tuples.
  """

  name, records = block

  # Atom label, atom index and 3x3 tensor
  def sitensor33(name):
     return lambda d: {'atom': {'label': data[0], 'index': int(data[1])}, name: tensor33(list(map(float, data[2:])))}

  # 2x(Atom label, atom index) and 3x3 tensor
  def sisitensor33(name):
     return lambda d: {'atom1': {'label': data[0], 'index': int(data[1])}, 'atom2': {'label': data[2], 'index': int(data[3])}, name: tensor33(list(map(float, data[4:])))}

  tags = {'ms': sitensor33('sigma'),
          'efg': sitensor33('V'),
          'efg_local': sitensor33('V'),
          'efg_nonlocal': sitensor33('V'),
          'isc': sisitensor33('K'),
          'isc_fc': sisitensor33('K'), 'isc_spin': sisitensor33('K'), 'isc_orbital_p': sisitensor33('K'), 'isc_orbital_d': sisitensor33('K'),
          'units': check_units}

  data_dict = {}

  for record in records:
    tag, data = record

    if tag not in data_dict:
      data_dict[tag] = []

    data_dict[tag].append(tags[tag](data))

  return data_dict

def write_units(data, out):
  if 'units' in data:
    for tag, units in data['units']:
      out.append("  units %s %s" % (tag, units))

def tensor_string(tensor):
  return " ".join([" ".join(map(str, xs)) for xs in tensor])

def write_magres_block(data):
  """
    Write out a <magres> block from its dictionary representation
  """

  out = []

  def siout(tag, tensor_name):
    if tag in data:
      for atom_si in data[tag]:
        out.append("  %s %s %d %s" % (tag, atom_si['atom']['label'], atom_si['atom']['index'], tensor_string(atom_si[tensor_name])))

  write_units(data, out)

  siout('ms', 'sigma')

  siout('efg_local', 'V')
  siout('efg_nonlocal', 'V')
  siout('efg', 'V')

  def sisiout(tag, tensor_name):
    if tag in data:
      for isc in data[tag]:
        out.append("  %s %s %d %s %d %s" % (tag, isc['atom1']['label'], isc['atom1']['index'], isc['atom2']['label'], isc['atom2']['index'], tensor_string(isc[tensor_name])))

  sisiout("isc_fc", 'K')
  sisiout("isc_orbital_p", 'K')
  sisiout("isc_orbital_d", 'K')
  sisiout("isc_spin", 'K')
  sisiout("isc", 'K')

  return "\n".join(out)

def parse_atoms_block(block):
  """
    Parse atoms block into data dictionary given list of record tuples.
  """

  name, records = block

  # Lattice record: a1, a2 a3, b1, b2, b3, c1, c2 c3
  def lattice(d):
    return tensor33(list(map(float, data)))

  # Atom record: label, index, x, y, z
  def atom(d):
    return {'species': data[0], 'label': data[1], 'index': int(data[2]), 'position': tensor31(list(map(float, data[3:])))}

  def symmetry(d):
    return " ".join(data)

  tags = {'lattice': lattice,
          'atom': atom,
          'units': check_units,
          'symmetry': symmetry}

  data_dict = {}

  for record in records:
    tag, data = record

    if tag not in data_dict:
      data_dict[tag] = []

    data_dict[tag].append(tags[tag](data))

  return data_dict

def write_atoms_block(data):
  out = []

  write_units(data, out)

  if 'lattice' in data:
    for lat in data['lattice']:
      out.append("  lattice %s" % tensor_string(lat))

  if 'symmetry' in data:
    for sym in data['symmetry']:
      out.append("  symmetry %s" % sym)

  if 'atom' in data:
    for a in data['atom']:
      out.append("  atom %s %s %s %s" % (a['species'], a['label'], a['index'], " ".join(map(str, a['position']))))

  return "\n".join(out)

def parse_generic_block(block):
  """
    Parse any other block into data dictionary given list of record tuples.
  """

  name, records = block

  data_dict = {}

  for record in records:
    tag, data = record

    if tag not in data_dict:
      data_dict[tag] = []

    data_dict[tag].append(data)

  return data_dict

def write_generic_block(data):
  out = []

  for tag, data in list(data.items()):
    for value in data:
      out.append("%s %s" % (tag, " ".join(map(str, value))))

  return "\n".join(out)

class MagresFile(object):
  block_parsers = {'magres': parse_magres_block,
                   'atoms': parse_atoms_block,
                   'calculation': parse_generic_block,}


  block_writers = {'magres': write_magres_block,
                   'atoms': write_atoms_block,
                   'calculation': write_generic_block,}

  version = (1,0)

  def __init__(self, data=None):
    if data is not None:
      self.load(data)

  def load(self, data):
    if type(data) == dict:
      self.data_dict = data
    else:
      if type(data) == str:
        file_contents = data
      else:
        try:
          file_contents = data.read()
        except:
          raise BadMagresFile("Can't load given magres file")

      self.parse(file_contents)

    magres.schema.validate.validate_magres(self.data_dict)

  def parse(self, data, clean=True, include_unrecognised=False):
    if type(data) == str:
      try:
        file_contents = open(data).read()
      except:
        file_contents = data
    elif type(data) == dict:
      self.data_dict = data
      return
    else:
      try:
        file_contents = data.read()
      except:
        raise BadMagresFile("Can't load given magres file")

    version = get_version(file_contents)

    if version is None:
      #pass # Emit a warning?
      raise BadVersion("Version string not present. Possibly not magres file.")
    else:
      if version[0] != self.version[0]: # this is a major version 1 parser
        raise BadVersion("Version %d.%d not recognised. This is a version 1.x parser" % version)

    self.blocks = parse_blocks(file_contents)

    if clean or not hasattr(self, 'data_dict'):
      self.data_dict = {}

    for block_data in self.blocks:
      block = parse_block(block_data)

      if block[0] in self.block_parsers:
        data_dict = self.block_parsers[block[0]](block)
        self.data_dict[block[0]] = data_dict
      else:
        #print "Block type \"%s\" not recognised" % block[0]

        # Throw in the text content of blocks we don't recognise
        if include_unrecognised:
          self.data_dict[block[0]] = block_data[1]

  @classmethod
  def load_json(klass, json_string):
    """
      Load from a json dictionary
    """

    return MagresFile(json.loads(json_string))

  @classmethod
  def merge(klass, magres_files):
    # Make the output magres file using the first magres file's data dictionary
    out_magres_file = MagresFile(dict(magres_files[0].data_dict))
    out_dict = out_magres_file.data_dict

    # We should check that the atoms and lattice are the same here. We don't.

    # Loop through the rest and add their information in
    for magres_file in magres_files[1:]:
      data_dict = magres_file.data_dict

      if 'magres' not in out_dict:
        out_dict['magres'] = {}

      for tag_key in data_dict['magres']:
        if tag_key not in out_dict['magres']:
          out_dict['magres'][tag_key] = data_dict['magres'][tag_key]
        else:
          out_dict['magres'][tag_key] += data_dict['magres'][tag_key]

    return out_magres_file

  def as_json(self):
    """
      Dump as a json dictionary for easy storage/transmission
    """

    magres.schema.validate.validate_magres(self.data_dict)

    json_out = json.dumps(self.data_dict)

    return json_out

  #def as_ase(self):
  #  from ase import Atoms, Atom

  #  if 'lattice' in self.data_dict['atoms']:
  #    atoms = Atoms(cell=self.data_dict['atoms']['lattice'][0], pbc=(1,1,1))
  #  else:
  #    atoms = Atoms()

  #  for s, label, i, pos in self.data_dict['atoms']['atom']:
  #    atom = Atom(s, (pos[0][0], pos[1][0], pos[1][0]), tag=i)
  #    atoms.append(atom)

  #  return atoms

  def __str__(self):
    """
      Convert the internal data dictionary representation to a magres-abinitio format file
    """

    out = []

    out.append("#$magres-abinitio-v%d.%d" % self.version)
    out.append("# Generated by format.py. For format definition and code samples see http://www.ccpnc.ac.uk/pmwiki.php/CCPNC/Fileformat")

    # Order of the blocks, lower weights are higher up. Default weight is 0.
    order = {'calculation': -2,
             'atoms': -1,}

    for block_type in sorted(self.data_dict, key=lambda x: order.get(x,0)):
      data = self.data_dict[block_type]

      textout = self.block_writers[block_type](data)

      out.append("[%s]\n" % block_type + textout + "\n[/%s]" % block_type)

    return "\n".join(out)


if __name__ == "__main__":
  f = open(sys.argv[1])

  magres_file = MagresFile(f)

  # Convert it to json and back again for giggles
  #print magres_file.as_json()
  magres_file.load_json(magres_file.as_json())

  print(magres_file)

