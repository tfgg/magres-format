#!python
import yaml
import pyaml
import sys, os
from numpy import mean, std, array
from magres.utils import find_all
from magres.atoms import MagresAtoms
from castepy.utils import calc_from_path

magres_files = []

for dir in sys.argv[3:]:
  magres_files += find_all(dir, '.magres')

data = yaml.load(open(sys.argv[1]))
dataset_name = sys.argv[2]

seed_map = {}

# make all the coupling expressions lists
for structure_name, structure in data['structures'].items():
  couplings = structure['couplings']
  for i, coupling in enumerate(couplings):
    if type(coupling['index1']) is not list:
      coupling['index1'] = [x.strip() for x in coupling['index1'].split(',')]

    if type(coupling['index2']) is not list:
      coupling['index2'] = [x.strip() for x in coupling['index2'].split(',')]

    if type(coupling['expr']) is not list:
      coupling['expr'] = [x.strip() for x in coupling['expr'].split(',')]

  if 'seeds' in structure:
    for seed in structure['seeds']:
      seed_map[seed] = structure_name

couplings_map = {}
to_mean = {}

for magres_file in magres_files:
  print >>sys.stderr, magres_file
  try:
    atoms = MagresAtoms.load_magres(magres_file)
  except:
    continue

  dir, name = calc_from_path(magres_file)

  structure = None
  if name in seed_map:
    structure = data['structures'][seed_map[name]]
  elif name in data['structures']:
    structure = data['structures'][name]

  if structure is not None:
    couplings = structure['couplings']

    # Go through each coupling in the structure and check if we want it.
    for atom in atoms:
      if hasattr(atom, 'isc'):
        for isc in atom.isc.values():
          idx1 = "%s%d" % (isc.atom1.species, isc.atom1.index)
          idx2 = "%s%d" % (isc.atom2.species, isc.atom2.index)

          # Don't extract an atom's coupling with itself. That's weird.
          if idx1 != idx2:
            atom1 = isc.atom1
            atom2 = isc.atom2

            print >>sys.stderr, idx1, idx2

            for i, coupling in enumerate(couplings):
              values = []

              # Do we match this coupling?
              if (idx1 in coupling['index1'] and idx2 in coupling['index2']) or \
                 (idx1 in coupling['index2'] and idx2 in coupling['index1']):

                for expr in coupling['expr']:
                  tensor, quantity = expr.strip().split('.')

                  if hasattr(atom1, tensor):
                    atom1_tensor = getattr(atom1, tensor)[atom2]

                    if hasattr(atom1_tensor, quantity):
                      value = float(getattr(atom1_tensor, quantity))
                    else:
                      value = None
                  else:
                    value = None

                  values.append(value)

                coupling_id = "%s-%d" % (name, i)

                if coupling_id not in to_mean:
                  couplings_map[coupling_id] = coupling
                  to_mean[coupling_id] = []

                to_mean[coupling_id].append(values)

if dataset_name in data['datasets']:
  dataset = data['datasets'][dataset_name]
else:
  dataset = {'name': 'Automatically extracted',
             'dirs': sys.argv[3:],
             'cwd': os.getcwd(),
             'slug': dataset_name}

data['datasets'][dataset_name] = dataset

for coupling_id in to_mean:
  values = array(to_mean[coupling_id])
  coupling = couplings_map[coupling_id]

  value_means = map(float, mean(values,axis=0).tolist())
  value_stdevs = map(float, std(values,axis=0).tolist())
  value_count = values.shape[0]

  if 'values' not in coupling:
    coupling['values'] = []

  value_obj = {'value': value_means,
               'dataset': dataset}
               
  if value_count != 1:
    value_obj['stdev'] = value_stdevs
    value_obj['num'] = value_count
 
  coupling['values'].append(value_obj)

print yaml.dump(data, default_flow_style=False)#, sys.stdout, vspacing=[2,1])

