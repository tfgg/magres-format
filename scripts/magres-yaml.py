#!python
import yaml
import sys, os
from numpy import mean, std
from magres.utils import find_all
from magres.atoms import MagresAtoms
from castepy.utils import calc_from_path

magres_files = []

for dir in sys.argv[3:]:
  magres_files += find_all(dir, '.magres')

data = yaml.load(open(sys.argv[1]))
dataset_name = yaml.load(sys.argv[2])

couplings_map = {}
to_mean = {}

for magres_file in magres_files:
  print >>sys.stderr, magres_file
  try:
    atoms = MagresAtoms.load_magres(magres_file)
  except:
    continue

  dir, name = calc_from_path(magres_file)

  if name in data['structures'] and hasattr(atoms, 'isc'):
    couplings = data['structures'][name]['couplings']

    for isc in atoms.isc:
      idx1 = "%s%d" % (isc.atom1.species, isc.atom1.index)
      idx2 = "%s%d" % (isc.atom2.species, isc.atom2.index)

      if idx1 != idx2:
        for i, coupling in enumerate(couplings):
          if (idx1 in coupling['index1'] and idx2 in coupling['index2']) or \
             (idx1 in coupling['index2'] and idx2 in coupling['index1']):
  
            value = float(getattr(isc, coupling['expr']))

            coupling_id = "%s-%d" % (name, i)

            if coupling_id not in to_mean:
              couplings_map[coupling_id] = coupling
              to_mean[coupling_id] = []

            to_mean[coupling_id].append(value)

if dataset_name in data['datasets']:
  dataset = data['datasets'][dataset_name]
else:
  dataset = {'name': 'Automatically extracted',
             'dirs': sys.argv[3:],
             'cwd': os.getcwd(),
             'slug': dataset_name}

data['datasets'][dataset_name] = dataset

for coupling_id in to_mean:
  values = to_mean[coupling_id]
  coupling = couplings_map[coupling_id]

  value_mean = float(mean(values))
  value_stdev = float(std(values))
  value_count = len(values)

  if 'values' not in coupling:
    coupling['values'] = []

  coupling['values'].append({'value': value_mean,
                             'stdev': value_stdev,
                             'num': value_count,
                             'dataset': dataset})

print yaml.dump(data, default_flow_style=False)
