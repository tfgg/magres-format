#!python
import sys
import yaml

# Name of dataset to compare against
main_dataset = sys.argv[2]
compare_datasets = sys.argv[3:]

data = yaml.load(open(sys.argv[1]))

for name, structure in data['structures'].items():
  print "#" + name

  for coupling in structure['couplings']:
    values = {}

    for value in coupling['values']:
      values[value['dataset']['slug']] = value['value']

    if type(coupling['index1']) is list:
      indices1 = coupling['index1']
    else:
      indices1 = [coupling['index1']]

    if type(coupling['index2']) is list:
      indices2 = coupling['index2']
    else:
      indices2 = [coupling['index2']]

    print "#  " + ",".join(indices1) + " --> " + ",".join(indices2)

    vals = []
    abs_errors = []
    rel_errors = []
    
    compare = values[main_dataset]

    for compare_dataset in compare_datasets:
      try:
        value = values[compare_dataset]
      except:
        print "# Missing coupling for", name, "from", compare_dataset
        value = None

      if value is not None:
        vals.append("%6s" % "%.3f" % value)
        abs_errors.append("%6s" % "%.3f" % abs(abs(value) - abs(compare)))
        rel_errors.append("%6s" % "%.3f" % abs(abs(value) - abs(compare)) / abs(compare))
      else:
        vals.append("?")
        abs_errors.append("?")
        rel_errors.append("?")

    print "    ", "%6s" % "%.3f" % compare, "\t", "\t".join(vals), "\t".join(abs_errors), "\t", "\t".join(rel_errors), name
