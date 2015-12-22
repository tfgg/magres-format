#!python
from __future__ import print_function
import sys
import yaml
import re
from magres.constants import K_to_J

# Expression and dataset to compare against
expr = sys.argv[2]
main_dataset = sys.argv[3]

# Other datasets
compare_datasets = sys.argv[4:]

data = yaml.load(open(sys.argv[1]))

# make all the coupling expressions lists
for structure_name in data['structures']:
    couplings = data['structures'][structure_name]['couplings']
    for i, coupling in enumerate(couplings):
        if type(coupling['expr']) is not list:
            coupling['expr'] = coupling['expr'].split(',')

tensor, quantity = expr.strip().split('.')

for name, structure in list(data['structures'].items()):
    print("#" + name)

    for coupling in structure['couplings']:
        values = {}

        expr_idx = -1

        if type(coupling['expr']) is list:
            try:
                expr_idx = coupling['expr'].index(expr)
            except:
                expr_idx = None
        elif coupling['expr'] != expr:
            expr_idx = None

        if type(coupling['index1']) is list:
            indices1 = coupling['index1']
        else:
            indices1 = [coupling['index1']]

        if type(coupling['index2']) is list:
            indices2 = coupling['index2']
        else:
            indices2 = [coupling['index2']]


        def s(x):
            return re.findall('(.*?)([0-9]+)', x)[0][0]


        bond_symbol = s(indices1[0]) + "-->" + s(indices2[0])

        print("#  " + bond_symbol)

        convert_K_to_J = False
        # Ok, so we don't have J_iso, but what about K_iso?
        # If we do, we can convert it to J_iso.
        # if quantity == "J_iso" and expr_idx is None:
        #  try:
        #    expr_idx = coupling['expr'].index("%s.%s" % (tensor, "K_iso"))
        #    convert_J_to_K = True
        #  except:
        #    expr_idx = None

        if expr_idx is None:
            print("# %s not present" % expr)
            continue
            # break

        for value in coupling['values']:
            if type(value['value']) is not list:
                vals = list(map(float, str(value['value']).split(',')))
            else:
                vals = value['value']

            if expr_idx >= 0 and expr_idx < len(vals):
                values[value['dataset']['slug']] = vals[expr_idx]
            else:
                values[value['dataset']['slug']] = vals[0]

                # if convert_J_to_K:
                #  values[value['dataset']['slug']] = K_to_J(vals[expr_idx], s1, s2)

        vals = []
        abs_errors = []
        rel_errors = []

        try:
            compare = values[main_dataset]
        except KeyError:
            print("# skipping, no %s value" % main_dataset)
            continue

        for compare_dataset in compare_datasets:
            try:
                value = values[compare_dataset]
            except:
                print("# Missing coupling for", name, "from", compare_dataset)
                value = None

            if value is not None:
                vals.append("%6s" % "%.3f" % value)
                abs_errors.append("%6s" % "%.3f" % (abs(abs(value) - abs(compare))))

                if abs(compare) != 0.0:
                    rel_errors.append("%6s" % "%.3f" % (abs(abs(value) - abs(compare)) / abs(compare)))
                else:
                    rel_errors.append("%6s" % "?")
            else:
                vals.append("?")
                abs_errors.append("%6s" % "?")
                rel_errors.append("%6s" % "?")

        print(
            "    ", "%6s" % "%.3f" % compare, "\t", "\t".join(vals), "\t".join(abs_errors), "\t", "\t".join(rel_errors),
            name, bond_symbol)
