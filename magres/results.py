import yaml
import sys, os
from itertools import product


def load_results(f):
    if type(f) is file:
        data = yaml.load(f)
    elif type(f) is str and os.path.isfile(f):
        data = yaml.load(open(f))
    else:
        data = f

    structures = {}

    for name, structure in list(data['structures'].items()):
        print(name)

        if name not in structures:
            structures[name] = {}

        for coupling in structure['couplings']:
            if type(coupling['index1']) is not list:
                index1 = [coupling['index1']]
            else:
                index1 = coupling['index1']

            if type(coupling['index2']) is not list:
                index2 = [coupling['index2']]
            else:
                index2 = coupling['index2']

            pairs = list(product(index1, index2))

            for pair in pairs:
                structures[name][pair] = coupling

                # Include swizzle
                idx2, idx1 = pair
                structures[name][(idx1, idx2)] = coupling

    return structures


if __name__ == "__main__":
    coupling = structures['TlCl'][('Tl1', 'Cl1')]

    for value in coupling['values']:
        print(value['value'])
