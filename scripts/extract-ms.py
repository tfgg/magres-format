#!python
from __future__ import print_function
import re
import os
import sys
import argparse
from magres.atoms import MagresAtoms
from magres.utils import load_all_magres, get_numeric, parse_atom_list

parser = argparse.ArgumentParser(description='Extract magnetic shielding parameters from a set of calculations.')
parser.add_argument('-N', '--numbers', action="store_const",
                    help="Parse numbers from path and print. This is useful for e.g. convergence calculations.",
                    default=False, const=True)
parser.add_argument('source_dir', help='Directory to look for calculations below.')
parser.add_argument('--iso', action='append')
parser.add_argument('atoms', nargs='?', type=str, default=None,
                    help='Which atoms to print shieldings of. Specify with atom list notation, e.g. "H1" or "H1,H2,'
                         'H3" or "H,C" or "H1-3".')

a = parser.parse_args(sys.argv[1:])

atoms_filter_str = a.atoms

if atoms_filter_str:
    atoms_filter = parse_atom_list(atoms_filter_str)
else:
    atoms_filter = lambda x: True

tensors = ['ms']

lines = []

print("# Number\tAtom\tIso\tAniso\tAsym\tPath")

if os.path.isfile(a.source_dir):
    magres_atoms = [MagresAtoms.load_magres(a.source_dir)]
else:
    magres_atoms = load_all_magres(a.source_dir)

isos = {}

for iso_ in a.iso:
    num = re.findall('([0-9]+)', iso_)[0]
    sym = re.findall('([A-Za-z]+)', iso_)[0]

    isos[sym] = int(num)

for i, atoms in enumerate(magres_atoms):
    num = get_numeric(atoms.magres_file.path)

    for s, iso in list(isos.items()):
        for atom in atoms.species(s):
            atom.isotope = iso

    if num:
        idx = num
    else:
        idx = [i]

    for atom in atoms:
        if atoms_filter(atom) and \
                hasattr(atom, 'ms'):

            lines.append((idx,
                          atoms.magres_file.path,
                          str(atom),
                          [atom.ms.iso, atom.ms.aniso, atom.ms.eta]))

lines = sorted(lines, key=lambda xs: xs[0])

for idx, path, atom, data in lines:
    if a.numbers:
        print(" ".join(map(str, idx)) + "\t" + atom + "\t" + "\t".join("{:.3f}".format(x) for x in data) + "\t" + path)
    else:
        print(atom + "\t" + "\t".join("{:.3f}".format(x) for x in data) + "\t" + path)
