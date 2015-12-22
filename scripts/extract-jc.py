#!python

"""
  Extract J-couplings from given .magres file. You can specify the perturbing nuclei and receiving nuclei.
"""
from __future__ import print_function
import os, os.path
import sys
import argparse
from magres.atoms import MagresAtoms
from magres.utils import load_all_magres, get_numeric, parse_atom_list

parser = argparse.ArgumentParser(description='Extract J-coupling parameters from one set of atoms to another.')
parser.add_argument('-J', '--J_tensor', action="store_const", help="Display J tensor.", default=False, const=True)
parser.add_argument('-S', '--sort', action="store_const", help="Sort by coupling strength before displaying.",
                    default=False, const=True)
parser.add_argument('-N', '--numbers', action="store_const",
                    help="Parse numbers from path and print. This is useful for e.g. convergence calculations.",
                    default=False, const=True)
parser.add_argument('source_dir', help='Directory to look for calculations in')

parser.add_argument('atoms1', nargs='?', type=str, default=None,
                    help='Only print couplings from these atoms. Specify with atom list notation, e.g. "H1", "H1,H2,'
                         'H3", "H,C", or "H1-3".')
parser.add_argument('atoms2', nargs='?', type=str, default=None,
                    help='Only print couplings to these atoms. Specify with atom list notation, e.g. "H1", "H1,H2,'
                         'H3", "H,C", or "H1-3".')

a = parser.parse_args(sys.argv[1:])

atoms1 = a.atoms1
atoms2 = a.atoms2

if atoms1:
    atoms1_filter = parse_atom_list(atoms1)
else:
    atoms1_filter = lambda x: True

if atoms2:
    atoms2_filter = parse_atom_list(atoms2)
else:
    atoms2_filter = lambda x: True

all_Js = {}

tensors = ['isc', 'isc_fc', 'isc_spin', 'isc_orbital_p', 'isc_orbital_d']

if a.J_tensor:
    print("# Showing in Hz (J)")
    property = "J_iso"
else:
    print("# Showing in 10e19.T^2.J^-1 (K)")
    property = "K_iso"

print("# Number\tAtom1\tAtom2\t{}\tDist\tPath".format("\t".join(tensors)))

lines = []

if os.path.isfile(a.source_dir):
    magres_atoms = [MagresAtoms.load_magres(a.source_dir)]
else:
    magres_atoms = load_all_magres(a.source_dir)

for i, atoms in enumerate(magres_atoms):
    num = get_numeric(atoms.magres_file.path)

    if num:
        idx = num
    else:
        idx = [i]

    for isc in atoms.isc.perturbing(atoms1_filter).receiving(atoms2_filter):
        atom1 = isc.atom1
        atom2 = isc.atom2

        sort_val = abs(getattr(isc, property))

        all_tensors = [getattr(atoms, tensor).perturbing(atom1).receiving(atom2)[0] for tensor in tensors]

        tensor_strs = ["{:.3f}".format(getattr(isc_, property)) for isc_ in all_tensors]

        dist, _ = atoms.least_mirror(atom2.position, atom1.position)

        lines.append((idx,
                      atoms.magres_file.path,
                      str(atom1),
                      str(atom2),
                      tensor_strs,
                      "{:.3F}".format(dist),
                      sort_val))

if a.sort:
    lines = sorted(lines, key=lambda xs: xs[-1])
else:
    lines = sorted(lines, key=lambda xs: xs[0])

for idx, path, atom1, atom2, data, dist, _ in lines:
    if a.numbers:
        print(" ".join(map(str, idx)), atom1, atom2, "\t".join(data), dist, path)
    else:
        print(atom1, atom2, "\t".join(data), dist, path)
