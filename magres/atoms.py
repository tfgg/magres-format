"""
  magres.atoms is a collection of user-friendly data structures for representing
  groups of atoms and their NMR parameters. The NMR parameter structures can
  calculate a variety of properties from the underlying tensors, according to
  a variety of conventions.
"""
from __future__ import print_function

VERSION = (0, 0)

import sys, os
import math
import numpy
import re
from . import constants
from . import html_repr
from numpy import (
    cross,
    dot,
    sqrt,
    rad2deg,
    arctan2
)
from numpy.linalg import norm
from .format import MagresFile
from .efg import MagresAtomEfg
from .isc import MagresAtomIsc
from .ms import MagresAtomMs
from .atom import MagresAtom, MagresAtomImage
from .view import ListPropertyView

element_colours = {
    'H': ("#EEEEEE", "#000000"),
    'C': ("#999999", "#000000"),
    'O': ("#FF0000", "#FFFFFF"),
    'N': ("#0000FF", "#FFFFFF"), }

min_dist = 1.0
max_dist = 2.0


def insideout():
    """
      Count up in positive numbers and down in negative numbers
    """

    yield 0

    i = 1

    while True:
        yield i
        yield -i
        i += 1


class LabelNotFound(Exception):
    pass


class SpeciesNotFound(Exception):
    pass


class AtomNotFound(Exception):
    pass


class MagresAtomsView(object):
    """
      A container for a collection of atoms with an optional lattice.
    """

    __slots__ = ["atoms", "lattice", "label_index", "species_index"]

    def __init__(self, atoms=None, lattice=None):
        if atoms is not None:
            self.atoms = atoms
        else:
            self.atoms = []

        if lattice is not None:
            self.lattice = lattice
        else:
            self.lattice = None

        self.label_index = {}
        self.species_index = {}

        self._build_index()

    def __repr__(self):
        return "<magres.atom.MagresAtomsView - {} atoms>".format(len(self))

    def add(self, atoms):
        """
          Add an extra atom or list of atoms.
        """

        if type(atoms) == Atom:
            self.atoms.append(atoms)
        elif type(atom) == list:
            for atom in atoms:
                self.atoms.append(atom)

        self._build_index()

    def _build_index(self):
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

    def filter(self, fn):
        """
          Filter atoms by some function fn.
        """
        return MagresAtomsView([x for x in self if fn(x)], self.lattice)

    def get(self, species, index):
        """
          Get a single atom of a particular species and index.

          >>> atoms.get('C', 2)
        """

        if species not in self.species_index:
            return SpeciesNotFound("Species {} does not exist in this system".format(species))
        else:
            if len(self.species_index[species]) >= index:
                return self.species_index[species][index - 1]
            else:
                raise AtomNotFound("Atom %s %d does not exist in this system. There are %d atoms in the %s species." % (
                species, index, len(self.species_index[species]), species))

    def species(self, *args):
        """
          Return a MagresAtomsView containing only atoms of the specified species.

          >>> atoms.species('C')
        """

        species = args

        rtn_atoms = []
        for s in species:
            if s in self.species_index:
                rtn_atoms += self.species_index[s]

        return MagresAtomsView(rtn_atoms, self.lattice)

    def set_reference(self, reference):
        for atom in self.atoms:
            atom.reference = float(reference)

    def within(self, pos, max_dr):
        """
          Return all atoms within max_dr Angstroms of pos, including all images.

          >>> atoms.within(p, 5.0)
        """

        if type(pos) is MagresAtom or type(pos) is MagresAtomImage:
            pos = pos.position

        atoms = []

        for atom in self.atoms:
            if type(atom) == MagresAtom:
                images = self._all_images_within(atom.position, pos, max_dr)

                for image_dist, image_pos in images:
                    if image_dist <= max_dr:
                        atoms.append(MagresAtomImage(atom, image_pos))
            elif type(atom) == MagresAtomImage:
                if atom.dist(pos) <= max_dr:
                    atoms.append(atom)

        return MagresAtomsView(atoms, self.lattice)

    def least_mirror(self, a, b):
        """
          Give the closest periodic image of a to b given the current lattice.
        """
        min = None
        min_p = None

        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    ap = numpy.add(a, numpy.dot(self.lattice.T, (float(i), float(j), float(k))))
                    r = numpy.subtract(ap, b)
                    d = numpy.dot(r, r)

                    if min is None or d < min:
                        min = d
                        min_p = ap

        return (math.sqrt(min), min_p)

    def dist(self, atom1, atom2):
        return self.least_mirror(atom1.position, atom2.position)[0]

    def angle(self, atom1, atom2, atom3, degrees=False):
        dr1 = self.least_mirror(atom1.position, atom2.position)[1] - atom2.position
        dr2 = self.least_mirror(atom3.position, atom2.position)[1] - atom2.position

        dr1 = dr1 / sqrt(dot(dr1, dr1))
        dr2 = dr2 / sqrt(dot(dr2, dr2))

        if degrees:
            return math.acos(dot(dr1, dr2)) * 180.0 / math.pi
        else:
            return math.acos(dot(dr1, dr2))

    def dihedral(self, atom1, atom2, atom3, atom4, degrees=False):
        # Use least_mirror for positions rather than atom.position to allow for a closer image
        dr12 = self.least_mirror(atom2.position, atom1.position)[1] - atom1.position
        dr23 = self.least_mirror(atom3.position, atom2.position)[1] - atom2.position
        dr34 = self.least_mirror(atom4.position, atom3.position)[1] - atom3.position

        # Make unit vectors
        dr12 = dr12 / norm(dr12)
        dr23 = dr23 / norm(dr23)
        dr34 = dr34 / norm(dr34)

        norm1 = cross(dr12, dr23)
        norm2 = cross(dr23, dr34)

        # Get m, where m, norm1 and bond2 form an orthogonal frame
        m = cross(dr23, norm1)

        # x = cos(theta) where theta is the angle between norm1 and norm2
        x = dot(norm1, norm2)

        # y = cos(phi) where phi is the angle between m and norm2
        y = dot(m, norm2)

        if degrees:
            return rad2deg(arctan2(y, x))
        else:
            return arctan2(y, x)

    def _all_images_within(self, a, b, r):
        """
          Give all images of a to b within distance r.
        """

        images = []

        for i in insideout():
            any_j = False

            for j in insideout():
                any_k = False

                for k in insideout():
                    R = numpy.dot(self.lattice.T, numpy.array([float(i), float(j), float(k)]))

                    if numpy.dot(a + R - b, a + R - b) > r * r:
                        if k < 0:
                            break
                        else:
                            continue

                    any_k = True
                    any_j = True

                    ap = numpy.add(a, R)
                    dr = numpy.subtract(ap, b)
                    d = numpy.dot(dr, dr)

                    images.append((math.sqrt(d), ap))

                if not any_k and j < 0:
                    break

            if not any_j and i < 0:
                break

        images = sorted(images, key=lambda d_p: d_p[0])

        return images

    re_species_index = re.compile('([A-Za-z]+)([0-9]+)')

    def __getattribute__(self, attr_name):
        try:
            return object.__getattribute__(self, attr_name)
        except AttributeError:
            try:
                s, i = self.re_species_index.findall(attr_name)[0]
                return self.get(s, int(i))
            except Exception:
                return getattr(ListPropertyView(self.atoms), attr_name)

    def __getitem__(self, idx):
        if type(idx) == tuple:
            s, i = idx
            return self.get(s, i)
        else:
            return self.atoms[idx]

    def __iter__(self):
        return self.atoms.__iter__()

    def __len__(self):
        return len(self.atoms)

    def __add__(self, b):
        if type(b) is MagresAtomsView and (self.lattice == b.lattice).all():
            new_atoms = set(self.atoms).union(set(b.atoms))

            return MagresAtomsView(list(new_atoms), self.lattice)

        elif type(b) is MagresAtom:
            new_atoms = set(self.atoms + [b])

            return MagresAtomsView(list(new_atoms), self.lattice)

    def __radd__(self, b):
        if type(b) is MagresAtom:
            new_atoms = set(self.atoms + [b])

            return MagresAtomsView(list(new_atoms), self.lattice)

    def _repr_png_(self):
        return None  # Disable for now

        try:
            import pydot
        except ImportError:
            return None

        dist_graph = pydot.Dot(graph_type='graph', size="100", prog='neato', dim=2)

        has_isc = numpy.array([hasattr(atom, 'isc') for atom in self.atoms]).any()
        has_ms = numpy.array([hasattr(atom, 'ms') for atom in self.atoms]).any()

        def lm_dist(atom1, atom2):
            return self.least_mirror(atom1.position, atom2.position)

        def strength_color(dist):
            x = min(max((abs(dist) - min_dist) / (max_dist - min_dist), 0.0), 1.0)
            y = (1.0 - x) * 255

            return "#000000{:02x}".format(int(y))

        for atom in self:
            fillcolor, fontcolor = element_colours.get(atom.species, ("#CCCCCC", "#000000"))

            # if has_ms:
            #  label = "{}\n{:.3f}".format(str(atom), atom.ms.iso)
            # else:
            label = str(atom)  # "{}\n{:.3f}".format(str(atom), 3.0)

            node = pydot.Node(str(atom),
                              label=label,
                              shape="circle",
                              style="filled",
                              size="0.01",
                              fillcolor=fillcolor,
                              fontcolor=fontcolor,
                              fontsize=10)

            dist_graph.add_node(node)

        bonds_done = set()
        atoms_set = set(self.atoms)

        for atom1 in self:
            for atom2 in atom1.bonded:
                if atom2 not in atoms_set:
                    continue

                dist, pos = lm_dist(atom1, atom2)

                idx1 = (atom2, atom1)
                idx2 = (atom1, atom2)

                if idx1 not in bonds_done and idx2 not in bonds_done:
                    bonds_done.add(idx2)

                    # Hide bonds if we have ISC, but keep for structure
                    # if has_isc:
                    #  color = "#000000{:02x}".format(64)
                    # else:

                    color = strength_color(dist)

                    edge = pydot.Edge(str(atom1),
                                      str(atom2),
                                      color=color,
                                      len=dist,
                                      fontsize=8)

                    dist_graph.add_edge(edge)

        def flatten(xs):
            rtn = xs[0]
            for x in xs[1:]:
                rtn += x
            return rtn

        if has_isc:
            isc_done = set()

            K_isos = [abs(isc.K_iso) for isc in self.isc
                      if isc.atom1 in atoms_set and isc.atom2 in atoms_set]

            min_isc = min(K_isos)
            max_isc = max(K_isos)

            print(K_isos)

            def strength_color_K(K_iso):
                x = min(max((abs(K_iso) - min_isc) / (max_isc - min_isc), 0.0), 1.0)
                y = x * 255

                if K_iso > 0.0:
                    return "#FF0000{:02X}".format(int(y))
                else:
                    return "#0000FF{:02X}".format(int(y))

            for isc in self.isc:
                print(isc)
                strength = min(max((abs(isc.K_iso) - min_isc) / (max_isc - min_isc), 0.0), 1.0)

                if strength > 0.01 and (isc.atom2, isc.atom1) not in isc_done \
                        and isc.atom1 in atoms_set and isc.atom2 in atoms_set:

                    isc_done.add((isc.atom1, isc.atom2))

                    edge = pydot.Edge(str(isc.atom1),
                                      str(isc.atom2),
                                      fontsize=8,
                                      color=strength_color_K(isc.K_iso),
                                      fontcolor=strength_color_K(isc.K_iso),
                                      label="{:.3f}".format(isc.K_iso),
                                      len=lm_dist(isc.atom1, isc.atom2)[0])

                    dist_graph.add_edge(edge)

        return dist_graph.create_png(prog='neato')


class MagresAtoms(MagresAtomsView):
    """
      A collection of atoms, including lattice parameters, and (if available) lists of NMR parameters.
    """

    def __init__(self, atoms=None, lattice=None):
        if atoms is not None:
            # We've been passed a list of MagresAtom or MagresAtomImage, just drop it in.
            if type(atoms) == list and all([type(atom) in [MagresAtom, MagresAtomImage] for atom in atoms]):
                pass
            # We've been passed a MagresFile, build all the atoms off it
            elif type(atoms) == MagresFile:
                atoms, lattice = self._from_magres(atoms)
            # We've been passed a list of MagresFiles, merge them and build the atoms
            elif type(atoms) == list and all([type(magres_file) is MagresFile for magres_file in atoms]):
                atoms, lattice = self._from_magres(MagresFile.merge(atoms))
        else:
            atoms = []

        super(MagresAtoms, self).__init__(atoms, lattice)

        if len(self) < 20:
            self.calculate_bonds()

    def __repr__(self):
        return "<magres.atom.MagresAtoms - {} atoms>".format(len(self))

    def _from_magres(self, magres_file):
        """
          Take a MagresFile and create the MagresAtoms structure.
        """

        self.magres_file = magres_file
        atoms = []

        if 'atoms' in magres_file.data_dict and 'lattice' in magres_file.data_dict['atoms'] and len(
                magres_file.data_dict['atoms']['lattice']) == 1:
            lattice = numpy.array(magres_file.data_dict['atoms']['lattice'][0])

        if 'atoms' in magres_file.data_dict and 'atom' in magres_file.data_dict['atoms']:
            temp_label_index = {}
            for magres_atom in magres_file.data_dict['atoms']['atom']:
                atom = MagresAtom(magres_atom)

                temp_label_index[(atom.label, atom.index)] = atom

                atoms.append(atom)

        if 'magres' in magres_file.data_dict:
            for tag in magres_file.data_dict['magres']:
                if not (tag.startswith("ms_") or tag == "ms"):
                    continue

                ms_type = tag

                for magres_ms in magres_file.data_dict['magres'][ms_type]:
                    atom = temp_label_index[(magres_ms['atom']['label'], magres_ms['atom']['index'])]
                    magres_atom_ms = MagresAtomMs(atom, magres_ms)
                    # getattr(self, ms_type).append(magres_atom_ms)
                    setattr(atom, ms_type, magres_atom_ms)

            for tag in magres_file.data_dict['magres']:
                if not (tag.startswith("efg_") or tag == "efg"):
                    continue

                efg_type = tag

                for magres_efg in magres_file.data_dict['magres'][efg_type]:
                    atom = temp_label_index[(magres_efg['atom']['label'], magres_efg['atom']['index'])]
                    magres_atom_efg = MagresAtomEfg(atom, magres_efg)
                    # getattr(self, efg_type).append(magres_atom_efg)
                    setattr(atom, efg_type, magres_atom_efg)

            for tag in magres_file.data_dict['magres']:
                if not (tag.startswith("isc_") or tag == "isc"):
                    continue

                isc_type = tag

                # setattr(self, isc_type, IscListPropertyView([]))

                for magres_isc in magres_file.data_dict['magres'][isc_type]:
                    atom1 = temp_label_index[(magres_isc['atom1']['label'], magres_isc['atom1']['index'])]
                    atom2 = temp_label_index[(magres_isc['atom2']['label'], magres_isc['atom2']['index'])]

                    # Don't bother with self couplings
                    if atom1 == atom2:
                        continue

                    magres_atom_isc = MagresAtomIsc(atom1, atom2, magres_isc)

                    # getattr(self, isc_type).append(magres_atom_isc)

                    if not hasattr(atom1, isc_type):
                        setattr(atom1, isc_type, ListPropertyView([]))

                    getattr(atom1, isc_type).append(magres_atom_isc)

        return (atoms, lattice)

    def load_bonds(self, castep_file, pop_tol=0.2):
        from castepy.output.bonds import parse_bonds
        from collections import Counter

        bonds = next(parse_bonds(castep_file))

        bonded_dict = dict([((atom.species, atom.index), []) for atom in self])

        for idx1, idx2, pop, length in bonds:
            if pop >= pop_tol:
                bonded_dict[idx1].append(idx2)
                bonded_dict[idx2].append(idx1)

        for idx1, idx2s in list(bonded_dict.items()):
            atom1 = self.get(*idx1)

            bonded_atoms = []

            for idx2 in idx2s:
                atom2 = self.get(*idx2)
                bonded_atoms.append(atom2)

            atom1.bonded = MagresAtomsView(list(bonded_atoms), self.lattice)

    def calculate_bonds(self, tol=2.0):
        for atom1 in self.atoms:
            bonded_atoms = []

            for atom2 in self.within(atom1, tol):
                if (atom1.position != atom2.position).any():
                    bonded_atoms.append(atom2)

            atom1.bonded = MagresAtomsView(list(bonded_atoms), self.lattice)

    @classmethod
    def load_magres(self, f):
        """
          A class method to easily load a :py:class:`magres.format.MagresFile` and return the corresponding MagresAtoms.

          >>> MagresAtoms.load_magres("path/to/magres/file.magres")

          or

          >>> MagresAtoms.load_magres(open("path/to/magres/file.magres"))
        """

        if type(f) == str:
            magres_file = MagresFile(open(f))
            magres_file.path = f
        elif type(f) == MagresFile:
            magres_file = f
        elif type(f) == list:
            magres_file = f
        else:
            magres_file = MagresFile(f)

        return MagresAtoms(magres_file)
