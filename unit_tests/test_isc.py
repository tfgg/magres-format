import math
import numpy
import unittest
import os
from magres.format import MagresFile
from magres.atoms import MagresAtoms

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data")


class MagresIscTest(unittest.TestCase):
    def test_isc(self):
        atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "ethanol-isc.magres"))

        self.assertTrue(hasattr(atoms.C2, "isc"))
        self.assertTrue(hasattr(atoms.C2, "isc_orbital_p"))
        self.assertTrue(hasattr(atoms.C2, "isc_orbital_d"))
        self.assertTrue(hasattr(atoms.C2, "isc_spin"))
        self.assertTrue(hasattr(atoms.C2, "isc_fc"))

        self.assertEqual(len(atoms.C2.isc), len(atoms) - 1)

        perturb_atom = atoms.get('C', 2)

        for isc in atoms.C2.isc:
            self.assertEqual(isc.atom1, perturb_atom)
            self.assertEqual(isc.K_iso, numpy.trace(isc.K) / 3.0)

        # Check that the principal components are ordered by the Haeberlen convention
        for isc in atoms.C2.isc:
            self.assertTrue(abs(isc.K_evals[2] - isc.K_iso) >= abs(isc.K_evals[0] - isc.K_iso))
            self.assertTrue(abs(isc.K_evals[0] - isc.K_iso) >= abs(isc.K_evals[1] - isc.K_iso))

    def test_replace_isc(self):
        atoms = MagresAtoms.load_magres(os.path.join(DATA_DIR, "ethanol-isc.magres"))

        atoms.C2.isc[0].K

        with self.assertRaises(Exception):
            atoms.C2.isc[0].K = numpy.array([1, 2, 3])

        orig_iso = atoms.C2.isc[0].K_iso
        new_K = atoms.C2.isc[0].K * 2.0

        atoms.C2.isc[0].K = new_K

        self.assertEqual(atoms.C2.isc[0].K_iso, orig_iso * 2.0)

    def test_full_isc(self):
        magres_files = [MagresFile(open(os.path.join(DATA_DIR, "ethanol", f))) for f in
                        os.listdir(os.path.join(DATA_DIR, "ethanol")) if f.startswith('ethanol-jc')]

        self.assertEqual(len(magres_files), 9)

        atoms = MagresAtoms.load_magres(magres_files)

        for atom in atoms:
            self.assertTrue(hasattr(atom, "isc"))

        self.assertEqual(len(atoms.isc), len(atoms) * (len(atoms) - 1))

        # Check every atom has couplings to every other atom
        for atom in atoms:
            self.assertEqual(len(atom.isc), len(atoms) - 1)


if __name__ == "__main__":
    unittest.main()
