import unittest

import magres.format as format

class FormatTest(unittest.TestCase):
  magres1_path = "test_data/ethanol/ethanol.atoms"

  def test_load(self):
    magres = format.load_magres(open(self.magres1_path).read())

    self.assertIn('lattice', magres, "Lattice not found") 
    self.assertIn('atom', magres, "Atoms not found")
    self.assertIn('ms', magres, "Magnetic shielding not found")

  def test_dump(self):
    c = calc.CastepCalc(self.calc1_path[0], self.calc1_path[1])
    c.load()

    magres_text = format.write_cell(c.cell)

    self.assertIn('lattice', magres_text, "Lattice not dumped")
    self.assertIn('atom', magres_text, "Atoms not dumped")
    self.assertIn('ms', magres_text, "Magnetic shielding not dumped")
    self.assertIn('bond', magres_text, "Bonds not dumped")

    magres = format.load_magres(magres_text)
    ions = format.load_into_ions(magres)

    self.assertEqual(len(ions), len(c.cell.ions), "Number of ions read back in does not match")

    # Bonds not loaded from magres yet
    # self.assertEqual(len(ions.bonds), len(c.cell.ions.bonds), "Number of bonds read back in does not match")

if __name__ == "__main__":
  unittest.main()

