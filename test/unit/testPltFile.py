import unittest
from PltFile import PltFile
import matplotlib.pyplot as plt

class TestPltFile(unittest.TestCase):
    DATA_FILE = 'data/plt000100.2d.hdf5'
    CHK_DATA_FILE = 'data/chk000100.2d.hdf5'

    def test_load(self):

        # Test loading a file that doesn't exist
        pf = PltFile('file that does not exist')
        self.assertEqual(pf.defined, False)

        # Test loading a file that does exist
        pf = PltFile(self.DATA_FILE, load_data=True)
        self.assertEqual(pf.num_levels, 3)
        self.assertEqual(pf.plot_prefix, 'plt')
        self.assertEqual(pf.frame, 100)

        # Test pretty printing object
        print(pf)

        # Test loading after data already loaded (should do nothing)
        pf.load_data()

        # Check data is correct
        pf.load_data()
        max_enthalpy = float(pf.get_level_data('Enthalpy', level=0).max())
        self.assertAlmostEqual(max_enthalpy, 6.302214, 6)
        max_enthalpy = float(pf.get_level_data('Enthalpy', level=2).max())
        self.assertAlmostEqual(max_enthalpy, 6.307355035152367, 6)

        # Test removing data
        pf.unload_data()
        self.assertEqual(pf.data_loaded, False)


        # Checkpoint files
        cf = PltFile(self.CHK_DATA_FILE, load_data=True)



    def test_plotting(self):
        pf = PltFile(self.DATA_FILE, load_data=True)
        porosity = pf.get_level_data('Porosity')
        plt.pcolormesh(porosity.x, porosity.y, porosity)



        self.assertEqual(pf.get_field_label('Porosity'), r'$\chi$')

    def test_diagnostics(self):
        pf = PltFile(self.DATA_FILE, load_data=True)

        pf.compute_diagnostic_vars()

        pf.compute_mush_liquid_interface()




if __name__ == "__main__":
    unittest.main()
