import unittest
from PltFile import PltFile, latexify, latexify2
import matplotlib.pyplot as plt

class TestPltFile(unittest.TestCase):
    DATA_FILE = 'data/plt000100.2d.hdf5'
    CHK_DATA_FILE = 'data/chk000100.2d.hdf5'

    DATA_FILE_3D = 'data/3D/plt000100.3d.hdf5'

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


        pf_no_name = PltFile('data/pltnoframe.2d.hdf5')
        self.assertEqual(pf_no_name.frame, -1)

        pf_no_inputs = PltFile('data/plt000100.2d.hdf5', inputs_file='does not exist')
        self.assertIsNone(pf_no_inputs.inputs)

        # Checkpoint files
        cf = PltFile(self.CHK_DATA_FILE, load_data=True)



    def test_plotting(self):
        pf = PltFile(self.DATA_FILE)
        pf.load_data(zero_x=True)
        porosity = pf.get_level_data('Porosity')
        latexify()

        fig = plt.figure()
        ax = fig.gca()
        ax.pcolormesh(porosity.x, porosity.y, porosity)

        pf.plot_outlines(ax)

        self.assertEqual(pf.get_field_label('Porosity'), r'$\chi$')

        latexify2(5.0, 3.0)
        pf.plot_field('Porosity')

    def test_diagnostics(self):
        pf = PltFile(self.DATA_FILE, load_data=True)

        pf.compute_diagnostic_vars()

        pf.compute_mush_liquid_interface()

        properties = pf.channel_properties()

        num_channels = pf.num_channels(0.9)
        self.assertEqual(num_channels, 2)


        permeability = pf.get_permeability()
        self.assertAlmostEqual(permeability.max(), pf.inputs['parameters.heleShawPermeability'], 5)

        permeability = pf.get_permeability('cubic')
        self.assertEqual(permeability.max(), 1.0)

        permeability = pf.get_permeability('anything else')
        self.assertEqual(permeability, 1.0)

    def test_3d(self):
        pf = PltFile(self.DATA_FILE_3D, load_data=True)

        porosity = pf.get_level_data('Porosity')
        coords = porosity.coords
        x = coords['x']
        y = coords['y']
        z = coords['y']

        self.assertEqual(pf.space_dim, 3)
        self.assertEqual(len(z), 32)
        self.assertEqual(len(x), 32)
        self.assertEqual(len(y), 32)



    def test_static_methods(self):
        pf = PltFile(self.DATA_FILE, load_data=True)

        x,y = pf.get_mesh_grid()
        self.assertEqual(len(x), 16)

        x,y = pf.get_mesh_grid(extend_grid=False)
        self.assertEqual(len(x), 15)

        porosity = pf.get_level_data('Porosity')
        x,y = PltFile.get_mesh_grid_n(porosity)
        self.assertEqual(len(x), 16)



if __name__ == "__main__":
    unittest.main()
