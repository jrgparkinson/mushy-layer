import unittest
from PltFile import PltFile

class TestPltFile(unittest.TestCase):
    DATA_FILE = 'data/plt000100.2d.hdf5'

    def test_load(self):
        pf = PltFile(self.DATA_FILE, load_data=True)

        self.assertEqual(pf.num_levels, 3)








if __name__ == "__main__":
    unittest.main()
