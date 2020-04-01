from regression.mushy_layer_test import MushyLayerTest
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def test_darcy():
    darcy_test = MushyLayerTest(os.path.join(THIS_DIR, 'AMRDarcy'))

    darcy_test.run()

    # Check latest diagnostics are the same
    darcy_test.compare_latest_diags()


def test_darcy_brinkman():
    db_test = MushyLayerTest(os.path.join(THIS_DIR, 'AMRDarcyBrinkman'))

    db_test.run()

    # Check latest diagnostics are the same
    db_test.compare_latest_diags()
