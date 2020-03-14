from regression.mushy_layer_test import MushyLayerTest
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def test_darcy():
    ml_test = MushyLayerTest(os.path.join(THIS_DIR, 'AMRRestartDarcy'))

    ml_test.run()

    # Check latest diagnostics are the same
    ml_test.compare_latest_diags()

def test_darcy_brinkman():
    ml_test = MushyLayerTest(os.path.join(THIS_DIR, 'AMRRestartDarcyBrinkman'))

    ml_test.run()

    # Check latest diagnostics are the same
    ml_test.compare_latest_diags()

