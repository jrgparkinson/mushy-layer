import os
import sys


def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


if __name__ == "__main__":

    # Add folders to PYTHONPATH
    required_paths = [os.path.join(get_script_path(), 'test'),
                      os.path.join(get_script_path(), 'plotting')]

    for required_path in required_paths:
        if 'PYTHONPATH' not in os.environ or required_path not in os.environ['PYTHONPATH']:
            os.environ['PYTHONPATH'] += ':' + required_path
    print('PYTHONPATH=%s' % os.environ['PYTHONPATH'])

    # Setup MUSHY_LAYER_DIR
    if 'MUSHY_LAYER_DIR' not in os.environ:
        os.environ['MUSHY_LAYER_DIR'] = get_script_path()
    print('MUSHY_LAYER_DIR=%s' % os.environ['MUSHY_LAYER_DIR'])



