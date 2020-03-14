import mushyLayerRunUtils
import os
import json
import logging
import subprocess
import pandas as pd
from chombopy.plotting import PltFile
from chombopy.inputs import write_inputs
import run_chombo_compare as chcompare

LOGGER = logging.getLogger(__name__)


class MushyLayerTest():
    EXPECTED = '.expected'
    PROPERTIES_FILE = 'properties.json'
    DIFF = 'diff-'
    DIFF_FOLDER = 'diffs'
    REQUIRED = ['inputs', PROPERTIES_FILE]

    def __init__(self, test_folder):
        self.test_folder = test_folder
        self.properties = {}
        self.has_run = False

        with open(os.path.join(self.test_folder, self.PROPERTIES_FILE)) as json_file:
            self.properties = json.load(json_file)

        try:
            self.mpi = subprocess.check_output(['which', 'mpiruna'])
            self.mpi_path = str(self.mpi.decode()).strip()
        except subprocess.CalledProcessError:
            self.mpi = None
            self.mpi_path = None

    def run(self):
        # Skip if parallel test and no mpirun
        if self.properties['proc'] > 1 and self.mpi is None:
            LOGGER.info('Parallel test but no MPI available, skipping test.')

        if self.mpi is not None:
            cmd = 'cd %s; %s -np %d %s inputs' % (self.test_folder, self.mpi_path, self.properties['proc'], self.get_mushy_layer_exec())
            os.system(cmd)

        else:
            cmd = 'cd %s; %s inputs > pout.0' % (self.test_folder, self.get_mushy_layer_exec())
            LOGGER.info(cmd)
            res = subprocess.check_output(cmd, shell=True)

            LOGGER.info('Response: "%s"' % str(res.decode()))

        self.has_run = True

    def get_diagnostics(self, file_name='diagnostics.csv'):
        """ Returns Pandas DataFrame of diagnostics """
        return pd.read_csv(os.path.join(self.test_folder, file_name))

    def get_latest_diagnostics(self):
        """ Returns Pandas DataFrame of latest diagnostics """
        return self.get_diagnostics('diagnosticsLatest.csv')

    def get_expected_latest_diags(self):
        return pd.read_csv(os.path.join(self.test_folder, 'diagnosticsLatest.csv.expected'))

    def get_plt_file(self, filename):
        return PltFile(os.path.join(self.test_folder, filename), load_data=True)

    def compare_latest_diags(self):
        return compare_diagnostics(self.get_latest_diagnostics(), self.get_expected_latest_diags())

    def get_mushy_layer_exec(self):
        exec_name = 'mushyLayer%dd' % self.properties['dim']
        exec_dir = os.path.join(mushyLayerRunUtils.get_mushy_layer_dir(), 'execSubcycle')
        mushy_layer_exec_path = mushyLayerRunUtils.get_executable_name(exec_dir, exec_name=exec_name,
                                                                       return_full_path=True)

        return mushy_layer_exec_path


# Static compare methods
def compare_diagnostics(data, data_expected, tolerance=1e-5):
    data_keys = data.keys()

    # diff = data.copy()
    # diff.to_csv(os.path.join(test_directory, 'diff-%s' % expected_file))

    for key in data_expected.keys():
        assert key in data_keys, 'Checking all diagnostics are present in both files'

        expected_value = data_expected[key].iloc[0]
        value_diff = float(data[key].iloc[0] - data_expected[key])

        # Compute relative different for large values
        if abs(expected_value) > 1.0:
            value_diff = value_diff / expected_value

        assert abs(value_diff) < tolerance, 'Comparing diagnostic: %s' % key


# Compare pout files
def compare_hdf5(computed_filename, exact_filename, test_directory, compare_exec, err_tolerance=1e-5):
    """ computed_filename and exact_filename relative to test_directory
    """
    diffs_folder = os.path.join(test_directory, 'diffs')
    if not os.path.exists(diffs_folder):
        os.makedirs(diffs_folder)
        LOGGER.info('Making folder %s' % diffs_folder)

    # chombo_dir = os.environ['CHOMBO_HOME']
    # compare_dir = os.path.join(chombo_dir, 'util', 'ChomboCompare')
    # compare_exec = mushyLayerRunUtils.get_executable_name(compare_dir, 'compare%dd' % properties['dim'])
    # compare_exec = os.path.join(compare_dir, compare_exec)
    # LOGGER.info('Found executable %s ' % compare_exec)

    if not os.path.exists(compare_exec):
        LOGGER.warning('Could not find Chombo compare executable %s' % compare_exec)
        LOGGER.warning('So cannot compare hdf5 output files')

    compare_params_file = os.path.join(test_directory, 'compare.inputs')

    computed_file = os.path.join(test_directory, computed_filename)
    exact_file = os.path.join(test_directory, exact_filename)
    error_file = os.path.join(diffs_folder, 'diff-' + computed_filename)
    compare_params = {'compare.sameSize': 1,
                      'compare.exactRoot': exact_file,
                      'compare.computedRoot': computed_file,
                      'compare.errorRoot': error_file,
                      'compare.doPlots': 1,
                      'compare.HOaverage': 0,
                      'compare.no_average_var': 'T err'}

    write_inputs(compare_params_file, compare_params)
    cmd = 'cd %s ; %s %s > pout.0' % (diffs_folder, compare_exec, compare_params_file)

    LOGGER.info('Executing: %s' % cmd)
    res = subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT)
    LOGGER.info('Compare command run, response: %s' % res)

    # Rename pout.0 in case we make lots
    old_pout_name = os.path.join(diffs_folder, 'pout.0')
    new_pout_name = os.path.join(diffs_folder, 'pout-%s.0' % computed_filename)
    LOGGER.debug('Rename %s to %s' % (old_pout_name, new_pout_name))
    if not os.path.exists(old_pout_name):
        LOGGER.warning('Cannot find %s' % old_pout_name)

    os.rename(old_pout_name, new_pout_name)

    # Now check the output diff
    errs = chcompare.load_error_file(new_pout_name)
    LOGGER.debug(errs)

    for field in errs.keys():
        field_errs = errs[field]
        for err_type in field_errs.keys():
            this_field_err = field_errs[err_type]
            assert abs(this_field_err) < err_tolerance, 'Comparing field: %s' % field