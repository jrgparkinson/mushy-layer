import os
import re
import json
import difflib
import subprocess
import traceback
import time
import getopt
import sys
from colorama import Fore
import pandas as pd

# Two imports from the directory above this one - make sure /path/to/mushy-layer/test is in your python path
import mushyLayerRunUtils
import run_chombo_compare as chcompare


# Error tolerance for two data files to be considered equal
err_tolerance = 1e-10


class Logger():

    def __init__(self, log_file_path):
        self.log_file = open(log_file_path, 'w+')
        self.verbose = False

    def set_verbose(self, v):
        self.verbose = v

    def log(self, text, console_display=False):
        if self.verbose or console_display:
            print(text)

        self.log_file.write(text + '\n')

    def logl(self, text, console_display=False):
        if self.verbose or console_display:
            sys.stdout.write(text)
            sys.stdout.flush()

        self.log_file.write(text)

    def log_failed(self):
        self.log_status(Fore.RED + '[Failed]' + Fore.RESET)

    def log_ok(self):
        self.log_status(Fore.LIGHTGREEN_EX + '[OK]' + Fore.RESET)

    def log_void(self):
        self.log_status(Fore.YELLOW + '[Void]' + Fore.RESET)

    def log_status(self, text):
        formatted = '%-8s' % text
        sys.stdout.write(formatted)
        sys.stdout.flush()

        self.log_file.write(formatted)


def filter_pout(text_lines):
    """
    Remove wallclock time reporting from pout files so they can be compared sensibly.
    The wallclock time will always vary, so we don't want to include it in the comparison.
    :param text_lines: array of lines of text in the file
    :return: filtered_output: array of lines of text
    """

    filtered_output = []
    for line in text_lines:
        new_line = re.sub(r'\s+wallclocktime = \d+\.\d+e[+-]\d+\s?\n', r'\n', line)
        filtered_output.append(new_line)

        # Once we've had the 'total number of points updated' line, stop processing
        if 'total number of points updated = ' in line:
            break

    return filtered_output


def test_folder(test_directory, verbose_output=False):
    """
    Run the regression test contained within a folder. Expects the folder o have the following files:
    - inputs
    - properties.json
    - at least one file ending with .expected to compare with the computed ouptut,
        e.g. pout.0.expected, or plt000010.2d.hdf5.expected
    :param test_directory:
    :return: success - if test was successful or not
    """

    logger.log('Processing folder "%s"' % test_directory)

    remove_existing_diffs_cmd = 'rm %s/diff-*' % test_directory
    logger.log('Remove existing diffs: ' + remove_existing_diffs_cmd)
    os.system(remove_existing_diffs_cmd)

    test_files = os.listdir(test_directory)

    logger.log('Initial files in folder ' + str(test_files))



    # Check the required files exist
    # if not, skip this test
    if not all(elem in test_files for elem in REQUIRED):
        test_name = test_directory.split('/')[-1]
        logger.logl('%-25s    ' % test_name, console_display=True)

        logger.log('  Required files for conducting a test (%s) not found in %s' % (REQUIRED, test_directory))
        logger.log('  Skipping test \n')

        logger.log_void()
        return False, 'Void'

    # Load properties
    with open(os.path.join(test_directory, PROPERTIES_FILE)) as json_file:
        properties = json.load(json_file)

    # logger.log('==Running test: %s==' % properties['name'])
    logger.logl('%-25s    ' % properties['name'], console_display=True)

    try:
        mpi = subprocess.check_output(['which', 'mpiruna'])
        mpi_path = str(mpi.decode()).strip()
    except subprocess.CalledProcessError:
        mpi = None
        mpi_path = None

    # Skip if parallel test and no mpirun
    if properties['proc'] > 1 and mpi is None:
        logger.log_void()
        return False, 'Void'

    # Get correct executable
    # mushy_layer_exec = mushyLayerRunUtils.get_executable(dim=properties['dim'])
    # mushy_layer_exec_path = os.path.join(mushyLayerRunUtils.get_mushy_layer_dir(), 'execSubcycle', mushy_layer_exec)
    exec_name = 'mushyLayer%dd' % properties['dim']
    exec_dir = os.path.join(mushyLayerRunUtils.get_mushy_layer_dir(), 'execSubcycle')
    mushy_layer_exec_path = mushyLayerRunUtils.get_executable_name(exec_dir, exec_name=exec_name, return_full_path=True)

    if not os.path.exists(mushy_layer_exec_path):
        logger.log('\n**Could not find mushy layer executable: %s' % mushy_layer_exec_path)
        logger.log('**Have you compiled the code for the right number of dimensions?')
        logger.logl('**Use \'make all DIM=3\' to compile in 3D    ')

        logger.log_failed()
        return False, 'Failed'

    # Run test
    # try:
    if mpi is not None:
        cmd = 'cd %s; %s -np %d %s inputs' % (test_directory, mpi_path, properties['proc'], mushy_layer_exec_path)
        os.system(cmd)
        # mpi_path = 'mpirun'
        # res = subprocess.run([mpi_path, '-n ', str(properties['proc']), ' inputs'], cwd=test_directory)

    else:
        cmd = 'cd %s; %s inputs > pout.0' % (test_directory, mushy_layer_exec_path)
        logger.log(cmd)
        res = subprocess.check_output(cmd, shell=True)

        logger.logl('Response: "%s"' % str(res.decode()))

    # os.system(cmd)

    # except subprocess.CalledProcessError:
    #     logger.logl('[Exception]')
    #     logger.log_failed()
    #     return False, 'Failed'



    # Compare output against the expected output
    # logger.log('Test files: %s' % test_files)

    expected_files = [f for f in test_files if EXPECTED in f and DIFF not in f]
    if not expected_files:
        logger.log('No expected files to compared against')
        logger.log_void()
        return False, 'Void'

    failed_test = False
    for expected_file in expected_files:
        logger.log('Expected file: %s' % expected_file)

        # Check an output file to compare against exists
        test_output_filename = expected_file.replace(EXPECTED, '')
        test_output_file_path = os.path.join(test_directory, test_output_filename)
        if not os.path.exists(test_output_file_path):
            logger.log('No output file generated to compare against: %s' % test_output_file_path)
            logger.log_failed()
            return False, 'Failed'

        if '.hdf5' in expected_file:
            # Do chombo compare
            # logger.log('Using chombo compare')

            # Make diffs folder if it doesn't exist
            diffs_folder = os.path.join(test_directory, DIFF_FOLDER)
            if not os.path.exists(diffs_folder):
                os.makedirs(diffs_folder)
                logger.log('Making folder %s' % diffs_folder)

            chombo_dir = os.environ['CHOMBO_HOME']
            compare_dir = os.path.join(chombo_dir, 'util', 'ChomboCompare')
            compare_exec = mushyLayerRunUtils.get_executable_name(compare_dir, 'compare%dd' % properties['dim'])
            compare_exec = os.path.join(compare_dir, compare_exec)
            # logger.log('Found executable %s ' % compare_exec)

            if not os.path.exists(compare_exec):
                logger.log('Could not find Chombo compare executable %s' % compare_exec)
                logger.log('So cannot compare hdf5 output files')

            compare_params_file = os.path.join(test_directory, 'compare.inputs')

            computed_file = os.path.join(test_directory, test_output_file_path)
            exact_file = os.path.join(test_directory, expected_file)
            error_file = os.path.join(diffs_folder, DIFF + test_output_filename)
            compare_params = {'compare.sameSize': 1,
                              'compare.exactRoot': exact_file,
                              'compare.computedRoot': computed_file,
                              'compare.errorRoot': error_file,
                              'compare.doPlots': 1,
                              'compare.HOaverage': 0,
                              'compare.no_average_var': 'T err'}

            mushyLayerRunUtils.write_inputs(compare_params_file, compare_params)
            cmd = 'cd %s ; %s %s > pout.0' % (diffs_folder, compare_exec, compare_params_file)

            logger.log('Executing: %s' % cmd)
            res = subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT)
            logger.log('Compare command run, response: %s' % res)

            # Rename pout.0 in case we make lots
            old_pout_name = os.path.join(diffs_folder, 'pout.0')
            new_pout_name = os.path.join(diffs_folder, 'pout-%s.0' % test_output_filename)
            # logger.log('Rename %s to %s' % (old_pout_name, new_pout_name))
            if not os.path.exists(old_pout_name):
                logger.log('Cannot find %s' % old_pout_name, console_display=True)

            os.rename(old_pout_name, new_pout_name)

            # Now check the output diff
            errs = chcompare.load_error_file(new_pout_name)
            # logger.log(errs)

            for field in errs.keys():
                field_errs = errs[field]
                # logger.log(field_errs)
                for err_type in field_errs.keys():
                    this_field_err = field_errs[err_type]
                    if abs(this_field_err) > err_tolerance:
                        logger.log('Error in field %s is non-zero' % field)
                        logger.log('See %s and %s for more info' % (new_pout_name, error_file))
                        # logger.log_failed()
                        # return False, 'Failed'
                        failed_test = True

        elif '.csv' in expected_file:
            data = pd.read_csv(os.path.join(test_directory, test_output_file_path))
            data_expected = pd.read_csv(os.path.join(test_directory, expected_file))

            data_keys = data.keys()
            
            diff = data.copy()

            for key in data_expected.keys():
                if key not in data_keys:
                    logger.log('Key not found: %s' % key, verbose_output)
                    failed_test = True
                    break

                tolerance = 1e-5
                expected_value = data_expected[key].iloc[0]
                value_diff = float(data[key].iloc[0] - data_expected[key])

                # Compute relative different for large values
                if abs(expected_value) > 1.0:
                    value_diff = value_diff/expected_value

                diff[key] = value_diff
                if abs(value_diff) > tolerance:
                    logger.log('Error in field %s = %.2g' % (key, value_diff), verbose_output)
                    failed_test = True

            diff.to_csv(os.path.join(test_directory, 'diff-%s' % expected_file))

        else:
            # Assume text file - do a diff
            text1 = open(os.path.join(test_directory, expected_file)).readlines()
            text2 = open(os.path.join(test_directory, test_output_file_path)).readlines()

            text1 = filter_pout(text1)
            text2 = filter_pout(text2)

            diff = difflib.unified_diff(text1, text2)
            # diff = difflib.ndiff(text1, text2)
            differences = [line for line in diff]

            diff_out_file = os.path.join(test_directory, DIFF + test_output_filename)

            with open(diff_out_file, 'w') as diff_file:
                diff_file.writelines(differences)

            logger.log('Diff: %s' % differences)

            if differences:
                logger.log('Differences found in %s' % test_output_file_path)
                logger.log('For details, see %s' % diff_out_file)
                # logger.log('** Test failed \n')
                # logger.log_failed()
                # return False, 'Failed'
                # failed_test = True

    if failed_test:
        logger.log_failed()
        return False, 'Failed'
    else:
        logger.log_ok()
        return True, 'OK'


def usage():
    logger.log('python run_regression_tests.py [-t <test_dir>] [-v verbose] ', console_display=True)


if __name__ == "__main__":

    timings = [time.time()]

    # Initial setup

    # Constants
    EXPECTED = '.expected'
    PROPERTIES_FILE = 'properties.json'
    DIFF = 'diff-'
    DIFF_FOLDER = 'diffs'
    REQUIRED = ['inputs', PROPERTIES_FILE]

    # Get full path to this script
    script_loc = os.path.dirname(os.path.realpath(__file__))
    logger = Logger(os.path.join(script_loc, 'regression.log'))

    # By default assume every subdirectory is a test
    test_dirs = [f for f in os.listdir(script_loc) if os.path.isdir(os.path.join(script_loc, f))]
    verbose = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "ht:v", ["help", "test="])
    except getopt.GetoptError as err:
        # print help information and exit:
        logger.log(str(err), console_display=True)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    for o, arg in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-t", "--test"):
            test_dirs = [arg]
        elif o in "-v":
            verbose = True
        else:
            logger.log('Unknown command line option: %s' % o, console_display=True)
            sys.exit(-1)

    logger.set_verbose(verbose)
    logger.log('Tests to run: ' + ','.join(test_dirs), console_display=True)

    # Containers to collate passed/failed tests
    failed_tests = []
    passed_tests = []
    test_results = {}

    # Timing
    timings.append(time.time())

    for test_dir in test_dirs:

        status = 'Failed'

        # noinspection PyBroadException
        try:
            full_dir = os.path.join(script_loc, test_dir)
            success, status = test_folder(full_dir, verbose)
        except Exception as e:
            logger.log(traceback.format_exc())
            logger.log_failed()
            success = False

        # Timing
        timings.append(time.time())
        logger.log(' (%.3g seconds)' % (timings[-1] - timings[-2]), console_display=True)

        if success:
            passed_tests.append(test_dir)
        else:
            failed_tests.append(test_dir)

        test_results[test_dir] = status

    logger.log('', console_display=True)
    logger.log('==========================================================', console_display=True)
    logger.log('------------------------ Summary -------------------------', console_display=True)

    # Report tests
    void_tests = [k for (k, v) in test_results.items() if v == 'Void']
    passed_tests = [k for (k, v) in test_results.items() if v == 'OK']
    failed_tests = [k for (k, v) in test_results.items() if v == 'Failed']
    num_passed = len(passed_tests)
    num_failed = len(failed_tests)
    num_tests = len(test_results)
    logger.log('Tests passed (%d/%d): %s' % (num_passed, num_tests, ','.join(passed_tests)), console_display=True)
    logger.log('Tests failed (%d/%d): %s' % (num_failed, num_tests, ','.join(failed_tests)), console_display=True)
    logger.log('Tests void (%d/%d): %s' % (len(void_tests), num_tests, ','.join(void_tests)),
               console_display=True)

    # Timing
    timings.append(time.time())
    logger.log('Total runtime: %.3g seconds' % (timings[-1] - timings[0]), console_display=True)

    # Give a bad exit if we failed a test
    if num_failed > 0:
        sys.exit(1)
