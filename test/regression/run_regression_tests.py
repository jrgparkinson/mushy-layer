import os
import re
import json
import difflib
import subprocess
import traceback
import logging
import time
import getopt
import sys
from colorama import Fore

# Two imports from the directory above this one - make sure /path/to/mushy-layer/test is in your python path
import mushyLayerRunUtils
import run_chombo_compare as chcompare


class Logger():

    def __init__(self, log_file_path):
        self.log_file = open(log_file_path, 'w+')

    def log(self, text):
        print(text)
        self.log_file.write(text + '\n')

    def logl(self, text):
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
        new_line = re.sub(r'wallclocktime = \d+\.\d+e[+-]\d+\s?\n', '', line)
        filtered_output.append(new_line)

    return filtered_output


def test_folder(test_directory, verbose_output=False):
    """
    Run the regression test contained within a folder. Expects the folder o have the following files:
    - inputs
    - properties.json
    - at least one file ending with .expected to compare with the computed ouptut,
        e.g. pout.0.expected, or plt000010.2d.hdf5.expected
    :param test_directory:
    :return: success - if test was succesful or not
    """
    test_files = os.listdir(test_directory)

    # Check the required files exist
    # if not, skip this test
    if not all(elem in test_files for elem in REQUIRED):
        if verbose_output:
            logger.log('  Required files for conducting a test (%s) not found in %s' % (REQUIRED, test_directory))
            logger.log('  Skipping test \n')
        else:
            test_name = test_directory.split('/')[-1]
            logger.logl('%-25s    ' % test_name)
            logger.log_void()
        return False, 'Void'

    # Load properties
    with open(os.path.join(test_directory, PROPERTIES_FILE)) as json_file:
        properties = json.load(json_file)

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

    # logger.log('==Running test: %s==' % properties['name'])
    logger.logl('%-25s    ' % properties['name'])

    # Get correct executable
    # mushy_layer_exec = mushyLayerRunUtils.get_executable(dim=properties['dim'])
    # mushy_layer_exec_path = os.path.join(mushyLayerRunUtils.get_mushy_layer_dir(), 'execSubcycle', mushy_layer_exec)
    exec_name = 'mushyLayer%dd' % properties['dim']
    exec_dir = os.path.join(mushyLayerRunUtils.get_mushy_layer_dir(), 'execSubcycle')
    mushy_layer_exec_path = mushyLayerRunUtils.get_executable_name(exec_dir, exec_name=exec_name, return_full_path=True)

    if not os.path.exists(mushy_layer_exec_path):
        if verbose_output:
            logger.log('\n**Could not find mushy layer executable: %s' % mushy_layer_exec_path)
            logger.log('**Have you compiled the code for the right number of dimensions?')
            logger.logl('**Use \'make all DIM=3\' to compile in 3D    ')
        else:
            logger.log_failed()
        return False, 'Failed'

    # Run test
    # try:
    if mpi is not None:
        cmd = 'cd %s; %s -np %d %s inputs' % (test_directory, mpi_path, properties['proc'], mushy_layer_exec_path)
        # mpi_path = 'mpirun'
        # res = subprocess.run([mpi_path, '-n ', str(properties['proc']), ' inputs'], cwd=test_directory)

    else:
        cmd = 'cd %s; %s inputs' % (test_directory, mushy_layer_exec_path)
        # res = subprocess.check_output([mushy_layer_exec_path, 'inputs'], cwd=test_directory)

    os.system(cmd)

    # except subprocess.CalledProcessError:
    #     logger.logl('[Exception]')
    #     logger.log_failed()
    #     return False, 'Failed'



    # Compare output against the expected output
    # logger.log('Test files: %s' % test_files)

    expected_files = [f for f in test_files if EXPECTED in f]
    if not expected_files:
        if verbose_output:
            logger.log('No expected files to compared against')
        else:
            logger.log_void()
        return False, 'Void'

    for expected_file in expected_files:

        # Check an output file to compare against exists
        test_output_filename = expected_file.replace(EXPECTED, '')
        test_output_file_path = os.path.join(test_directory, test_output_filename)
        if not os.path.exists(test_output_file_path):
            if verbose_output:
                logger.log('No output file generated to compare against: %s' % test_output_file_path)
            else:
                logger.log_failed()
            return False, 'Failed'

        if '.hdf5' in expected_file:
            # Do chombo compare
            # logger.log('Using chombo compare')

            # Make diffs folder if it doesn't exist
            diffs_folder = os.path.join(test_directory, DIFF_FOLDER)
            if not os.path.exists(diffs_folder):
                os.makedirs(diffs_folder)

            chombo_dir = os.environ['CHOMBO_HOME']
            compare_dir = os.path.join(chombo_dir, 'util', 'ChomboCompare')
            compare_exec = mushyLayerRunUtils.get_executable_name(compare_dir, 'compare%dd' % properties['dim'])
            compare_exec = os.path.join(compare_dir, compare_exec)
            # logger.log('Found executable %s ' % compare_exec)

            if not os.path.exists(compare_exec) and verbose_output:
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
            cmd = 'cd %s ; %s %s  > /dev/null' % (diffs_folder, compare_exec, compare_params_file)

            if verbose_output:
                logger.log('Executing: %s' % cmd)
            # This prints the console output, which can be messy:
            # os.system(cmd)

            # This doesn't print the console output, which is cleaner
            with open(os.devnull, 'wb') as devnull:
                subprocess.check_call(cmd, shell=True, stdout=devnull, stderr=subprocess.STDOUT)

            # Rename pout.0 in case we make lots
            old_pout_name = os.path.join(diffs_folder, 'pout.0')
            new_pout_name = os.path.join(diffs_folder, 'pout-%s.0' % test_output_filename)
            # logger.log('Rename %s to %s' % (old_pout_name, new_pout_name))
            os.rename(old_pout_name, new_pout_name)

            # Now check the output diff
            errs = chcompare.load_error_file(new_pout_name)
            # logger.log(errs)

            for field in errs.keys():
                field_errs = errs[field]
                # logger.log(field_errs)
                for err_type in field_errs.keys():
                    this_field_err = field_errs[err_type]
                    if abs(this_field_err) > 1e-10:
                        if verbose_output:
                            logger.log('Error in field %s is non-zero' % field)
                            logger.log('See %s and %s for more info' % (new_pout_name, error_file))
                        else:
                            logger.log_failed()
                        return False, 'Failed'
        else:
            # Assume text file - do a diff
            text1 = open(os.path.join(test_directory, expected_file)).readlines()
            text2 = open(os.path.join(test_directory, test_output_file_path)).readlines()

            text1 = filter_pout(text1)
            text2 = filter_pout(text2)

            diff = difflib.unified_diff(text1, text2)
            # diff = difflib.ndiff(text1, text2)

            diff_out_file = os.path.join(test_directory, DIFF + test_output_filename)

            with open(diff_out_file, 'w') as diff_file:
                diff_file.writelines(diff)

            differences = [line for line in diff]
            # logger.log('Diff: %s' % differences)

            if differences:
                if verbose_output:
                    logger.log('Differences found in %s' % test_output_file_path)
                    logger.log('For details, see %s' % diff_out_file)
                    logger.log('** Test failed \n')
                else:
                    logger.log_failed()
                return False, 'Failed'

    logger.log_ok()
    return True, 'OK'


def usage():
    logger.log('python run_regression_tests.py [-t <test_dir>] [-v verbose] ')


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
        logger.log(str(err))  # will print something like "option -a not recognized"
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
            logger.log('Unknown command line option: %s' % o)
            sys.exit(-1)

    logger.log('Tests to run: ' + ','.join(test_dirs))

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

            if verbose:
                logger.log('**Error running test**')
                logging.error(traceback.format_exc())
            else:
                logger.log_failed()
            success = False

        # Timing
        timings.append(time.time())
        logger.log(' (%.3g seconds)' % (timings[-1] - timings[-2]))

        if success:
            passed_tests.append(test_dir)
        else:
            failed_tests.append(test_dir)

        test_results[test_dir] = status

    logger.log('')
    logger.log('==========================================================')
    logger.log('------------------------ Summary -------------------------')

    # Report tests
    num_passed = len(passed_tests)
    num_failed = len(failed_tests)
    logger.log('Tests passed (%d/%d): %s' % (num_passed, num_passed + num_failed, ','.join(passed_tests)))
    logger.log('Tests failed (%d/%d): %s' % (num_failed, num_passed + num_failed, ','.join(failed_tests)))

    # Timing
    timings.append(time.time())
    logger.log('Total runtime: %.3g seconds' % (timings[-1] - timings[0]))

    # Give a bad exit if we failed a test
    if num_failed > 0:
        sys.exit(1)
