import os
import re
import json
import difflib
import subprocess
import traceback
import logging
import time
import getopt, sys

# Two imports from the directory above this one - make sure /path/to/mushy-layer/test is in your python path
import mushyLayerRunUtils as utils
import run_chombo_compare as chcompare

def filter_pout(text_lines):
    '''
    Remove wallclock time reporting from pout files so they can be compared sensibly.
    The wallclock time will always vary, so we don't want to include it in the comparison.
    :param text_lines: array of lines of text in the file
    :return: filtered_output: array of lines of text
    '''

    filtered_output = []
    for line in text_lines:
        new_line =  re.sub('wallclocktime = \d+\.\d+e[+-]\d+\s?\n', '', line)
        filtered_output.append(new_line)

    return filtered_output


def test_folder(test_dir):
    '''
    Run the regression test contained within a folder. Expects the folder o have the following files:
    - inputs
    - properties.json
    - at least one file ending with .expected to compare with the computed ouptut,
        e.g. pout.0.expected, or plt000010.2d.hdf5.expected
    :param test_dir:
    :return: success - if test was succesful or not
    '''
    test_files = os.listdir(test_dir)

    # Check the required files exist
    # if not, skip this test
    if not all(elem in test_files for elem in REQUIRED):
        print('Required files for conducting a test (%s) not found in %s' % (REQUIRED, test_dir))
        print('Skipping test \n')
        return False

    # Load properties
    with open(os.path.join(test_dir, PROPERTIES_FILE)) as json_file:
        properties = json.load(json_file)

    print('==Running test: %s==' % properties['name'])

    # Get correct executable
    mushy_layer_exec = utils.get_executable(dim=properties['dim'])
    mushy_layer_exec_path = os.path.join(utils.get_mushy_layer_dir(), 'execSubcycle', mushy_layer_exec)

    if not os.path.exists(mushy_layer_exec_path):
        print('Could not find mushy layer executable: %s' % mushy_layer_exec_path)
        print('Have you compiled the code for the right number of dimensions?')
        print('Use \'make all DIM=3\' to compile in 3D')
        return False

    # Run test
    cmd = 'cd %s; mpirun -np %d %s inputs' % (test_dir, properties['proc'], mushy_layer_exec_path)
    # print('Executing: %s' % cmd)
    os.system(cmd)

    # Compare output against the expected output
    # print('Test files: %s' % test_files)

    expected_files = [f for f in test_files if EXPECTED in f]
    if not expected_files:
        print('No expected files to compared against')
        return False

    for expected_file in expected_files:

        # Check an output file to compare against exists
        test_output_filename = expected_file.replace(EXPECTED, '')
        test_output_file_path = os.path.join(test_dir, test_output_filename)
        if not os.path.exists(test_output_file_path):
            print('No output file generated to compare against: %s' % test_output_file_path)
            return False

        if '.hdf5' in expected_file:
            # Do chombo compare
            # print('Using chombo compare')

            # Make diffs folder if it doesn't exist
            diffs_folder = os.path.join(test_dir, DIFF_FOLDER)
            if not os.path.exists(diffs_folder):
                os.makedirs(diffs_folder)

            chombo_dir = os.environ['CHOMBO_HOME']
            compare_dir = os.path.join(chombo_dir, 'lib', 'util', 'ChomboCompare')
            compare_exec = utils.get_executable_name(compare_dir, 'compare%dd' % properties['dim'])
            compare_exec = os.path.join(compare_dir, compare_exec)
            # print('Found executable %s ' % compare_exec)

            if not os.path.exists(compare_exec):
                print('Could not find Chombo compare executable %s' % compare_exec)
                print('So cannot compare hdf5 output files')

            compare_params_file = os.path.join(test_dir, 'compare.inputs')

            computed_file = os.path.join(test_dir, test_output_file_path)
            exact_file = os.path.join(test_dir, expected_file)
            error_file = os.path.join(diffs_folder, DIFF + test_output_filename)
            compare_params = {'compare.sameSize': 0,
                              'compare.exactRoot': exact_file,
                              'compare.computedRoot': computed_file,
                              'compare.errorRoot': error_file,
                              'compare.doPlots': 1,
                              'compare.HOaverage': 0,
                              'compare.no_average_var': 'T err'}

            utils.write_inputs(compare_params_file, compare_params)
            cmd = 'cd %s ; %s %s  > /dev/null' % (diffs_folder, compare_exec, compare_params_file)

            # print('Executing: %s' % cmd)
            # This prints the console output, which can be messy:
            # os.system(cmd)

            # This doesn't print the console output, which is cleaner
            with open(os.devnull, 'wb') as devnull:
                subprocess.check_call(cmd, shell=True, stdout=devnull, stderr=subprocess.STDOUT)

            # Rename pout.0 in case we make lots
            old_pout_name = os.path.join(diffs_folder, 'pout.0')
            new_pout_name = os.path.join(diffs_folder, 'pout-%s.0' % test_output_filename)
            # print('Rename %s to %s' % (old_pout_name, new_pout_name))
            os.rename(old_pout_name, new_pout_name)

            # Now check the output diff
            errs = chcompare.load_error_file(new_pout_name)
            # print(errs)

            for field in errs.keys():
                field_errs = errs[field]
                # print(field_errs)
                for err_type in field_errs.keys():
                    err = field_errs[err_type]
                    if abs(err) > 1e-10:
                        print('Error in field %s is non-zero' % field)
                        print('See %s and %s for more info' % (new_pout_name, error_file))
                        return False




        else:
            # Assume text file - do a diff


            text1 = open(os.path.join(test_dir, expected_file)).readlines()
            text2 = open(os.path.join(test_dir, test_output_file_path)).readlines()

            text1 = filter_pout(text1)
            text2 = filter_pout(text2)

            diff = difflib.unified_diff(text1, text2)
            # diff = difflib.ndiff(text1, text2)

            diff_out_file = os.path.join(test_dir, DIFF + test_output_filename)

            with open(diff_out_file, 'w') as diff_file:
                diff_file.writelines(diff)

            differences = [line for line in diff]
            # print('Diff: %s' % differences)

            if differences:
                print('Differences found in %s' % test_output_file_path)
                print('For details, see %s' % diff_out_file)
                print('** Test failed \n')
                return False

    return True


def usage():
    print('python run_regression_tests.py [-t test_dir]')

if __name__ == "__main__":

    timings = []
    timings.append(time.time())

    ## Initial setup

    # Constants
    EXPECTED = '.expected'
    PROPERTIES_FILE = 'properties.json'
    DIFF = 'diff-'
    DIFF_FOLDER = 'diffs'
    REQUIRED = ['inputs', PROPERTIES_FILE]

    # Get full path to this script
    script_loc = os.path.dirname(os.path.realpath(__file__))

    # By default assume every subdirectory is a test
    test_dirs = [f for f in os.listdir(script_loc) if os.path.isdir(os.path.join(script_loc, f))]

    try:
        opts, args = getopt.getopt(sys.argv[1:], "ht:", ["help", "test="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    for o, arg in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-t", "--test"):
            test_dirs = [arg]
        else:
            assert False, "unhandled option"
            
    print('Tests to run: ' + ','.join(test_dirs))



    # Containers to collate passed/failed tests
    failed_tests = []
    passed_tests = []

    #Timing
    timings.append(time.time())
    print('Setup time: %.3g seconds' % (timings[-1] - timings[-2]))

    for test_dir in test_dirs:

        try:
            full_dir = os.path.join(script_loc, test_dir)
            success = test_folder(full_dir)

            # Timing
            timings.append(time.time())
            print('Test runtime: %.3g seconds' % (timings[-1] - timings[-2]))

        except Exception as e:

            print('**Error running test**')
            logging.error(traceback.format_exc())
            success = False

        if success:
            print('**Test passed** \n')
            passed_tests.append(test_dir)
        else:
            print('**Test failed** \n')
            failed_tests.append(test_dir)


    print('==========================================================')
    print('------------------------ Summary -------------------------')

    # Report tests
    num_passed = len(passed_tests)
    num_failed = len(failed_tests)
    print('Tests passed (%d/%d): %s' % (num_passed, num_passed+num_failed, ','.join(passed_tests)))
    print('Tests failed (%d/%d): %s' % (num_failed, num_passed + num_failed, ','.join(failed_tests)))

    # Timing
    timings.append(time.time())
    print('Total runtime: %.3g seconds' % (timings[-1] - timings[0]))