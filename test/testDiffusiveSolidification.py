# Run all convergence tests for the methods paper
import os, sys
from colorama import Fore, Style

from runAMRConvergenceTest import runTest
from mushyLayerRunUtils import get_base_output_dir, get_matlab_base_command

######################################
# 1) Diffusive solidification problem
#######################################
def testDiffusiveSolidification():

    base_output_dir = get_base_output_dir()
    matlab_command = get_matlab_base_command()

    print(Fore.GREEN + 'Setup tests for solidification without flow' + Style.RESET_ALL)

    physicalProblem = 'noFlow'
    AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': [8, 16, 32, 64, 128, 256]},
                {'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [8, 16, 32, 64, 128]},
                {'max_level': 2, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [8, 16, 32, 64]}]
    # While testing:
    # AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform']}];

    # Nzs 	  = [16, 32, 64]
    num_procs = [1] * 6  # Needs to be as long as the longest Nzs

    # Setup up the post processing command
    dataFolder = os.path.join(base_output_dir, 'NoFlow')
    figureName = os.path.join(dataFolder, 'noFlow.pdf')
    analysis_command = matlab_command + ' "Fig4NoFlow(\'' + dataFolder + '\', \'' + figureName + '\'); exit;"'

    # Run
    extra_params = {'main.debug': 'true'}
    runTest(dataFolder, physicalProblem, AMRSetup, num_procs, analysis_command, extra_params)



def main(argv):
    # try:
    #    opts, args = getopt.getopt(argv,"n:f:H:")
    # except getopt.GetoptError:
    #    print 'runAMRConvergenceTest.py -n <num processors>'
    #    sys.exit(2)
    #
    # for opt, arg in opts:
    #     if opt in ("-n"):
    #         num_proc = int(arg)

    testDiffusiveSolidification()

if __name__ == "__main__":
    main(sys.argv[1:])