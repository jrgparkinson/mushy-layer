# Run all convergence tests for the methods paper
import os, sys
from colorama import Fore, Style
import getopt

from runAMRConvergenceTest import runTest
from mushyLayerRunUtils import getBaseOutputDir, getMatlabBaseCommand

##########################################################################
# 3) Convection in a mushy layer with an initial porous hole
##########################################################################

def testPorousMushyHole(max_time=0.5):

    base_output_dir = getBaseOutputDir()
    matlab_command = getMatlabBaseCommand()

    print(Fore.GREEN + 'Setup tests for porous mushy hole' + Style.RESET_ALL)
    physicalProblem = 'PorousMushyHole'
    dataFolder = os.path.join(base_output_dir, 'PorousMushyHole-t' + str(max_time) )

    Nz_uniform = [16, 32, 64, 128, 256, 512]
    AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': Nz_uniform},
                {'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [16, 32, 64, 128]},
                {'max_level': 2, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [8, 16, 32, 64]},
                {'max_level': 1, 'ref_rat': 4, 'run_types': ['amr'], 'Nzs': [8, 16, 32, 64]}]

    # While testing:
    Nz_uniform = [16, 32, 64]
    AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': Nz_uniform}]

    # Nzs 	  = [16, 32, 64]
    num_procs = [1, 1, 1, 4, 4, 4]  # Needs to be as long as the longest Nzs

    # Setup up the post processing command

    # figureName = os.path.join(dataFolder, 'noFlow.pdf')
    fine_res_folder = 'Uniform-' + physicalProblem + '-' + str(Nz_uniform[-1]) + '--0'
    analysis_command = matlab_command + ' "analyseVariablePorosityTest(\'' + dataFolder + '\', ' + str(
        Nz_uniform[-1]) + ', true, true, \'' + fine_res_folder + '\'); exit;"'

    # Run
    extra_params = {'main.max_time':max_time}
    runTest(dataFolder, physicalProblem, AMRSetup, num_procs, analysis_command, extra_params)

def main(argv):
    # Unknown values:
    max_time = -1

    try:
        opts, args = getopt.getopt(argv, "t:")
    except getopt.GetoptError as err:
        print(str(err))
        print('testPorousMushyHole.py -t<max time>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-t"):
            max_time = float(arg)

    # Only pass in the defined arguments
    if max_time < 0:
        testPorousMushyHole()
    else:
        testPorousMushyHole(max_time)

if __name__ == "__main__":
    main(sys.argv[1:])

