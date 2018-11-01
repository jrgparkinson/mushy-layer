# Run all convergence tests for the methods paper
import os, sys
from colorama import Fore, Style

from runAMRConvergenceTest import runTest
from mushyLayerRunUtils import getBaseOutputDir, getMatlabBaseCommand

##########################################################################
# 3) Convection in a mushy layer with an initial porous hole
##########################################################################

def testPorousMushyHole():

    base_output_dir = getBaseOutputDir()
    matlab_command = getMatlabBaseCommand()

    print(Fore.GREEN + 'Setup tests for porous mushy hole' + Style.RESET_ALL)
    physicalProblem = 'PorousMushyHole'
    dataFolder = os.path.join(base_output_dir, 'PorousMushyHole')

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
    extra_params = {}
    runTest(dataFolder, physicalProblem, AMRSetup, num_procs, analysis_command, extra_params)

def main(argv):
    testPorousMushyHole()

if __name__ == "__main__":
    main(sys.argv[1:])

