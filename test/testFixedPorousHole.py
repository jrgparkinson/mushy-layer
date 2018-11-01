# Run all convergence tests for the methods paper
import os, sys
from colorama import Fore, Style

from runAMRConvergenceTest import runTest
from mushyLayerRunUtils import getBaseOutputDir, getMatlabBaseCommand

##########################################################################
# 3) Convection in a fixed porous medium with variable porosity
##########################################################################

def testFixedPorousHole():

    base_output_dir = getBaseOutputDir()
    matlab_command = getMatlabBaseCommand()

    print(Fore.GREEN + 'Setup tests for (fixed) porous hole' + Style.RESET_ALL)
    physicalProblem = 'DBVariablePorosity'
    dataFolder = os.path.join(base_output_dir, 'FixedPorousHole')

    Nz_uniform = [16, 32, 64, 128, 256, 512]
    AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': Nz_uniform},
                {'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [16, 32, 64, 128]},
                {'max_level': 2, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [8, 16, 32, 64]},
                {'max_level': 1, 'ref_rat': 4, 'run_types': ['amr'], 'Nzs': [8, 16, 32, 64]}]

    # While testing:


    # Nzs 	  = [16, 32, 64]
    num_procs = [1, 1, 1, 4, 4, 4]  # Needs to be as long as the longest Nzs

    # Setup up the post processing command

    # figureName = os.path.join(dataFolder, 'noFlow.pdf')
    uniform_prefix = 'Uniform-DBVariablePorosity-'
    analysis_command = matlab_command + ' "analyseVariablePorosityTest(\'' + dataFolder + '\', ' + ',['.join([str(a) for a in Nz_uniform]) + '], true, true, \'' + uniform_prefix + '\'); exit;"'

    # Run
    extra_params = {}
    runTest(dataFolder, physicalProblem, AMRSetup, num_procs, analysis_command, extra_params)

def main(argv):
    testFixedPorousHole()

if __name__ == "__main__":
    main(sys.argv[1:])

