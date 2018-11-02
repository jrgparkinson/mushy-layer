# Run all convergence tests for the methods paper
import os, sys
from colorama import Fore, Style
import getopt

from runAMRConvergenceTest import runTest
from mushyLayerRunUtils import getBaseOutputDir, getMatlabBaseCommand

##########################################################################
# 3) Convection in a mushy layer with an initial porous hole
##########################################################################

def testHeleShawFixedChill(max_time=1e6, Ra=4e3, Da=2e-3, C=2.0):

    base_output_dir = getBaseOutputDir()
    matlab_command = getMatlabBaseCommand()

    print(Fore.GREEN + 'Setup tests for fixed chill in a Hele-Shaw cell' + Style.RESET_ALL)
    physicalProblem = 'FixedChill'
    folderName = "FixedChill-t%1.1f-Ra%.0e-Da%1.1e-C%1.2f" % (max_time, Ra, Da, C)
    dataFolder = os.path.join(base_output_dir, folderName)

    Nz_uniform = 256
    Nz_amr_2 = int(float(Nz_uniform) / 2)
    Nz_amr_4 = int(float(Nz_uniform) / 4)
    AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': [Nz_uniform]},
                {'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [Nz_amr_2]},
                {'max_level': 2, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [Nz_amr_4]},
                {'max_level': 1, 'ref_rat': 4, 'run_types': ['amr'], 'Nzs': [Nz_amr_4]}]

    # While testing:
    Nz_uniform = [128]
    AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': Nz_uniform}]

    # Nzs 	  = [16, 32, 64]
    num_procs = [1]  # Needs to be as long as the longest Nzs

    # Setup up the post processing command

    # figureName = os.path.join(dataFolder, 'noFlow.pdf')
    #fine_res_folder = 'Uniform-' + physicalProblem + '-' + str(Nz_uniform[-1]) + '--0'
    analysis_command = matlab_command + ' "processFixedChill(\'' + dataFolder + '\'); exit;"'
    #analysis_command = '' # No analysis yet. Eventually should collate run times, make some plots, and maybe compute differences between solutions

    # Run
    extra_params = {'main.max_time':max_time, 'parameters.rayleighComp':Ra, 'parameters.darcy': Da, 'parameters.compositionRatio':C }
    extra_params['main.max_dt'] = 10.0 / (Da*Ra)
    runTest(dataFolder, physicalProblem, AMRSetup, num_procs, analysis_command, extra_params)

def main(argv):
    # Unknown values:
    max_time = -1
    Ra = -1
    Da = -1
    C = -1

    try:
        opts, args = getopt.getopt(argv, "t:R:D:C:")
    except getopt.GetoptError as err:
        print(str(err))
        print('testPorousMushyHole.py -t<max time> -R<Compositional Rayleigh Number> -D<Darcy number> -C<Composition Ratio>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-t"):
            max_time = float(arg)
        elif opt in ("-R"):
            Ra = float(arg)
        elif opt in ("-D"):
            Da = float(arg)
        elif opt in ("-C"):
            C = float(arg)

    # Only pass in the defined arguments
    if max_time < 0:
        textFixedChillHeleShaw()
    elif Ra < 0:
        textFixedChillHeleShaw(max_time)
    elif Da < 0:
        textFixedChillHeleShaw(max_time, Ra)
    elif C < 0:
        textFixedChillHeleShaw(max_time, Ra, Da)
    else:
        textFixedChillHeleShaw(max_time, Ra, Da, C)

if __name__ == "__main__":
    main(sys.argv[1:])

