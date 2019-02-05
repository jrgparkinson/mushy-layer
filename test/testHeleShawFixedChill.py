# Run all convergence tests for the methods paper
import os, sys
from colorama import Fore, Style
import getopt

from runAMRConvergenceTest import runTest
from mushyLayerRunUtils import get_base_output_dir, get_matlab_base_command, read_inputs, get_executable_name
from MushyLayerRunSimple import MushyLayerRunSimple
from SlurmTask import  SlurmTask

##########################################################################
# 3) Convection in a mushy layer with an initial porous hole
##########################################################################

def testHeleShawFixedChill(argv):
    mushyLayerBaseDir = os.path.abspath(os.pardir)
    # params_file = mushyLayerBaseDir + '/params/convergenceTest/FixedChill.parameters'
    params_file = mushyLayerBaseDir + '/params/convergenceTest/FixedChill.parameters'

    # Defaults
    defaultParams = read_inputs(params_file)
    extra_params = {'main.max_time': float(defaultParams['main.max_time']),
                    'parameters.rayleighComp': float(defaultParams['parameters.rayleighComp']),
                    'parameters.darcy': float(defaultParams['parameters.darcy']),
                    'parameters.compositionRatio': float(defaultParams['parameters.compositionRatio']),
                    'parameters.nonDimReluctance': float(defaultParams['parameters.nonDimReluctance']),
                    'parameters.prandtl': float(defaultParams['parameters.prandtl'])}

    doAMR = True

    #Pr = 10.0  # fix this for now
    periodic = False

    try:

        opts, args = getopt.getopt(argv, "t:R:D:C:P:N:A")
    except getopt.GetoptError as err:
        print(str(err))
        print(
            'testPorousMushyHole.py -t<max time> -R<Compositional Rayleigh Number> -D<Darcy number> -C<Composition Ratio> -A<do amr> -P<true(1)/false(0) do plot files> -N<non dimensional reluctance>')

        sys.exit(2)

    for opt, arg in opts:
        if opt in "-t":
            extra_params['main.max_time'] = float(arg)
        elif opt in ("-R"):
            extra_params['parameters.rayleighComp'] = float(arg)
        elif opt in ("-D"):
            extra_params['parameters.darcy'] = float(arg)
        elif opt in ("-C"):
            extra_params['parameters.compositionRatio'] = float(arg)
        elif opt in ("-N"):
            extra_params['parameters.nonDimReluctance'] = float(arg)
        elif opt in ("-P"):
            doPlotFiles = bool(int(arg))
            print('Do plot files: ' + str(doPlotFiles) +  ', arg = ' + str(arg))
            if not doPlotFiles:
                extra_params['main.plot_interval'] = -1
        elif opt in ("-A"):
            doAMR = True

    base_output_dir = get_base_output_dir()
    matlab_command = get_matlab_base_command()

    print(Fore.GREEN + 'Setup tests for fixed chill in a Hele-Shaw cell' + Style.RESET_ALL)
    physicalProblem = 'FixedChill'
    folderName = "FixedChill-t%1.1e-Ra%.0e-Da%1.1e-C%1.2f-Rel%1.1e" % (extra_params['main.max_time'], extra_params['parameters.rayleighComp'], extra_params['parameters.darcy'], extra_params['parameters.compositionRatio'], extra_params['parameters.nonDimReluctance'])
    if periodic:
        folderName = folderName + '-periodic'
    dataFolder = os.path.join(base_output_dir, folderName)

    analysis_command = matlab_command + ' "processFixedChill(\'' + dataFolder + '\', [3000, 4800, 17000]); exit;"'


    singleRun = True

    if singleRun:

        num_proc = 1
        defaultParams['concise_run_name'] = folderName

        allowRestarts = False
        s = SlurmTask('', defaultParams['concise_run_name'], '', num_proc)
        ml_run = MushyLayerRunSimple(base_output_dir, num_proc, defaultParams, s, allowRestarts, get_executable_name())
        ml_run.single_run(folderName)


    else:

        Nz_uniform = 256
        Nz_amr_2 = int(float(Nz_uniform) / 2)
        Nz_amr_4 = int(float(Nz_uniform) / 4)
        AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': [Nz_uniform]},
                    {'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [Nz_amr_2]},
                    {'max_level': 2, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [Nz_amr_4]},
                    {'max_level': 1, 'ref_rat': 4, 'run_types': ['amr'], 'Nzs': [Nz_amr_4]}]

        # While testing:
        Nz_uniform = [64]
       # AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': Nz_uniform},
        #            {'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [64]}]
        AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': Nz_uniform}]

        if doAMR:
            AMRSetup.append({'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': Nz_uniform})
            AMRSetup.append({'max_level': 2, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': Nz_uniform})


        num_procs = [1]  # Needs to be as long as the longest Nzs

        # Setup up the post processing command

        # figureName = os.path.join(dataFolder, 'noFlow.pdf')
        #fine_res_folder = 'Uniform-' + physicalProblem + '-' + str(Nz_uniform[-1]) + '--0'

        #analysis_command = '' # No analysis yet. Eventually should collate run times, make some plots, and maybe compute differences between solutions

        runTest(dataFolder, physicalProblem, AMRSetup, num_procs, analysis_command, extra_params, 0, params_file)

if __name__ == "__main__":
    testHeleShawFixedChill(sys.argv[1:])

