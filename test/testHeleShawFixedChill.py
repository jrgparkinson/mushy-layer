# Run all convergence tests for the methods paper
import os, sys
from colorama import Fore, Style
import getopt

from runAMRConvergenceTest import runTest
from mushyLayerRunUtils import getBaseOutputDir, getMatlabBaseCommand, readInputs

##########################################################################
# 3) Convection in a mushy layer with an initial porous hole
##########################################################################

def testHeleShawFixedChill(argv):
    mushyLayerBaseDir = os.path.abspath(os.pardir)
    # params_file = mushyLayerBaseDir + '/params/convergenceTest/FixedChill.parameters'
    params_file = mushyLayerBaseDir + '/params/convergenceTest/fixedChillNew.parameters'

    # Defaults
    defaultParams = readInputs(params_file)
    extra_params = {}

    extra_params['main.max_time']=float(defaultParams['main.max_time'])
    extra_params['parameters.rayleighComp']=float(defaultParams['parameters.rayleighComp'])
    extra_params['parameters.darcy']=float(defaultParams['parameters.darcy'])
    extra_params['parameters.compositionRatio']=float(defaultParams['parameters.compositionRatio'])
    extra_params['parameters.nonDimReluctance'] = float(defaultParams['parameters.nonDimReluctance'])
    extra_params['parameters.prandtl'] = float(defaultParams['parameters.prandtl'])

    doAMR = False

    #Pr = 10.0  # fix this for now
    periodic = True


    try:
        opts, args = getopt.getopt(argv, "t:R:D:C:P:A")
    except getopt.GetoptError as err:
        print(str(err))
        print(
            'testPorousMushyHole.py -t<max time> -R<Compositional Rayleigh Number> -D<Darcy number> -C<Composition Ratio> -A<do amr> -P<true(1)/false(0) do plot files>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-t"):
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

    base_output_dir = getBaseOutputDir()
    matlab_command = getMatlabBaseCommand()

    print(Fore.GREEN + 'Setup tests for fixed chill in a Hele-Shaw cell' + Style.RESET_ALL)
    physicalProblem = 'FixedChill'
    folderName = "FixedChill-t%1.1e-Ra%.0e-Da%1.1e-C%1.2f-Rel%1.1e" % (extra_params['main.max_time'], extra_params['parameters.rayleighComp'], extra_params['parameters.darcy'], extra_params['parameters.compositionRatio'], extra_params['parameters.nonDimReluctance'])
    if periodic:
        folderName = folderName + '-periodic'
    dataFolder = os.path.join(base_output_dir, folderName)



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

    # Nzs 	  = [16, 32, 64]
    num_procs = [4]  # Needs to be as long as the longest Nzs

    # Setup up the post processing command

    # figureName = os.path.join(dataFolder, 'noFlow.pdf')
    #fine_res_folder = 'Uniform-' + physicalProblem + '-' + str(Nz_uniform[-1]) + '--0'
    analysis_command = matlab_command + ' "processFixedChill(\'' + dataFolder + '\'); exit;"'
    #analysis_command = '' # No analysis yet. Eventually should collate run times, make some plots, and maybe compute differences between solutions


    # Run
    #extra_params = {'main.max_time':max_time, 'parameters.rayleighComp':Ra, 'parameters.darcy': Da, 'parameters.compositionRatio':C }
    extra_params['main.max_dt'] = 0.1*extra_params['parameters.prandtl'] / (extra_params['parameters.darcy']*extra_params['parameters.rayleighComp'])

    if 'main.plot_interval' in extra_params and int(extra_params['main.plot_interval']) > 0: 
    	extra_params['main.plot_interval'] = int(1e-4 / float(extra_params['parameters.darcy']))
    	extra_params['main.checkpoint_interval'] = extra_params['main.plot_interval']

    extra_params['main.useAccelDt'] = 0

    extra_params['regrid.plume_salinity']=-1.0
    extra_params['regrid.plume_vel']= 0.1*extra_params['parameters.darcy']*extra_params['parameters.rayleighComp']*extra_params['parameters.rayleighComp']*extra_params['parameters.prandtl']
    if periodic:
        extra_params['main.periodic_bc'] = '1 0 0'

    runTest(dataFolder, physicalProblem, AMRSetup, num_procs, analysis_command, extra_params, 0, params_file)

# def main(argv):
    # Only pass in the defined arguments
    # if max_time < 0:
    #     testHeleShawFixedChill()
    # elif Ra < 0:
    #     testHeleShawFixedChill(max_time)
    # elif Da < 0:
    #     testHeleShawFixedChill(max_time, Ra)
    # elif C < 0:
    #     testHeleShawFixedChill(max_time, Ra, Da)
    # else:
    #     testHeleShawFixedChill(max_time, Ra, Da, C)

if __name__ == "__main__":
    testHeleShawFixedChill(sys.argv[1:])

