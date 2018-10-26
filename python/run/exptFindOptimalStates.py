# Script to find optimal states, varying the domain
# width to find the optimal flux for fixed Ra, then increasing Ra and 
# varying the width again

from MushyLayerRun import MushyLayerRun
import os
#import numpy
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile

CHKFILE = 'chkfile'
WIDTH = 'width'
FLUX = 'flux'
#################################
# These shouldn't change
#################################
#cwd = '/home/parkinsonj/convection-in-sea-ice/MushyLayer'
cwd = os.getcwd()
mushyLayerBaseDir = cwd.replace('/run', '')
params = readInputs(mushyLayerBaseDir + '/params/mushyLayerDarcyLowC.parameters')
os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
subcycled = True

###################################
# Parameters to change
###################################
#output_dir = 'mushyLayerLowC-lowerBranch'
output_dir = 'mushyLayerLowC-upperBranch'

#params['main.params_file'] = mushyLayerBaseDir + '/params/mushyLayerDarcyLowC.parameters'
#params['main.restart_file'] = cwd + '/' + output_base_dir + '/.2d.hdf5'
#params['main.restart_file'] = cwd + '/' + output_base_dir + '/chk-init-44.2d.hdf5'

# Restart
#params['main.restart_file'] = cwd + '/mushyLayer-periodic/128restart.2d.hdf5'
#params['main.restart_perturbation'] = 0.02

params['main.multiCompSolve'] = 'true'
params['picardSolver.max_iter'] = '1'


# --------------------------------------------
# Grid setup

# Quite possibly oscillating modes of convection so should be periodic
params['main.periodic_bc'] = '0 0 1'

run_type = MushyLayerRun.UNIFORM_GRID
params['main.max_level'] = '0'

params['main.steady_state']='1e-4'

# -----------------------------------------
# Parallel stuff
num_proc = 8

# Make output dir name 
if params['main.periodic_bc'] == '0 0 1':
	output_dir = output_dir + '-insulating'
	    
	# To force channels to form at edge of domain
	# params['main.initialPerturbation'] = '0.1'
else:
	output_dir = output_dir + '-periodic'

	# This appends "fixedGrid" or "adaptive" if appropriate
if run_type != MushyLayerRun.UNIFORM_GRID and int(params['main.max_level']) > 0:
	output_dir = output_dir + '-' + run_type

# Construct name for run

params['parameters.compositionRatio'] = 1.1
params['main.max_step'] = '9999999'
maxStep = 170138


params['main.checkpoint_interval'] = '2000'
params['main.plot_interval'] = '2000'

params['main.max_grid_size'] = '32'
params['main.max_dt_growth'] = '1.02'

prevChkFile = '/home/parkinsonj/convection-in-sea-ice/MushyLayer/run/mushyLayerLowC-upperBranch-insulating/CR1.25RaC400Le200ChiCubedPermeabilitypts36-0/chkmushyLayerLowC-upperBranch-insulating-CR1.25RaC400Le200ChiCubedPermeabilitypts36-188663.2d.hdf5'

# +1 or -1, to denote whether the channel is on the left(-1) or right(+1)
channelSide = -1

params['main.restart_perturbation'] = '0.0'

# Initial number of horizontal grid cells
Nx_old = 32
Nx_new = Nx_old
# how many grid cells to change by. 2 is the smallest allowed by multigrid
deltaNx = 2

Ra_arr = [600, 700, 800, 900, 1000]
Ra_arr = [350, 300, 350, 200, 150, 100]
#Ra_arr = [400, 500, 600]

for Ra in Ra_arr:
    # need to keep track of fluxes and checkpoint files for this Ra to see if we find the maximum
    completedRuns = []

    params['parameters.rayleighComp'] = Ra
    decreaseWidth = (Ra_arr[1] > Ra_arr[0]) # initially try decreasing the width if we're increasing the wavelength
    
    foundMaxFlux = False
    while foundMaxFlux == False: 
        maxStep = maxStep + 2
        #params['main.max_step'] = str(maxStep) #testing

        # Create restart file for this run
        processParams = readInputs(mushyLayerBaseDir + '/setupNewRun/inputs')
        
        # Work out what side the previous checkpoint file channel is on
        chkFile = ChkFile(prevChkFile)
        [numChannels, relPos] = chkFile.channelProperties()
        
        if numChannels == 1:
            
            if relPos[0] < 0.5:
                channelSide = -1
                sideStr = 'left'
            else:
                channelSide = +1
                sideStr = 'right'
                
            
            print('Determined that the channel is on the ' + sideStr + ' side. Adding/subtracting cells from the opposite side.')
                
        elif numChannels > 1:
            print('More than one channel - stop simulations')
            quit() # stop if we've started creating multiple channels
            

        bits = prevChkFile.split('/')
        base_output_dir = '/'.join(bits[:-1])
        runInputsFile = base_output_dir + '/inputs'
        processParams['inFile'] = prevChkFile
        processParams['run_inputs'] = runInputsFile
        # If the channel is on the left (channelSide < 0) then remove cells from the right
        if channelSide < 0:
            processParams['deltaNxRight'] = str(Nx_new-Nx_old)
            processParams['deltaNxLeft'] = '0'
        else:
            processParams['deltaNxLeft'] = str(Nx_new-Nx_old)
            processParams['deltaNxRight'] = '0'
        processParams['box_size'] = params['main.max_grid_size']
        
        this_dir = '/'.join(bits[:-2])
        restartFile = this_dir + '/' +  str(Nx_new) + "-restart.2d.hdf5"

        processParams['outFile'] = restartFile 

        inputsFile = 'processInputsTemp'
        writeParamsFile(inputsFile, processParams)
        processProgram = mushyLayerBaseDir + '/setupNewRun/setupnewrun2d.Linux.64.mpiCC.gfortran.DEBUG.MPI.ex'
        os.system(processProgram + ' ' + inputsFile);


	# Set params for this run
	params['main.num_cells'] = str(Nx_new) + ' 256 64'
	params['main.domain_length'] = 0.2*float(Nx_new)/float(256)

        if restartFile:
            params['main.restart_file'] = restartFile


	run_name = constructRunName(params)
	#run_name = run_name + '-Nx' + str(Nx)
	params['main.plot_prefix'] = output_dir + '-' + run_name + '-'
	params['main.chk_prefix'] =  'chk' + params['main.plot_prefix']

	# Do run
	ml_run = MushyLayerRun(output_dir, num_proc, params, subcycled)
	ml_run.single_run(run_name)

        # Get the flux from this run
        prevChkFile =  ml_run.finalChkFile()
        flux = ml_run.getFinalFlux() 
        thisRun = {WIDTH: Nx_new, FLUX: flux, CHKFILE:prevChkFile}
        completedRuns.append(thisRun)

        # Get ready for next run
        # This ensures that we modify the checkpoint file before the next run by adding/subtracting some cells
        Nx_old = Nx_new
        if decreaseWidth:
            Nx_new = Nx_new - deltaNx
        else:
            Nx_new = Nx_new + deltaNx

        # Work out if we've found the maximum
        # and if so, specify which checkpoint file we should restart from for the next run
         
        # First, need to sort completedRuns by width, as they may not have been all run in the same order
        completedRuns = sorted(completedRuns, key=lambda k: k[WIDTH])
 
        if len(completedRuns) > 2:
            for i in range(1, len(completedRuns)-1):
                if completedRuns[i][FLUX] > completedRuns[i-1][FLUX] and completedRuns[i][FLUX] > completedRuns[i+1][FLUX]:
                    foundMaxFlux = True
                    if decreaseWidth:
                        prevChkFile = completedRuns[i+1][CHKFILE]
                        Nx_new = completedRuns[i+1][WIDTH]
                    else:
                        prevChkFile = completedRuns[i-1][CHKFILE]
                        Nx_new = completedRuns[i-1][WIDTH]

                    # Make sure we don't do any processing of the checkpoint file - it's already the correct width
                    Nx_old = Nx_new

                    # Produce some logging:
                    print('Found a max flux, after runs: ')
                    print(completedRuns)
                    print('Max flux was: ' + completedRuns[i][FLUX])
                    print('Restarting with new Ra and Nx = ' + Nx_new)

        # If we haven't found a maximum, check we're searching in the right direction
        # remember that completedRuns has been sorted by width
        if foundMaxFlux == False and len(completedRuns) > 1:

            # Step through simulations in the same direction we've been going
            if decreaseWidth:
                r = range(len(completedRuns), 1)
                sign = -1
            else:
                r = range(1, len(completedRuns))
                sign = 1

            
            for i in r:
                prev_i = i - sign
                if completedRuns[i][FLUX] < completedRuns[prev_i][FLUX]:
                    # If this flux is less than the previous one, change direction!
                    decreaseWidth = (decreaseWidth == False)
                    
                    # Also, restart from the maximum flux found so far
                    Nx_old = completedRuns[prev_i][WIDTH]
                    Nx_new = Nx_old - sign*deltaNx
                    prevChkFile = completedRuns[prev_i][CHKFILE]


                    # Finally, break out of this for loop
                    break
                    
	

