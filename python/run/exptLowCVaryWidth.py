# Script to run mushy layer simulations for a variety of parameters
# at fixed resolution

from MushyLayerRun import MushyLayerRun
import os
#import numpy
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile

#################################
# These shouldn't change
#################################
cwd = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer'
mushyLayerBaseDir = cwd.replace('/run', '')
params = readInputs(mushyLayerBaseDir + '/params/mushyLayerDarcyLowC.parameters')
os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
subcycled = True

###################################
# Parameters to change
###################################
output_dir = 'mushyLayerLowC-lowerBranch'

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

# If choosing fixed grid above, specify grid file here:
# params['main.gridfile'] = cwd + '/grids/grids-mush-insulating-32-64'

# If choosing AMR, refinement criteria needs to change with grid resolution
# else there will be no refinement for finer base grids
# refinement criteria on coarsest grid resolution
#params['main.refine_thresh'] = 0.3
#params['main.ref_ratio'] = '4 4 4 2'
#params['main.regrid_interval'] = '16 16 16'

# -----------------------------------------
# Parallel stuff
num_proc = 4

# Make output dir name 
if params['main.periodic_bc'] == '0 0 1':
	output_dir = output_dir + '-insulating'
	    
	# To force channels to form at edge of domain
	params['main.initialPerturbation'] = '0.1'
else:
	output_dir = output_dir + '-periodic'

	# This appends "fixedGrid" or "adaptive" if appropriate
if run_type != MushyLayerRun.UNIFORM_GRID and int(params['main.max_level']) > 0:
	output_dir = output_dir + '-' + run_type

# Construct name for run
params['parameters.rayleighComp'] = 250
maxStep = 99999

params['main.checkpoint_interval'] = '100'
params['main.plot_interval'] = '100'

params['main.max_grid_size'] = '64'
params['main.max_dt_growth'] = '1.02'

restartFile = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/run/mushyLayerLowC-lowerBranch-insulating/chkmushyLayerLowC-insulating-CR1.25RaC250Le200ChiCubedPermeabilitypts16-003453.2d.hdf5'
params['main.restart_perturbation'] = '0.001'

#Nx_arr = [16, 20, 24, 28, 32, 40, 48, 56, 64] # 122
Nx_arr = [16, 20]
Nx_i = 0

for Nx in Nx_arr:
	# for testing
	params['main.max_step'] = str(maxStep)

	# See if we can get away with this
	params['main.num_cells'] = str(Nx) + ' 128 64'
	params['main.domain_length'] = 0.2*float(Nx)/float(128)

	#params['main.restart_file'] = mushyLayerBaseDir + '/run/mushyLayerLowC-periodic/chkmushyLayerLowC-periodic-CR1.25RaC' + str(params['parameters.rayleighComp']) + 'Le200ChiCubedPermeabilitypts'+str(Nx)+'x128-centred.2d.hdf5'
	#params['main.restart_file'] = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/run/mushyLayerLowC-periodic/CR1.25RaC250Le200ChiCubedPermeabilitypts116-0/chkmushyLayerLowC-periodic-CR1.25RaC250Le200ChiCubedPermeabilitypts116-009650.2d.hdf5'

        if restartFile:
            params['main.restart_file'] = restartFile


	run_name = constructRunName(params)
	#run_name = run_name + '-Nx' + str(Nx)
	params['main.plot_prefix'] = output_dir + '-' + run_name + '-'
	params['main.chk_prefix'] =  'chk' + params['main.plot_prefix']

	# Do run
	ml_run = MushyLayerRun(output_dir, num_proc, params, subcycled)
	ml_run.single_run(run_name)


        # Process output file to give inputs for next run
        if Nx_i+1 >= len(Nx_arr):
            continue
        nextNx = Nx_arr[Nx_i+1]
        processParams = readInputs(mushyLayerBaseDir + '/setupNewRun/inputs')
        prevChkFile =  ml_run.finalChkFile()
        bits = prevChkFile.split('/')
        base_output_dir = '/'.join(bits[:-1])
        runInputsFile = base_output_dir + '/inputs'
        processParams['inFile'] = prevChkFile
        processParams['run_inputs'] = runInputsFile
        processParams['deltaNxLeft'] = str(nextNx - Nx)
        processParams['deltaNxRight'] = '0'
        processParams['box_size'] = params['main.max_grid_size']

        this_dir = '/'.join(bits[:-2])
        restartFile = this_dir + '/' +  str(nextNx) + "-restart.2d.hdf5"

        processParams['outFile'] =restartFile 

        inputsFile = 'processInputsTemp'
        writeParamsFile(inputsFile, processParams)
        processProgram = mushyLayerBaseDir + '/setupNewRun/setupnewrun2d.Linux.64.mpiCC.gfortran.DEBUG.MPI.ex'
        os.system(processProgram + ' ' + inputsFile);

        Nx_i = Nx_i + 1


	
		

