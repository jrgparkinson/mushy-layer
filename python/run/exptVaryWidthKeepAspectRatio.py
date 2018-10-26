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
output_dir = 'mushyLayerLowC-upperBranch'

params['main.multiCompSolve'] = 'true'
params['picardSolver.max_iter'] = '1'

# --------------------------------------------
# Grid setup

# Quite possibly oscillating modes of convection so should be periodic
params['main.periodic_bc'] = '0 0 1'

run_type = MushyLayerRun.UNIFORM_GRID
params['main.max_level'] = '0'


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
params['main.max_step'] = '999999'

params['main.checkpoint_interval'] = '100'
params['main.plot_interval'] = '100'

params['main.max_grid_size'] = '32'
params['main.max_dt_growth'] = '1.02'
params['main.steady_state'] = '1e-4'

restartFile = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/run/mushyLayerLowC-upperBranch-insulating/chkmushyLayerLowC-insulating-CR1.25RaC250Le200ChiCubedPermeabilitypts32-037322.2d.hdf5'
params['main.restart_perturbation'] = '0.0'

params['main.num_cells'] = '32 128 64'

#Nx_arr = [16, 20, 24, 28, 32, 40, 48, 56, 64] # 122
firstWidth = 0.05
increment = -0.0005
numSteps = 20
width_arr = [(firstWidth + increment*i) for i in range(1, numSteps)]
i = 0

for width in width_arr:
	params['main.domain_length'] = str(width)

        if restartFile:
            params['main.restart_file'] = restartFile

	run_name = constructRunName(params)
	#run_name = run_name + '-Nx' + str(Nx)
	params['main.plot_prefix'] = output_dir + '-' + run_name + '-'
	params['main.chk_prefix'] =  'chk' + params['main.plot_prefix']

	# Do run
	ml_run = MushyLayerRun(output_dir, num_proc, params, subcycled)
	ml_run.single_run(run_name)

        # Get restart file for next run
        if i >= len(width_arr):
            continue
        
        restartFile =  ml_run.finalChkFile()
        
        i = i + 1


