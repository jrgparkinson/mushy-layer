# Script to run mushy layer simulations for a variety of parameters
# at fixed resolution

from MushyLayerRun import MushyLayerRun
import os
#import numpy
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs

#################################
# These shouldn't change
#################################
cwd = os.getcwd()
mushyLayerBaseDir = cwd.replace('/run', '')
params = readInputs(mushyLayerBaseDir + '/params/mushyLayerDarcyLowC.parameters')
os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
subcycled = True

###################################
# Parameters to change
###################################
output_dir = 'mushyLayerLowC-wideDomain'

#params['main.params_file'] = mushyLayerBaseDir + '/params/mushyLayerDarcyLowC.parameters'
#params['main.restart_file'] = cwd + '/' + output_base_dir + '/.2d.hdf5'
#params['main.restart_file'] = cwd + '/' + output_base_dir + '/chk-init-44.2d.hdf5'

# Restart
#params['main.restart_file'] = cwd + '/mushyLayer-periodic/128restart.2d.hdf5'
#params['main.restart_perturbation'] = 0.02


# --------------------------------------------
# Grid setup
# See if we can get away with this
params['main.num_cells'] = '512 128 64'

# Quite possibly oscillating modes of convection so should be periodic
params['main.periodic_bc'] = '1 0 1'
params['main.domain_length'] = '0.8'

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
Ra_arr = [300]
prevFinalChkFile = ''
stepInterval = 3000
maxStep = 40000

params['main.checkpoint_interval'] = '100'
params['main.plot_interval'] = '100'

#params['main.restart_file'] = mushyLayerBaseDir + '/run/mushyLayerLowC-periodic/chkmushyLayerLowC-periodic-CR1.25RaC200Le200ChiCubedPermeabilitypts64-000750.2d.hdf5'

for Ra in Ra_arr:
	# for testing
	params['main.max_step'] = str(maxStep)

	params['parameters.rayleighComp'] = Ra
	if prevFinalChkFile:
		params['main.restart_file'] = prevFinalChkFile

	run_name = constructRunName(params)
	params['main.plot_prefix'] = output_dir + '-' + run_name + '-'
	params['main.chk_prefix'] =  'chk' + params['main.plot_prefix']

	# Do run
	ml_run = MushyLayerRun(output_dir, num_proc, params, subcycled)
	ml_run.single_run(run_name)

	prevFinalChkFile = ml_run.finalChkFile()
	if prevFinalChkFile:
		name = prevFinalChkFile.replace('.2d.hdf5', '')
		number = name[-5:]
		finalStep = int(number)
		maxStep = finalStep + stepInterval

		

