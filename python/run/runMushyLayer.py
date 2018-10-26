# Script to run mushy layer simulations for a variety of parameters
# at fixed resolution

from MushyLayerRun import MushyLayerRun
import os
import numpy
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName

#################################
# These shouldn't change
#################################
params = {}
cwd = os.getcwd()
mushyLayerBaseDir = cwd.replace('/run', '')
os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
subcycled = True

###################################
# Parameters to change
###################################

params['main.params_file'] = mushyLayerBaseDir + '/params/mushyLayerDarcy.parameters'
#params['main.restart_file'] = cwd + '/' + output_base_dir + '/.2d.hdf5'
#params['main.restart_file'] = cwd + '/' + output_base_dir + '/chk-init-44.2d.hdf5'

# Restart
#params['main.restart_file'] = cwd + '/mushyLayer-periodic/128restart.2d.hdf5'
#params['main.restart_perturbation'] = 0.02


# --------------------------------------------
# Grid setup
grid_res = 16
ncells = [grid_res grid_res grid_res]
params['main.num_cells'] = str(ncells).join(' ')

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
num_proc = 1

# Construct name for run
output_dir = 'mushyLayer'
if params['main.periodic_bc'] == '0 0 1':
	output_dir = output_dir + '-insulating'
	    
	# To force channels to form at edge of domain
	params['main.initialPerturbation'] = '0.1'
else:
	output_dir = output_dir + '-periodic'

	# This appends "fixedGrid" or "adaptive" if appropriate
if run_type != MushyLayerRun.UNIFORM_GRID and int(params['main.max_level']) > 0:
	output_dir = output_dir + '-' + run_type

run_name = constructRunName(params)
params['main.plot_prefix'] = output_dir + '-' + run_name + '-'
params['main.chk_prefix'] =  'chk' + params['main.plot_prefix']

# Do run
ml_run = MushyLayerRun(output_dir, num_proc, params, subcycled)
ml_run.single_run(run_name)

