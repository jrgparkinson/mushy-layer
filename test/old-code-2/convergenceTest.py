# Script to run a convergence test for a specified param file

from MushyLayerRun import MushyLayerRun
import os
import numpy
import math

# Define different types of convergence tests
# Uniform - just a series of uniform grids
# Fixed Grid - a series of fixed grids
# AMR - AMR with specified number of levels
TEST_UNIFORM = 0
TEST_FIXED_GRID = 1
TEST_AMR = 2

params = {}

###################################
# Parameters to change
###################################
params['params_file'] = 'params/hrl.parameters'
output_base_dir = 'subcycled-hrl'

#--------------------------------------------
## Grid setup
grid_res = [ 2**j for j in range(5,10) ] # 2^4 = 16, 2^10 = 1024
run_type = MushyLayerRun.UNIFORM_GRID
gridType = MushyLayerRun.SQUARE
max_level = '0'

# If choosing fixed grid above, specify grid files here:
base_grid_file = 'grids/grids-viscous-'

# If choosing AMR, refinement criteria needs to change with grid resolution 
# else there will be no refinement for finer base grids
coarse_refinement_criteria = 0.2 # refinement criteria on coarsest grid resolution

#-----------------------------------------
## Parallel stuff
num_proc = 1

#################################
# These shouldn't change
#################################
os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
subcycled = True

plot_interval = 100
check_interval = 100
params['main.verbosity'] = '1'
params['main.max_time'] = '99999.0'
params['main.max_step'] = '99999'
params['main.steady_state'] = '1e-3'
params['main.enforceAnalyticSoln'] = 'true' # this will speed things up a bit

flags = ''
if run_type != MushyLayerRun.UNIFORM_GRID and max_level > 0:
    flags = run_type

for res in grid_res:
    grids_file = ''

    # Sort out a few things that need to change for each run
    if run_type == MushyLayerRun.FIXED_GRID and max_level > 0:
        params['main.gridfile'] = (base_grid_file  + str(res) + "-"
                              + str(res*2))
    if run_type == MushyLayerRun.ADAPTIVE_GRID and max_level > 0:
        params['main.refine_thresh'] = str(coarse_refinement_criteria * grid_res[0]/res)

    # Want to choose a max_box_size to make the parallel computations most efficient
    min_box_size = numpy.sqrt(res*res/num_proc)

    # Round this up to the nearest power of 2
    min_box_size = 2**(math.ceil(math.log(min_box_size, 2)))
    params['main.max_box_size'] = str(min_box_size)

    output_dir = output_base_dir
    if flags:
        output_dir = output_base_dir+'-' + flags
    plot_prefix = output_dir + '-'
    check_prefix = plot_prefix

    ml_run = MushyLayerRun(output_dir,  plot_prefix,  check_prefix,
                           max_level,  num_proc,  plot_interval,
                           check_interval, gridType, params, subcycled)

    ml_run.run([res])

    # Copy data to dropbox so we can analyse it
    #if ml_run.machine_specific.has_dropbox():
    #    ml_run.copy_final_output(res, dropbox_data_folder)

    # Copy data from gyre servers to cloud
    if ml_run.machine_specific.is_gyre():
        ml_run.copy_all_output_to_cloud(res)

    # ml_run.continueRun('64-2')
