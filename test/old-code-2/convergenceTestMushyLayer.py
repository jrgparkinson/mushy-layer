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
#################################
# These shouldn't change
#################################
os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
subcycled = True

params['main.plot_interval'] = '10'
params['main.checkpoint_interval'] = '10'
params['main.verbosity'] = '2'
params['main.max_time'] = '99999.0'
params['main.max_step'] = '99999'
params['main.steady_state'] = '5e-2'
params['main.enforceAnalyticSoln'] = 'true'  # this will speed things up a bit

##################################
# Overwriting default options:
##################################
# based on other runs, after t=20 it's more or less steady state
params['main.max_time'] = '5.0'

###################################
# Parameters to change
###################################
cwd = os.getcwd()
params['main.params_file'] = cwd + '/params/mushyLayer2.parameters'
output_base_dir = 'mushyLayer-insulating'
params['main.restart_file'] = cwd + \
   '/mushyLayer-insulating/chk-32pts-finished.2d.hdf5'
params['main.max_dt_growth'] = '1.005'
# params['main.cfl'] = '0.7'
# params['parameters.rayleighComp'] = 100.0
params['parameters.topEnthalpy'] = 0.0
params['parameters.lewis'] = 500.0

# --------------------------------------------
# Grid setup
# grid_res = [ 2**j for j in range(6,7) ] # 2^4 = 16, 2^10 = 1024
grid_res = [32]
run_type = MushyLayerRun.UNIFORM_GRID
gridType = MushyLayerRun.SQUARE
params['main.max_level'] = '0'

# If choosing fixed grid above, specify grid files here:
base_grid_file = 'grids/grids-viscous-'

# If choosing AMR, refinement criteria needs to change with grid resolution
# else there will be no refinement for finer base grids
# refinement criteria on coarsest grid resolution
coarse_refinement_criteria = 0.2

# -----------------------------------------
# Parallel stuff
num_proc = 1

flags = ''
if run_type != MushyLayerRun.UNIFORM_GRID and params['main.max_level'] > 0:
    flags = run_type

for res in grid_res:
    grids_file = ''

    params['main.num_cells'] = str(res) + ' ' + str(res) + ' ' + str(res)

    # Sort out a few things that need to change for each run
    if run_type == MushyLayerRun.FIXED_GRID and params['main.max_level'] > 0:
        params['main.gridfile'] = (base_grid_file + str(res) + "-"
                                   + str(res * 2))
    if run_type == MushyLayerRun.ADAPTIVE_GRID and params['main.max_level'] > 0:
        params['main.refine_thresh'] = str(
            coarse_refinement_criteria * grid_res[0] / res)

    # Want to choose a max_box_size to make the parallel computations most
    # efficient
    min_box_size = numpy.sqrt(res * res / num_proc)

    # Round this up to the nearest power of 2
    min_box_size = 2**(math.ceil(math.log(min_box_size, 2)))
    params['main.max_box_size'] = str(min_box_size)

    print('Max box size = ' + str(min_box_size) + '\n')

    output_dir = output_base_dir
    if flags:
        output_dir = output_base_dir + '-' + flags
    params['main.plot_prefix'] = output_dir + '-'
    params['main.chk_prefix'] = 'chk-' + params['main.plot_prefix']

    ml_run = MushyLayerRun(output_dir, num_proc, gridType, params, subcycled)

    run_name = str(res) + 'pts-Le' + str(params['parameters.lewis'])

    ml_run.single_run(run_name)


