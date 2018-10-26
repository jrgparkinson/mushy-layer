# Script to run mushy layer simulations for a variety of parameters
# at fixed resolution

from MushyLayerRun import MushyLayerRun
import os
import numpy
import math
from subprocess import Popen

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

params['main.plot_interval'] = '100'
params['main.check_interval'] = '5'
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
output_base_dir = 'mushyLayer-varyWavelength'

params['main.max_dt_growth'] = '1.005'
#params['main.cfl'] = '0.7'
#params['parameters.rayleighComp'] = 100.0
params['parameters.topEnthalpy'] = 0.0
params['parameters.lewis'] = 100.0

# For testing
params['main.max_step'] = '5'

#--------------------------------------------
# Grid setup
numCellsVertical = 32
horizontal_res = range(32, 40, 4)

run_type = MushyLayerRun.UNIFORM_GRID
gridType = MushyLayerRun.CUSTOM
params['main.max_level'] = '0'

# If choosing fixed grid above, specify grid files here:
base_grid_file = 'grids/grids-viscous-'

# If choosing AMR, refinement criteria needs to change with grid resolution
# else there will be no refinement for finer base grids
# refinement criteria on coarsest grid resolution
params['main.refine_thresh'] = 0.2

#-----------------------------------------
# Parallel stuff
num_proc = 4

grids_file = ''

horizontal_i = 0

for numCellsHorizontal in horizontal_res:

    params['main.num_cells'] = str(numCellsHorizontal) + ' ' + str(numCellsVertical) + ' ' + str(numCellsVertical)

    flags = ''
    if run_type != MushyLayerRun.UNIFORM_GRID and params['main.max_level'] > 0:
        flags = run_type

    # Want to choose a max_box_size to make the parallel computations most
    # efficient
    min_box_size = numpy.sqrt(numCellsHorizontal * numCellsVertical / num_proc)

    # Round this up to the nearest power of 2
    min_box_size = 2**(math.ceil(math.log(min_box_size, 2)))
    params['main.max_box_size'] = str(min_box_size)

    print('Max box size = ' + str(min_box_size) + '\n')

    run_name = 'width' + str(numCellsHorizontal)

    output_dir = output_base_dir
    if flags:
        output_dir = output_base_dir + '-' + flags
    params['main.plot_prefix'] = output_dir + '-' + run_name + '-'
    # check_prefix = plot_prefix
    params['main.chk_prefix'] =  'chk' + params['main.plot_prefix']

    ml_run = MushyLayerRun(output_dir, num_proc, gridType, params, subcycled)

    # Not sure this is necessary any more
    ml_run.custom_grid = [numCellsHorizontal, numCellsVertical,
                          numCellsHorizontal]

    
    ml_run.single_run(run_name)

    # after we've completed a run, create initial conditions for next run
    if len(horizontal_res) > horizontal_i + 1:
        run_dir = ml_run.get_most_recent_directory(run_name)
        final_check_file = ml_run.get_final_output(
            run_dir,  params['main.chk_prefix'])

        final_check_file = os.path.join(run_dir, final_check_file)
        
        if len(final_check_file) == 0:
            print("Error - no checkpoint file!")

        new_file_name = 'chk-init-' + str(horizontal_res[horizontal_i + 1]) + '.2d.hdf5'
        new_restart_file = os.path.join(ml_run.output_dir, new_file_name)

        
        command = ('matlab -nodisplay -nojvm -r \'cd ' + cwd + ' ;ChomboExtend ' +
                   final_check_file + '  ' + new_restart_file +
                   '  2; exit;\'  ')
        print(command)
        
        process = Popen(command, shell=True)
        exit_code = process.wait()

        params['main.restart_file'] = os.path.join(cwd, new_restart_file)

    horizontal_i = horizontal_i + 1

# Copy data to dropbox so we can analyse it
#if ml_run.machine_specific.has_dropbox():
#    ml_run.copy_final_output(res, dropbox_data_folder)

# Copy data from gyre servers to cloud
#if ml_run.machine_specific.is_gyre():
#    ml_run.copy_all_output_to_cloud(res)

# ml_run.continueRun('64-2')
