from MushyLayerRun import MushyLayerRun
import os
import numpy

os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"

params = {}

aspect_ratio = 6 # 6.05 stick with integer - shouldn't affect answer too much
Ra_arr = [2000, 2500, 3000, 5000, 10000, 20000, 30000, 50000]
Pr = 0.7

Ra_arr = [5000, 10000, 20000, 30000, 50000]
Ra_arr = [30000, 50000]
Ra_arr = [4000]
params['main.cfl'] = '0.5' # wondering if a smaller cfl will fix issues with high Ra convection
params['main.initial_cfl'] = '0.1'
params['main.time_integration_order'] = '1' # try more accurate integration in time?

params['main.max_step'] = '150' # make sure we conclude the run in a sensible period of time

#Ra_arr = [2000]

# These are the options to specify
run_type = MushyLayerRun.ADAPTIVE_GRID
max_level = 1
num_proc = 1
plot_interval = 100
check_interval = 100
params['main.verbosity'] = '3'
gridType = MushyLayerRun.CUSTOM
grid_res = 8

subcycled = True

# These parameter should be the same for every run
params['main.params_file'] = 'params/rayleighBenard.parameters'

# for testing
#params['main.restart_file'] = 'rayleighBenard-Ra2000/8-4/chk000500.2d.hdf5'


flags = 'rayleighBenard-'
if run_type != MushyLayerRun.UNIFORM_GRID and max_level > 0:
    flags = flags + run_type + '-'


for RaT in Ra_arr:

    numCellsVertical = int(grid_res)
    numCellsHorizontal = int(round(grid_res*aspect_ratio))

    #params['main.max_box_size'] = str(MushyLayerRun.computeMinBoxSize(grid_res, numCellsHorizontal, numCellsVertical))
    params['main.max_box_size'] = '124'

    grids_file = ''
    
    # Only support unifomr grids for now
    #if run_type == MushyLayerRun.FIXED_GRID and max_level > 0:
    #    params['main.gridsFile'] = ('grids/grids-minBox' +
    #                                params['main.max_box_size'] +
    #                                '-res-' + str(res) + "-" + str(res*2))
    #    params['main.regrid_interval'] = '-1'  # make sure we don't adapt the grid

    if run_type == MushyLayerRun.ADAPTIVE_GRID:
        params['main.tagging_val'] = str(0.3 * float(8/grid_res)) #0.08
        params['main.taggingVar'] = '0'
        #params['main.tags_grow'] = str(int(numpy.ceil(4*(float(res)/float(128)))))
        params['main.regrid_interval'] = '16'

    params['parameters.rayleighTemp'] = str(RaT*Pr)
    params['parameters.darcy'] = str(Pr)
    params['parameters.reynolds'] = str(1/Pr)

    output_dir = flags+'Ra'+str(RaT)
    plot_prefix = output_dir+'-'
    check_prefix = plot_prefix + 'chk-'

    # parameters = ' '.join([key+'='+value for key, value in params.items()])

    ml_run = MushyLayerRun(output_dir,  plot_prefix,  check_prefix,
                           max_level,  num_proc,  plot_interval,
                           check_interval, gridType, params, subcycled)
    
    ml_run.custom_grid = [numCellsHorizontal, numCellsVertical, numCellsHorizontal]


    ml_run.run([grid_res])
    
   # ml_run.continue_run('128-0', 'bm2-Ra200-chk-128pts1800.2d.hdf5')

    # Copy data to dropbox so we can analyse it
    if ml_run.machine_specific.has_dropbox():
        ml_run.copy_final_output(grid_res, dropbox_data_folder)

    # Copy data from gyre servers to cloud
    if ml_run.machine_specific.is_gyre():
        ml_run.copy_all_output_to_cloud(grid_res)

    # ml_run.continueRun('64-2')
