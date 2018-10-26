from MushyLayerRun import MushyLayerRun
import os

os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
output_base_dir = 'subcycled-solidificationNoFlow'

subcycled = True

params = {}

# These are the options to specify
run_type = MushyLayerRun.ADAPTIVE_GRID
max_level = 1
num_proc = 1
plot_interval = 100
check_interval = 100
params['main.max_box_size'] = '1024'
params['main.verbosity'] = '1'
gridType = MushyLayerRun.SQUARE
grid_res = [64]

# Parameters for this benchmark problem
params['parameters.darcy'] = 0.0
params['parameters.rayleighTemp'] = 0.0
params['parameters.stefan'] = 5
params['parameters.compositionRatio'] = 2
params['parameters.K'] = 1
params['parameters.specificHeatRatio'] = 1
#params['parameters.lewis'] = 1e300
params['parameters.V'] = 1e-6
params['parameters.topEnthalpy'] = 2.5 #need to ensure this is the top of the eutectic
params['parameters.topEnthalpy'] = 7   # make sure this is something above the liquidus
params['parameters.topBulkConc'] = 1.0
params['parameters.bottomBulkConc'] = 1.0

params['main.max_time'] = '99999.0'
params['main.max_step'] = '99999'
params['main.steady_state'] = '1e-3'
params['main.enforceAnalyticSoln'] = 'true' # this will speed things up a bit

params['params_file'] = '../MushyLayer/execSubcycle/solidificationNoFlow.parameters'

flags = ''
if run_type != MushyLayerRun.UNIFORM_GRID and max_level > 0:
    flags = run_type

for res in grid_res:
    grids_file = ''
    if run_type == MushyLayerRun.FIXED_GRID and max_level > 0:
#         params['gridfile'] = ('grids/grids-minBox' +
#                                     params['main.max_box_size'] +
#                                     '-res-' + str(res) + "-" + str(res*2))
        params['main.gridfile'] = ('grids/grids-viscous-'  + str(res) + "-"
                              + str(res*2))
        
    if run_type == MushyLayerRun.ADAPTIVE_GRID:
        params['main.tagging_val'] = str(0.05*16/res)

    output_dir = output_base_dir
    if flags:
        output_dir = output_base_dir+'-' + flags
    plot_prefix = output_dir + '-'
    check_prefix = plot_prefix
    # parameters = ' '.join([key+'='+value for key, value in params.items()])

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
