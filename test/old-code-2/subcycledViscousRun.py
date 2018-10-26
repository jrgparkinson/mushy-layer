from MushyLayerRun import MushyLayerRun
import os

os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
output_base_dir = 'poiseuille'

subcycled = True

params = {}

# These are the options to specify
run_type = MushyLayerRun.ADAPTIVE_GRID
max_level = 1
num_proc = 4
plot_interval = 50
check_interval = 50
params['main.max_box_size'] = '8'
params['main.verbosity'] = '3'
gridType = MushyLayerRun.SQUARE
grid_res = [8, 16, 32]

# These parameter should be the same for every run
# params['parameters.darcy'] = 1.0
# params['parameters.rayleighTemp'] = 1.0
# params['parameters.stefan'] = 0  # 0 means H = theta
# params['parameters.K'] = 1
# params['parameters.lewis'] = 10e200
# params['parameters.V'] = 0

params['main.max_time'] = '100.0'
params['main.steady_state'] = '1e-9'
params['parameters.porosityFunction'] = 2
params['parameters.fixed_porosity'] = 0.5

params['params_file'] = '../MushyLayer/execSubcycle/poiseuille.parameters'
params['main.enforceAnalyticSoln'] = 'true'

flags = ''

if run_type != MushyLayerRun.UNIFORM_GRID and max_level > 0:
    flags = run_type



#params['periodic_bc'] = '0 1 1'  # periodic in vertical direction

for res in grid_res:
    grids_file = ''
    if run_type == MushyLayerRun.FIXED_GRID and max_level > 0:
#         params['gridfile'] = ('grids/grids-minBox' +
#                                     params['main.max_box_size'] +
#                                     '-res-' + str(res) + "-" + str(res*2))
        params['main.gridfile'] = ('grids/grids-viscous-'  + str(res) + "-"
                              + str(res*2))
        
#    if run_type == MushyLayerRun.ADAPTIVE_GRID:
#        params['main.refine_thresh'] = str(0.05*32/res)

    output_dir = output_base_dir
    if flags:
        output_dir = output_base_dir+'-' + flags
    plot_prefix = output_dir + '-'
    check_prefix = plot_prefix
    #parameters = ' '.join([key+'='+value for key, value in params.items()])
    
    #print(parameters)

    ml_run = MushyLayerRun(output_dir,  plot_prefix,  check_prefix,
                           max_level,  num_proc,  plot_interval,
                           check_interval, gridType, params, subcycled)

    ml_run.run([res])

    # Copy data to dropbox so we can analyse it
    if ml_run.machine_specific.has_dropbox():
        ml_run.copy_final_output(res, dropbox_data_folder)

    # Copy data from gyre servers to cloud
    if ml_run.machine_specific.is_gyre():
        ml_run.copy_all_output_to_cloud(res)

    # ml_run.continueRun('64-2')
