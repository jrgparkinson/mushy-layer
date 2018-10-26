from MushyLayerRun import MushyLayerRun
import os
import math

SIN = 1
PERIODIC = 2
burgersProblem = SIN

os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"

if burgersProblem == SIN:
    output_base_dir = 'burgers-sin-noreflux'
else:
    output_base_dir = 'burgers-periodic'


subcycled = True

params = {}

# These are the options to specify
run_type = MushyLayerRun.FIXED_GRID
max_level = 1
num_proc = 1
params['main.max_box_size'] = '256'
params['main.verbosity'] = '3'
gridType = MushyLayerRun.SQUARE
resolutions = [16, 32, 64, 128]

# These parameter should be the same for every run
if burgersProblem == SIN:
    params['main.periodic_bc'] = '0 1 1'
    params['main.domain_length'] = '2.0'
    params['main.symmetric_domain'] = 'true'
   # params['parameters.darcy'] = 0.01/math.pi
    params['main.params_file'] = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/burgersSin.parameters'
else:
    params['main.periodic_bc'] = '1 1 1'
    params['main.domain_length'] = str(2*math.pi)
    params['main.symmetric_domain'] = 'false'
    params['main.params_file'] = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/burgersPeriodic.parameters'


# params['main.max_step'] = '1'
params['main.max_time'] = str(1e-1)
params['main.use_limiting'] = 'true'
params['main.fixed_dt'] = str(1e-3)

flags = ''
if run_type != MushyLayerRun.SINGLE_LEVEL and max_level > 0:
    flags = run_type

if run_type == MushyLayerRun.ADAPTIVE_GRID:
    params['main.tagging_val'] = '100'


output_dir = output_base_dir
if flags:
    output_dir = output_base_dir + '-' + flags

plot_prefix = output_dir + '-'
check_prefix = plot_prefix + 'chk-'

for res in resolutions:

    grids_file = ''
    if run_type == MushyLayerRun.FIXED_GRID and max_level > 0:
        params['main.gridfile'] = ('grids/grids-burgers-'  + str(res) + "-"
                              + str(res*2))


    # parameters = ' '.join([key+'='+value for key, value in params.items()])
    
    plot_interval = 1
    check_interval = 50

    ml_run = MushyLayerRun(output_dir,  plot_prefix,  check_prefix,
                           max_level,  num_proc,  plot_interval,
                           check_interval, gridType, params, subcycled)

    ml_run.run([res])

    # Copy data to dropbox so we can analyse it
    if ml_run.machine_specific.has_dropbox():
        ml_run.copy_final_output(res, dropbox_data_folder)

    # Copy data from gyre servers to cloud
    #if ml_run.machine_specific.is_gyre():
    #ml_run.copy_all_output_to_cloud(res)

    # ml_run.continueRun('64-2')
