# Run all the parameter files listed below for 20 timesteps on an AMR Hierarchy
# and just check that nothing breaks
from MushyLayerRun import MushyLayerRun
import os

param_files = ['poiseuille', 'hrl', 'solidificationNoFlow', 'mushyLayer']
#param_files = ['mushyLayer']
max_level = 2
num_proc = 4
plot_interval = 1
check_interval = 20
gridType = MushyLayerRun.SQUARE
subcycled = True

execDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle'
output_base_dir = 'test-output'

if not (os.path.exists(output_base_dir)):
    os.mkdir(output_base_dir)

params = {}

params['main.max_step'] = '20'
params['main.enforceAnalyticSoln'] = 'false'
params['main.gridfile'] = ('grids/grids-minBox8-res-16-32-64')
params['main.max_box_size'] = '8'


for filename in param_files:

    output_dir = output_base_dir + '/' + filename
    plot_prefix = filename
    check_prefix = filename

    params['params_file'] = '../MushyLayer/execSubcycle/'+filename+'.parameters'

    ml_run = MushyLayerRun(output_dir,  plot_prefix,  check_prefix,
                           max_level,  num_proc,  plot_interval,
                           check_interval, gridType, params, subcycled)

    exitStatus = ml_run.run([16])
