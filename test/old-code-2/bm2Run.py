from MushyLayerRun import MushyLayerRun
import os
import numpy

os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"

# different options within the general idea of convection in a fixed porous medium
DEFAULT = 'bm2-'
SIDEWALL = 'bm2-sidewall-'
SIDEWALL_PERIODIC = 'bm2-sidewall-periodic-'

params = {}

# These are the options to specify
run_type = MushyLayerRun.UNIFORM_GRID
bm2_subtype = DEFAULT
max_levels = [0]
num_proc = 8
plot_interval = 100
check_interval = 100
params['main.verbosity'] = '3'
gridType = MushyLayerRun.SQUARE
grid_res = [128]
params['main.max_box_size'] = '32'
RaT_arr = [250, 300]  # 50,100,200,250,300

# These parameter should be the same for every run
# params['parameters.rayleighTemp'] = 200
params['parameters.rayleighComposition'] = 0
params['parameters.stefan'] = 0  # 0 means H = theta
params['parameters.K'] = 1
params['parameters.lewis'] = 10e200
params['parameters.V'] = 0
params['parameters.topTemp'] = 0
params['parameters.bottomTemp'] = 10
params['parameters.eutecticTemp'] = 0
params['parameters.deltaTemp'] = 10
params['parameters.initialComposition'] = 0.7
params['parameters.problem_type'] = 2
params['main.enforceAnalyticSoln'] = 'false'
#params['parameters.porosityTimescale'] = 0

#params['main.restart_file'] = ''


for max_level in max_levels:
    
    flags = ''
    if run_type != MushyLayerRun.UNIFORM_GRID and max_level > 0:
        flags = run_type + '-'
    
    
    
    if bm2_subtype == DEFAULT:
        params['parameters.fixedTempDirection'] = '1'
    elif bm2_subtype == SIDEWALL:
        params['parameters.fixedTempDirection'] = '0'
    elif bm2_subtype == SIDEWALL_PERIODIC:
        params['parameters.fixedTempDirection'] = '0'
        params['main.is_periodic'] = '0 1 0'  # periodic in vertical direction
    
    
    for RaT in RaT_arr:
        for res in grid_res:
            
            grids_file = ''
            if run_type == MushyLayerRun.FIXED_GRID and max_level > 0:
                params['main.gridsFile'] = ('grids/grids-minBox' +
                                            params['main.max_box_size'] +
                                            '-res-' + str(res) + "-" + str(res*2))
                params['main.regrid_interval'] = '-1'  # make sure we don't adapt the grid
                
            if run_type == MushyLayerRun.ADAPTIVE_GRID:
                params['main.tagging_val'] = str(0.08 * float(32/res)) #0.08
                params['main.tagging_scalar_var'] = '8'
                params['main.tags_grow'] = str(int(numpy.ceil(4*(float(res)/float(128)))))
                params['main.regrid_interval'] = '1'
                 
    
            params['parameters.rayleighTemp'] = str(RaT)
    
            output_dir = bm2_subtype+flags+'Ra'+str(RaT)
            plot_prefix = bm2_subtype+flags+'Ra'+str(RaT)+'-'
            check_prefix = plot_prefix + 'chk-'
    
            # parameters = ' '.join([key+'='+value for key, value in params.items()])
    
            ml_run = MushyLayerRun(output_dir,  plot_prefix,  check_prefix,
                                   max_level,  num_proc,  plot_interval,
                                   check_interval, gridType, params)
            
            
    
            ml_run.run([res])
            
           # ml_run.continue_run('128-0', 'bm2-Ra200-chk-128pts1800.2d.hdf5')
    
            # Copy data to dropbox so we can analyse it
            if ml_run.machine_specific.has_dropbox():
                ml_run.copy_final_output(res, dropbox_data_folder)
    
            # Copy data from gyre servers to cloud
            if ml_run.machine_specific.is_gyre():
                ml_run.copy_all_output_to_cloud(res)
    
            # ml_run.continueRun('64-2')
