# Script to run mushy layer simulations for a variety of parameters
# at fixed resolution


import os
import sys
parDir = os.path.abspath(os.pardir)
pythonDir = os.path.join(parDir, 'python')
sys.path.append(pythonDir)

from MushyLayerRun import MushyLayerRun
#import numpy
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs

#################################
# These shouldn't change
#################################
#cwd = '/home/parkinsonj/convection-in-sea-ice/MushyLayer'
cwd = os.getcwd()
mushyLayerBaseDir = cwd.replace('/python', '')
#params = readInputs(mushyLayerBaseDir + '/params/wettlaufer1997.parameters')
#params = readInputs(mushyLayerBaseDir + '/params/middleton2016.parameters')
params = readInputs(mushyLayerBaseDir + '/params/mushyLayerDarcy.parameters')
os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
subcycled = True
base_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/3D/'
execFile = 'mushyLayer3d.Linux.64.mpiCC.gfortran.OPTHIGH.MPI.ex'

###################################
# Parameters to change
###################################

params['main.multiCompSolve'] = 'true'
params['picardSolver.max_iter'] = '1'

params['main.initial_cfl'] = '0.001'
params['main.ignoreSolverFails'] = 'true'

# --------------------------------------------
# Grid setup
params['main.num_cells'] = '64 64 64'
params['main.domain_length'] = '0.5'
params['main.periodic_bc'] = '1 1 0'

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

params['main.plot_prefix'] = output_dir + '-' + run_name + '-'
params['main.chk_prefix'] =  'chk' + params['main.plot_prefix']

params['main.initialPerturbation'] = '0.0'

params['main.checkpoint_interval'] = '20'
params['main.plot_interval'] = '20'
params['main.max_step'] = '999999'

params['parameters.stefan'] 			= '5' 
params['parameters.compositionRatio'] = '3.0'
params['parameters.rayleighComp']				= '200.0' 
params['parameters.darcy']                    = '0.0'

params['bc.scalarLo'] = '1 1 2'
params['bc.scalarHi'] = '1 1 0'

params['bc.velLo'] = '0 0 4' # x-dir, y-dir, z-dir
params['bc.velHi'] = '0 0 0' # x-dir, y-dir, z-dir

# -----------------------------------------
# Parallel stuff
num_proc = 8

# Construct name for run
run_name = constructRunName(params) 

# Do run
ml_run = MushyLayerRun(base_dir, num_proc, params, subcycled, execFile)
ml_run.single_run(run_name)

