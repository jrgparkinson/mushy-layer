# Script to run mushy layer simulations for a variety of 
# parameters to try and understand what sets the channel spacing

import os
import sys
parDir = os.path.abspath(os.pardir)
pythonDir = os.path.join(parDir, 'python')
sys.path.append(pythonDir)

from MushyLayerRun import MushyLayerRun
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile

#################################
# These shouldn't change
#################################
cwd = '/home/parkinsonj/convection-in-sea-ice/MushyLayer'
mushyLayerBaseDir = cwd.replace('/run', '')

#params = readInputs(mushyLayerBaseDir + '/params/mushyLayerDarcyLowC.parameters')

base_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/channelSpacing/'
#prevDir = os.path.join(base_dir, 'CR1.25RaC0Le200ChiCubedPermeabilitypts1024-1/')
#restartFile = os.path.join(prevDir,'chkchannelSpacing-CR1.25RaC0Le200ChiCubedPermeabilitypts1024-056000.2d.hdf5')

#prevDir = os.path.join(base_dir, 'CR1.25RaC0Le200ChiCubedPermeabilitypts2048-unsteady/')
#restartFile = os.path.join(prevDir, 'chkchannelSpacing-CR1.25RaC0Le200ChiCubedPermeabilitypts2048.2d.hdf5')

#prevDir = os.path.join(base_dir, 'CR1.25RaC0Le200ChiCubedPermeabilitypts2048-insulating/')
#restartFile = os.path.join(prevDir, 'CR1.25RaC0Le200ChiCubedPermeabilitypts2048-insulating.2d.hdf5')

prevDir = os.path.join(base_dir, 'CR20RaC50Le200ChiCubedPermeability2048-0/')
restartFile = os.path.join(prevDir, 'chkCR20RaC50Le200ChiCubedPermeability.2d.hdf5')

params = readInputs(os.path.join(prevDir, 'inputs'))
params['main.restart_file'] = restartFile

os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
subcycled = True

#output_dir = base_dir
#output_dir = '/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/channelSpacing'

# -----------------------------------------
# Parallel stuff
num_proc = 4

params['main.checkpoint_interval'] = '1000'
params['main.plot_interval'] = '100'

params['main.wavenumber'] = '0'
params['main.maxRestartWavenumber'] = '50'


#params['main.verbosity'] = '10'

# First just get steady state
Ra_arr = [300, 200, 150]
Ra_arr = [40, 50, 60]

for Ra in Ra_arr:
    # for testing
    params['parameters.rayleighComp'] = str(Ra)
    
    
    # Try Larger perturbation for smaller Ra
    pertSize = 0.1*250/Ra
    pertSize = 0.0
    params['main.restart_perturbation'] = str(pertSize)
    
    
    run_name = constructRunName(params)
    run_name = run_name + '-noPerturbation'
    #run_name = run_name + '-Nx' + str(Nx)
    params['main.plot_prefix'] = run_name + '-'
    params['main.chk_prefix'] =  'chk' + params['main.plot_prefix']
    
    # Do run
    ml_run = MushyLayerRun(base_dir, num_proc, params, subcycled)
    ml_run.allowMultipleOutputDirs = False # If the run ending -0 has already completed, don't do a run ending -1
    ml_run.single_run(run_name)



	
		

