# Script to run mushy layer simulations for a variety of parameters
# at fixed resolution
import os
import sys
import math
from subprocess import Popen
import getopt

pythonDir = os.path.abspath(os.pardir)
sys.path.append(pythonDir)
from util.mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile, str2arr
from classesMushyLayerRun import MushyLayerRun
from classe.sSlurmTask import SlurmTask

def main(argv):
    num_proc = 12
    runOnLegacy = False
 
    try:
       opts, args = getopt.getopt(argv,"ln:")
    except getopt.GetoptError:
       print 'fixedChill.py  -l (run on legacy) -n <num processors>'
       sys.exit(2)
    for opt, arg in opts:
        if opt in ("-n"):
            num_proc = int(arg)
        elif opt in ("-l"):
            runOnLegacy = True
        
          
    print 'Num processors: ', str(num_proc)
    
    #################################
    # These shouldn't change
    #################################
    cwd = '/home/parkinsonj/convection-in-sea-ice/MushyLayer'
    mushyLayerBaseDir = cwd.replace('/run', '')
    
    os.environ["CH_TIMER"] = "1"
    dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
    subcycled = True
    base_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/fixedChill/'

    params = readInputs(os.path.join(base_dir ,'inputs'))
    

    # Typical params to vary:
    #params['parameters.lewis'] = '200'
    params['parameters.darcy'] = 1e-2 #1e-2
    params['parameters.rayleighComp'] = 50
    params['main.num_cells']='64 128'
    params['main.num_cells']='128 256'
    params['main.num_cells']='256 512'

    params['main.num_cells']='256 768'
    #params['main.num_cells']='512 1536'
    params['main.num_cells']='512 1024'

    #params['main.num_cells']='64 192'

    if num_proc < 1:
        ncells = str2arr(params['main.num_cells'])
        totalCells = ncells[0]*ncells[1]

        max_num_proc = 32
        params['main.max_grid_size'] = 32

        # Ideal number of processors - divide domain
        num_proc = math.ceil(totalCells/(params['main.max_grid_size']*params['main.max_grid_size']))

        while num_proc > max_num_proc:
            num_proc = math.ceil(num_proc/2)

        print('Using ' + str(num_proc) + ' processors')

    #params['main.domain_width']=1.7
    params['parameters.compositionRatio'] = 1.18
    params['parameters.nonDimReluctance'] = 1e-1

    # Darcy flow:
    darcyFlow = True
    if darcyFlow:

        params['parameters.nonDimReluctance'] = 0.05
        params['main.domain_width']=0.02
        params['parameters.rayleighComp'] = 1000

        params['main.domain_width']=1.0
        params['parameters.rayleighComp'] = 85
        params['parameters.nonDimReluctance'] = 0.1

        params['main.initial_cfl']=2e-5
        params['main.max_dt_growth']=1.05

        params['main.time_integration_order']=2 # need this (or very small init cfl)

        params['main.doEuler'] = 0
        params['parameters.darcy'] = 0.0 # doesn't matter
        params['parameters.prandtl'] = 0 # da = pr
        #params['parameters.rayleighComp'] = 1000
        #params['parameters.rayleighComp'] = params['parameters.rayleighComp']/params['parameters.darcy']
        params['main.plot_interval'] = 50
        params['main.nondimensionalisation'] = 0

        params['bc.enthalpyHiVal']='0 0.9'
        params['main.periodic_bc']='1 0 0'



    #params['parameters.darcy'] = 1e-8
    #params['parameters.rayleighComp'] = 1e6
    #params['main.num_cells']='128 256'
    #params['main.domain_width']=0.1
    #params['parameters.compositionRatio'] =  1.2 #20.0

    restDir = 'CR1.200RaC1.0e+06Le100KozenyPermeabilityDa1.0e-08R1.0e-03pts128-St20-domWidth0.1-stoppedEarly'
    restChkFile = 'chkfixedChill-periodic-CR1.200RaC1.0e+06Le100KozenyPermeabilityDa1.0e-08R1.0e-03pts128-St20-domWidth0.1-001200.2d.hdf5'
    #params['main.restart_file'] = os.path.join(base_dir, restDir, restChkFile)
    #params['main.verbosity'] = 10
    
    params['parameters.stefan'] = 5.0

    # Bottom temp: 1.03
    bottom_temperature = 1.3
    params['bc.enthalpyLoVal']='0.0 ' + str(float(params['parameters.stefan']) + bottom_temperature)


    # New params for AMR:
    AMR = True
    if AMR:
        params['main.max_level']=1
        params['main.refine_thresh']=0.15
        params['main.regrid_interval']='64 64 64'
        params['regrid.plume_vel'] = 8
        params['regrid.plume_salinity'] = -0.95
        params['main.tag_buffer_size']=4
        params['main.fill_ratio']=0.85
        params['main.grid_buffer_size']=[0, 0, 0]
        params['projection.eta']=0.95

        # General stuff:
        params['main.domain_width']=1.0
        params['main.num_cells']='64 128'
        params['main.num_cells']='128 256'
        params['main.num_cells']='256 512'

		# Physics
        params['parameters.rayleighComp']     = 350 #100
        params['parameters.nonDimReluctance'] = 0.5 #1e-5?
        params['parameters.compositionRatio'] = 1.4

		# Perturbation:
        params['main.perturbationWavenumber'] = 2.0
        params['main.perturbationTime']       = 0.1
        params['main.delayedPerturbaation']   = 0.1
        params['main.perturbationPhaseShift'] = 0.25



    # Seems to work:
    if float(params['parameters.darcy']) > 1e-10:
        params['main.max_dt']=1/float(params['parameters.darcy'])
    else:
        params['main.max_dt']=params['parameters.nonDimReluctance']

    params['main.max_time'] = 10000*params['main.max_dt']


    #params['main.restart_file'] = base_dir + 'CR1.15RaC400Le200ChiCubedPermeabilitypts64-domWidth5.0-stoppedEarly/chkfixedChill-periodic-CR1.15RaC400Le200ChiCubedPermeabilitypts64-domWidth5.0-007000.2d.hdf5'
    
    if 'main.max_level' not in params or params['main.max_level'] == '0':
        run_type = MushyLayerRun.UNIFORM_GRID
        params['main.max_level'] == '0'
    else:
        run_type = MushyLayerRun.ADAPTIVE_GRID
        params['main.max_level'] == '1'
     
    # Construct name for run
    output_dir = 'fixedChill'
    if params['main.periodic_bc'] == '0 0 1':
        output_dir = output_dir + '-insulating'
            
        # To force channels to form at edge of domain
        #params['main.initialPerturbation'] = '0.1'
    else:
        output_dir = output_dir + '-periodic'

        # This appends "fixedGrid" or "adaptive" if appropriate
    if run_type != MushyLayerRun.UNIFORM_GRID and int(params['main.max_level']) > 0:
        output_dir = output_dir + '-' + run_type


    run_name = constructRunName(params) 
    run_name = run_name + '-TBottom' + str(bottom_temperature) + '-domWidth' + str(params['main.domain_width'])

    if not params['main.max_level'] == '0':
        run_name = run_name + 'AMR' + str(params['main.max_level'])

    params['main.plot_prefix'] = 'plt' #output_dir + '-' + run_name + '-'
    params['main.chk_prefix'] =  'chk' # + params['main.plot_prefix']
    
       
    # Do run
    s = SlurmTask('', run_name, '', num_proc)

    if runOnLegacy:
        s.partitions='legacy'

    ml_run = MushyLayerRun(base_dir, num_proc, params, subcycled)
    ml_run.allowMultipleOutputDirs = False
    ml_run.makeSlurmJob = True
    ml_run.slurmJob = s
    ml_run.single_run(run_name)


if __name__ == "__main__":
    print('Running script')
    main(sys.argv[1:])
     
