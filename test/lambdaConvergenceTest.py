# Script to test convergence of uniform and AMR solutions

import os
import sys
parDir = os.path.abspath(os.pardir)
pythonDir = os.path.join(parDir, 'python')
sys.path.append(pythonDir)

from MushyLayerRun import MushyLayerRun
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile
import time

import getopt
from SlurmTask import SlurmTask

def str2arr(string):
    #print(string)
    parts = string.split(' ')
    array = [int(i) for i in parts]
    return array

def arr2str(array):
    str_array = [str(a) for a in array]
    string = ' '.join(str_array)
    return string

def isPowerOfTwo(n):
    test = math.log(n)/math.log(2)
    
    if round(test) == test:
        return True
    else:
        return False

def main(argv):
    
    # Defaults:
    num_proc = 1
    
    try:
       opts, args = getopt.getopt(argv,"n:f:H:")
    except getopt.GetoptError:
       print 'lambdaConvergenceTest.py -n num_proc'
       sys.exit(2)
    for opt, arg in opts:
        if opt in ("-n"):
            num_proc = int(arg)
        
    #print 'Num processors: ', str(num_proc)

    #################################
    # These shouldn't change
    #################################
    cwd = os.getcwd()
    mushyLayerBaseDir = cwd.replace('/run', '')

    os.environ["CH_TIMER"] = "1"
    dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
    subcycled = True # subcycled code

    # default:
    paramsToScale = ['main.block_factor', 'main.max_grid_size', 'main.tag_buffer_size']

    # These are the main options
    # Assume params file contains coarsest params desired, and compute refinements
    
    #execFile = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex'
    #base_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/'

    base_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/lambdaConvergenceTest/'

    #params_file = 'inputsUniformPorosity'
    original_params = readInputs(base_dir + 'inputs')
    #original_params['regrid.fixed_grid_time']=0.001

    # We should always be doing this now
    original_params['main.advectionMethod'] = 1

    original_params['projector.verbosity'] = 3
    

    run_name = 'coupled'
 
    # Options:
    original_params['projection.eta'] = 0.95
    original_params['main.turn_off_buoyancy_time'] = -1 #0.29
    original_params['regrid.fixed_grid_time'] = 0.3 
    original_params['projection.scaleSyncCorrection'] = 0
    original_params['main.max_eta'] = original_params['projection.eta']
    original_params['regrid.min_regrid_time'] = 0.3


    allPorous = False
    if allPorous:
        original_params['bc.enthalpyHiVal'] = '0.0 2.0'
        original_params['bc.enthalpyLoVal'] = '0.0 4.2'
        original_params['parameters.stefan'] = 5.0
        original_params['parameters.rayleighComp'] = 6e6

    tagML = True
    if tagML:
        
        original_params['regrid.fixed_grid_time'] = -1
        original_params['regrid.tagMushLiquidBoundary'] = 1
        original_params['main.block_factor'] = 4
        original_params['main.tag_buffer_size'] = 1
        original_params['main.plot_interval'] =10
        original_params['main.vel_refine_thresh'] = 10000 # just to ensure we don't do any other tagging
        original_params['main.grid_buffer_size'] ='2 2'
        original_params['main.regrid_interval'] = '64 128'
        paramsToScale = ['main.max_grid_size', 'main.tag_buffer_size', 'main.regrid_interval']


    outflow = False
    if outflow:
        original_params['bc.scalarLo']='1 2 0'
        original_params['bc.velLo']='6 3 0'

    # Construct final run name:
    run_name = run_name + 'Eta' + str(original_params['projection.eta'])

    if 'main.turn_off_buoyancy_time' in original_params and original_params['main.turn_off_buoyancy_time'] >= 0.0:
        run_name = run_name + 'NoBuoyancy'

    if 'projection.scaleSyncCorrection' in original_params and original_params['projection.scaleSyncCorrection'] == 1:
        run_name = run_name + 'ScaleSync'

    if 'regrid.min_regrid_time' in original_params and original_params['regrid.min_regrid_time'] >= 0:
        run_name = run_name + 'MinRegridTime' + str(original_params['regrid.min_regrid_time'])

    if allPorous:
        run_name = run_name + 'AllPorous'

    if tagML:
        run_name = run_name + 'TagMLV3'

    if outflow:
        run_name = run_name + 'Outflow'



    coupled = True

    # Don't want to have to uncomment all these lines, so just turn them off
    if not coupled:

       # run_name = 'uniformMiddleGridDarcyPorosity0.5NoEta'
        original_params['bc.porosityLoVal'] = '0.5 0.5'
        original_params['parameters.problem_type']=14
        original_params['parameters.porosityFunction']=0
        #original_params['regrid.fixed_grid_time']=0.03

        #run_name = 'coupledMiddleGridNoEtaNoBuoyancy'
        #run_name = 'uniformMiddleGridDarcyPorosity0.5NoEtaNoBuoyancy'
        run_name = 'Porosity1NoEtaNoBuoyancy'
        original_params['bc.porosityLoVal'] = '1.0 1.0'
        original_params['main.turn_off_buoyancy_time'] = 0.018
        original_params['regrid.fixed_grid_time']=0.02

        run_name = 'Porosity0.5NoEtaNoBuoyancy'
        original_params['bc.porosityLoVal'] = '0.5 0.5'

        run_name = 'Porosity0.5NoEtaNoBuoyancyScalePhi'
        #original_params['projection.phiScale'] = 0.5

        run_name = 'Porosity0.5NoEtaNoBuoyancyNoPorosityAdv'
        original_params['main.advectionMethod'] = 2

        run_name = 'Porosity0.5NoEtaNoBuoyancyNoPorosityAdvNoReflux'
        #original_params['main.reflux_momentum'] = 0

        run_name = 'Porosity0.5NoEtaNoBuoyancyNoPorosityAdvNoViscSrc'
        #original_params['advSrc.viscous'] = 0

        run_name = 'Porosity0.5NoEtaNoPorosityAdvection'
        original_params['main.darcyCoeff'] = 1
        original_params['projection.scalePressureWithPorosity'] = 1
        original_params['bc.porosityLoVal'] = '0.5 0.5'
        original_params['parameters.rayleighTemp'] = 2e5

        run_name = 'Porosity0.5NoEta'
        original_params['main.advectionMethod'] = 0

        run_name = 'Porosity0.5NoEtaCorrectedAdvection'
        original_params['main.advectionMethod'] = 0

        run_name = 'Porosity0.5NoEta'
        original_params['main.turn_off_buoyancy_time'] = -0.018
        original_params['main.advectionMethod'] = 1
        original_params['parameters.rayleighTemp'] = 5e4
        original_params['regrid.fixed_grid_time']=0.5
        original_params['main.max_time'] = 1.0

        run_name = 'Porosity1.0NoEta'
        original_params['bc.porosityLoVal'] = '1.0 1.0'
        original_params['parameters.rayleighTemp'] = 2.5e4


        run_name = 'VariablePorosity0.5NoEta'
        original_params['parameters.rayleighTemp'] = 1e3
        original_params['parameters.porosityFunction']=3
        original_params['main.maxChi']=0.5
        original_params['main.stdev'] = 0.03
        original_params['regrid.fixed_grid_time']=0.2
        original_params['main.max_time'] = 1.0

        #run_name = 'Porosity1.0NoDependence'
        #original_params['bc.porosityLoVal'] = '1.0 1.0'
        #original_params['parameters.rayleighTemp'] = 1e5

        #run_name = 'Porosity0.5NoEtaNoPScale'
        #original_params['projection.scalePressureWithPorosity'] = 0

        #run_name = 'Porosity0.5NoEtaNoPorosityAdvection'
        #original_params['main.advectionMethod'] = 2


        #run_name = 'Porosity0.5NoEtaScaledAdvectionSrc'
        #original_params['main.scaleAdvectionSrcExtraChi'] = 1
        #original_params['main.advectionMethod'] = 0

        #run_name = 'Porosity0.5NoEtaScaledProjBC'
        #run_name = 'Porosity0.5NoEtaNew'
        #original_params['projection.phiScale'] = 0.5
        #run_name = 'uniformMiddleGridPorosity0.5NoEtaNoBuoyancyNoDarcy'
        #original_params['main.darcyCoeff'] = 0.0

        #run_name = 'Porosity0.5NoEtaNoBuoyancyNoDarcyNoMomReflux'
        #original_params['main.reflux_momentum'] = 0

        #run_name = 'Porosity0.5NoEtaNoViscosity'
        #original_params['main.viscosityCoeff'] = 0.0

    
    



    #gridType = 'top'
    #run_name = 'coupledTopGrid'
    

    refines = [1, 2, 4]

    

    original_params['main.plot_interval'] = 50
    original_params['main.checkpoint_interval'] = 200
    

    for refinement in refines:

        params = dict(original_params)

        #refinement = int(pow(2, i))

        # Params that need to be scaled: 
        # block factor, grid size, tag_buffer_size, grid_buffer_size
        
        for par in paramsToScale:
            params[par]  = refinement*int(original_params[par])
       
        vectorParamsToScale = ['main.num_cells', 'main.grid_buffer_size']
        for par in vectorParamsToScale:  
            #print(par + '=' + original_params[par])          
            old_par = str2arr(original_params[par])
            new_par = [int(n * refinement) for n in old_par]

            #print(par + '=' + arr2str(new_par))
            #params[par] = arr2str(new_par)
            params[par] = new_par
      
        resolution = params['main.num_cells'] #str2arr(params['main.num_cells'])
        #params['main.regrid_gridfile'] = base_dir + '/' + str(resolution[0]) + 'x' + str(resolution[1]) + 'channel'
        params['main.regrid_gridfile'] = '/home/parkinsonj/convection-in-sea-ice/MushyLayer/grids/middle/' + str(resolution[0]) + 'x' + str(resolution[1])

        params['main.output_folder'] = run_name #'uniformChiTopGridLonger'

        res_run_name = params['main.output_folder'] + '-Nx' + str(params['main.num_cells'][0])
        output_folder = base_dir + res_run_name

        params['main.output_folder'] = output_folder

        if not os.path.exists(output_folder):
            #os.makedirs(output_folder)
            #time.sleep(0.1) # wait long enough to ensure folder is made
            params_file = os.path.join(output_folder, 'inputs')
            print('Params file: ' + params_file)           

            # Do the run
            #writeParamsFile(params_file, params)
            #cmd = 'cd ' + output_folder + '; ' + execFile + ' inputs'
            #os.system(cmd)
            # Do run
            s = SlurmTask('', res_run_name, '', num_proc)

            ml_run = MushyLayerRun(base_dir, num_proc, params, subcycled)
            ml_run.allowMultipleOutputDirs = False
            ml_run.makeSlurmJob = True
            ml_run.slurmJob = s
            ml_run.single_run(res_run_name)

        

if __name__ == "__main__":
    print('Running script')
    main(sys.argv[1:])
                 


