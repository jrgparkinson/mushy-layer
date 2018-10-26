# Script to test the lambda correction for a hierarchy
# of problems
# May later extend to consider convergence

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
       print 'simpleConvergenceTest.py'
       sys.exit(2)
    for opt, arg in opts:
        if opt in ("-n"):
            num_proc = int(arg)
        
    print 'Num processors: ', str(num_proc)

    #################################
    # These shouldn't change
    #################################
    cwd = os.getcwd()
    mushyLayerBaseDir = cwd.replace('/run', '')

    os.environ["CH_TIMER"] = "1"
    dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
    subcycled = True # subcycled code

    # These are the main options
    # Assume params file contains coarsest params desired, and compute refinements
    execFile = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex'
    base_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/'
    params_file = 'inputsUniformPorosity'
    num_refine = 1

    # Set debug true for lots of extra output
    debug = False

    original_params = readInputs(base_dir + params_file)
    
    original_params['main.max_time'] = 0.2
    original_params['regrid.fixed_grid_time'] = 0.05
    original_params['main.plot_interval'] = 25 #50
    original_params['main.checkpoint_interval'] = 200


    uniformPorosity = {'parameters.rayleighTemp':'1e5',
    'main.output_folder':'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputUniformPorosity1.0'}
    uniformLowerPorosity = {'bc.porosityLoVal':'0.1 0.1',
    'main.output_folder':'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputUniformPorosity0.1',
    'parameters.rayleighTemp':'1e6'}
    variablePorosity = {'parameters.porosityFunction':1,
                        'bc.porosityLoVal':'0.5 0.5',
                        'bc.enthalpyHiVal':'0.0 6.0',
                        'bc.enthalpyLoVal':'0.0 8.05',
                        'parameters.problem_type':14,
                        'parameters.rayleighTemp':'1e5',
                        'main.output_folder':'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputVariablePorosity'}   

    timeDependentPorosity = {'porosityTimescale':2/original_params['main.max_time'],
                            'parameters.porosityFunction':4,
                            'parameters.problem_type':14,
                            'parameters.rayleighTemp':'2e5',
                            'bc.porosityHiVal':'0.1 0.1',
                            'main.output_folder':'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputTimeDependentPorosity'}

    timeDepSimpleAdvection = dict(timeDependentPorosity)
    timeDepSimpleAdvection['main.advectionMethod'] = 2
    timeDepSimpleAdvection['main.output_folder'] = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputTimeDependentPorositySimpleAdvection'

    nonTimeDep = dict(timeDependentPorosity)
    nonTimeDep['main.porosityEndTime'] = original_params['regrid.fixed_grid_time']*0.9
    nonTimeDep['main.output_folder'] = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputTimeDependentBeforeRegrid'

    coupledPorosity = dict(uniformPorosity)
    coupledPorosity['main.output_folder'] = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputCoupledPr2'
    coupledPorosity['parameters.stefan'] = 0.5 # try and make lots of ice
    coupledPorosity['parameters.compositionRatio'] = 20
    coupledPorosity['parameters.rayleighTemp'] = 5e5
    coupledPorosity['parameters.prandtl'] = 2
    coupledPorosity['bc.enthalpyHiVal']='0.0 0.8'
    coupledPorosity['bc.enthalpyLoVal']='0.0 1.7'

    coupledPorosityMoreMush = dict(coupledPorosity)
    coupledPorosityMoreMush['main.output_folder'] = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputCoupledPr1'
    coupledPorosityMoreMush['parameters.prandtl'] = 1

    coupledPorosityLowerC = dict(coupledPorosityMoreMush)
    coupledPorosityLowerC['main.output_folder'] = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputCoupledCR5'
    coupledPorosityLowerC['parameters.compositionRatio'] = 5
    #coupledPorosityLowerC['parameters.darcy'] = 1e-2

    coupledPorosityLowC = dict(coupledPorosityMoreMush)
    coupledPorosityLowC['main.output_folder'] = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputCoupledCR3'
    coupledPorosityLowC['parameters.compositionRatio'] = 3

    threeLevels = dict(coupledPorosityMoreMush)
    threeLevels['main.output_folder'] = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputThreeLevel'
    threeLevels['parameters.compositionRatio'] = 5
    threeLevels['main.max_level'] = 2
    threeLevels['regrid.plume_vel'] = 1.0
    #threeLevels.pop('regrid.fixed_grid_time') # let's do automatic refinement

    threeLevels = dict(coupledPorosityMoreMush)
    threeLevels['main.output_folder'] = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputThreeLevel'
    threeLevels['parameters.compositionRatio'] = 5
    threeLevels['main.max_level'] = 2
    threeLevels['regrid.plume_vel'] = 1.0

    runsToDo = [uniformPorosity, uniformLowerPorosity, variablePorosity, timeDependentPorosity, nonTimeDep, 
    coupledPorosity, coupledPorosityMoreMush, coupledPorosityLowerC, coupledPorosityLowC, threeLevels]
    #runsToDo = [timeDependentPorosity]

    for new_params in runsToDo:

        params = dict(original_params)

        for key in new_params.keys():
            params[key] = new_params[key]

        if params['main.output_folder'] == '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/outputThreeLevel':
            params.pop('regrid.fixed_grid_time')
            params.pop('main.vel_refine_thresh')


        # Debugging options:
        if debug:
            params['main.plot_interval']=1
            params['main.subcycled_plots']=1
            params['main.regrid_plots']=1
            params['main.writePressureInitFields']=1
      

        output_folder = params['main.output_folder']

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

            #time.sleep(0.1) # wait long enough to ensure folder is made
            params_file = os.path.join(output_folder, 'inputs')
            print('Params file: ' + params_file)
            writeParamsFile(params_file, params)

            # Do the run
            cmd = 'cd ' + output_folder + '; ' + execFile + ' inputs'

            os.system(cmd)

        

if __name__ == "__main__":
    print('Running script')
    main(sys.argv[1:])
                 


