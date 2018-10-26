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
    params_file = 'inputs'
    num_refine = 1

    original_params = readInputs(base_dir + params_file)

    for i in range(0, num_refine+1):

        params = dict(original_params)

        refinement = int(pow(2, i))

        # Params that need to be scaled: 
        # block factor, grid size, tag_buffer_size, grid_buffer_size
        paramsToScale = ['main.block_factor', 'main.max_grid_size', 'main.tag_buffer_size']
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
      

        output_name = params['main.output_folder'] + '-Nx' + str(params['main.num_cells'][0])
        output_folder = base_dir + output_name

        params['main.output_folder'] = output_folder

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
                 


