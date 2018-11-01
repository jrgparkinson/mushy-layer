import re
import os
import time
import math
import sys

def getBaseOutputDir():
    #base_output_dir = '/home/parkinsonjl/mushy-layer/test/output/'
    base_output_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/Test/'

    if not os.path.exists(base_output_dir):
        os.makedirs(base_output_dir)

    return base_output_dir


def getMatlabBaseCommand():
    parDir = os.path.abspath(os.pardir)
    matlabFolder = os.path.join(parDir, 'matlab', 'MethodsPaper')
    matlab_command = 'cd ' + matlabFolder + '; \n \n matlab -nodisplay -nosplash -nodesktop -r'

    return matlab_command




def getExecName():
    mushyLayerDir = os.path.dirname(os.path.dirname(__file__))
    exec_dir = os.path.join(mushyLayerDir, 'execSubcycle')

    # Get files in exec dir starting with mushyLayer2d and ending with ex
    possible_exec_files = [f for f in os.listdir(exec_dir) if f[:12] == 'mushyLayer2d' and f[-2:] == 'ex']

    if len(possible_exec_files) == 0:
        print('Cannot find any executable files - have you compiled the code?')
        sys.exit(0)

    # Choose optimised execs over DEBUG execs as they will be quicker
    opthigh_exec = ''
    opt_exec = ''

    for f in possible_exec_files:
        if 'OPTHIGH' in f:
            opthigh_exec = f
        elif 'OPT' in f:
            opt_exec = f

    if opthigh_exec:
        exec_file = opthigh_exec
    elif opt_exec:
        exec_file = opt_exec
    else:
        exec_file = possible_exec_files[1]

    # Can also specify the file manually
    #exec = 'mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex'


    # Sanity check
    if not os.path.exists(os.path.join(exec_dir, exec_file)):
        print('Executable ' + exec_file +' not found in directory ' + exec_dir)
        sys.exit(0)

    return exec_file

def constructRunName(params):
    # run_name = 'CR' + str(params['parameters.compositionRatio']) + 'RaC' + str(params['parameters.rayleighComp'])

    longToShortName = {'parameters.compositionRatio': 'CR',
                       'parameters.rayleighComp': 'RaC',
                       'parameters.lewis': 'Le',
                       'parameters.darcy': 'Da',
                       'parameters.nonDimReluctance': 'R',
                       'parameters.bottomEnthalpy': 'HBottom'}

    paramsToInclude = ['parameters.compositionRatio',
                       'parameters.rayleighComp',
                       'parameters.lewis',
                       'parameters.permeabilityFunction']
                       
                       
                       #'parameters.bottomEnthalpy',
                       #'parameters.nonDimReluctance',
    
    paramFormat = {'parameters.compositionRatio': '%1.3f',
                       'parameters.lewis': '%d',
                       'parameters.darcy': '%.1e',
                       'parameters.nonDimReluctance': '%.1e',
                       'parameters.bottomEnthalpy': '%1.2f'}
                       
    if 'parameters.rayleighComp' in params:
        if float(params['parameters.rayleighComp']) > 1000:
            paramFormat['parameters.rayleighComp'] = '%.1e'
        else:
            paramFormat['parameters.rayleighComp'] = '%d'
                
             
    permeabilityFunctions = ['UniformPermeability',
                             'ChiCubedPermeability',
                             'KozenyPermeability',
                             'LogPermeability',
                             'CustomPermeability']

    run_name = ''

    for p in paramsToInclude:
        if p in longToShortName.keys() and p in params:
            run_name = run_name + longToShortName[p] + paramFormat[p] % float(params[p])
        elif p == 'parameters.permeabilityFunction' and p in params:
            permFunc = int(params[p])
            run_name = run_name + permeabilityFunctions[permFunc]
        else:
            # Don't know how to handle this. Just do nothing.
            pass
            
    if float(params['parameters.darcy']) > 0:
            run_name = run_name + longToShortName['parameters.darcy'] + paramFormat['parameters.darcy'] % float(params['parameters.darcy'])
            run_name = run_name + longToShortName['parameters.nonDimReluctance'] + paramFormat['parameters.nonDimReluctance'] % float(params['parameters.nonDimReluctance'])

    # params which don't follow the same format

    # if we're using a sponge
    if ('main.spongeHeight' in params 
        and 'main.spongeRelaxationCoeff' in params
        and float(params['main.spongeHeight']) > 0.0
        and float(params['main.spongeRelaxationCoeff']) > 0.0):
        run_name = run_name + '_sponge_'

    # grid resolution
    cells =params['main.num_cells'].split(' ')
    grid_res = int(cells[0])
    run_name = run_name + "pts" + str(grid_res)

        
    return run_name

# Load up an inputs file and parse it into a dictionary
def readInputs(inputs_file):

	params = {}

        # Read in the file
        filedata = None
        with open(inputs_file, 'r') as f:
            filedata = f.readlines()

        for line in filedata:
            # print(line)
            # Ignore lines that are commented out
            if line.startswith('#'):
                continue

            # Remove anything after a #
            line = re.sub('#[^\n]*[\n]', '', line)
            # print(line)

            parts = line.split('=')
            if len(parts) > 1:
                key = parts[0].strip()
                val = parts[1].strip()
                params[key] = val

       # print(params)
        return params

    

def writeParamsFile(location, params, ignoreList = [], doSort=True):
    output_file = ''


    
    keyList = params.keys()
    if doSort:
        keyList.sort()    
    
    for key in keyList:
        if key not in ignoreList:
            
            if isinstance(params[key], list):
                keyVal = ' '.join([str(a) for a in params[key]])
                
            else:
                keyVal = str(params[key])



            output_file = output_file + '\n' + \
                key + '=' + keyVal

    with open(location, 'w') as f:
        f.write(output_file)
        
        
def hasReachedSteadyState(folder):
    timeTableFile = os.path.join(folder, 'time.table.0')
    #print(timeTableFile)
    if os.path.exists(timeTableFile):
        return True
        

def timeSinceFolderUpdated(directory):
    smallest_t = 1e100
    for filename in os.listdir(directory):
        this_time_diff = timeSinceFileUpdated(os.path.join(directory, filename))
        smallest_t = min(smallest_t, this_time_diff)
        
    return smallest_t
        
def timeSinceFileUpdated(filename):
    current_t = time.time()
    t = os.path.getmtime(filename)
        
    this_time_diff = abs(current_t - t)
    
    return this_time_diff
    
def getRestartFile(most_recent_path):

    restartFile = ''
    
    chkFiles = [f for f in os.listdir(most_recent_path) if len(f) > 4 and f[:3] == 'chk']
    
    #print('Chk files in ' + most_recent_path + ': ')
    #print(chkFiles)

    if len(chkFiles) > 0:
        restartFile = chkFiles[-1]
        
    print('Found ' + str(len(chkFiles)) + ' checkpoint files, most recent: ' + restartFile)

    return restartFile

def getFinalPlotFile(directory):
    filesDir = [f for f in os.listdir(directory)  if (os.path.isfile(os.path.join(directory, f))) ]
    pltFiles = []
    #print(filesDir)
    for f in filesDir:
        #print(f)
        #print(f[-5:])
        if len(f) > 5 and not 'chk' in f and f[-5:] == '.hdf5':
            pltFiles.append(f)

    pltFiles = sorted(pltFiles)

    if pltFiles:
        return pltFiles[-1]
    else:
        return None


def getFinalChkFile(directory):
    filesDir = [f for f in os.listdir(directory)  if (os.path.isfile(os.path.join(directory, f))) ]
    pltFiles = []
    #print(filesDir)
    for f in filesDir:
        #print(f)
        #print(f[-5:])
        if len(f) > 5 and 'chk' in f and f[-5:] == '.hdf5':
            pltFiles.append(f)

    pltFiles = sorted(pltFiles)

    if pltFiles:
        return pltFiles[-1]
    else:
        return None

def isPowerOfTwo(n):
    test = math.log(n)/math.log(2)
    
    if round(test) == test:
        return True
    else:
        return False

def str2arr(string):
    #print(string)
    parts = string.split(' ')
    array = [int(i) for i in parts]
    return array

def arr2str(array):
    str_array = [str(a) for a in array]
    string = ' '.join(str_array)
    return string