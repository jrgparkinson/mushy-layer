import re
import os
import time
import math
import sys
import subprocess
import socket


def get_base_output_dir():
    #base_output_dir = '/home/parkinsonjl/mushy-layer/test/output/'
    #base_output_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/TestDiffusiveTimescale/'
    base_output_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/TestFinal/'


    if not os.path.exists(base_output_dir):
        os.makedirs(base_output_dir)

    return base_output_dir

def add_params(defaultParams, extra_params):
    """ Add params from extra_params to defaultParams """

    for k, v in extra_params.iteritems():
        defaultParams[k] = v

    return defaultParams


def get_matlab_base_command():
    parent_dir = os.path.abspath(os.pardir)
    matlab_folder = os.path.join(parent_dir, 'matlab', 'MethodsPaper')
    matlab_command = 'cd ' + matlab_folder + '; \n \n matlab -nodisplay -nosplash -nodesktop -r'

    return matlab_command


def get_current_vcs_revision():
    # For mercurial:
    # pipe = subprocess.Popen(
    #    ["hg", "identify", "--num"],
    #    stdout=subprocess.PIPE
    # )
    # repoVersion = pipe.stdout.read()
    # repoVersion = repoVersion.replace('\n', '')  # clean up a bit

    # For git:
    repo_version = subprocess.check_output(['git', 'rev-parse', 'HEAD'])

    return repo_version

def get_mushy_layer_dir():

    #if 'MUSHY_LAYER_DIR' in os
    this_file_path = os.path.realpath(__file__)
    par_dir = os.path.dirname(this_file_path)
    mushy_layer_dir = os.path.dirname(par_dir)

    #dir = os.environ['MUSHY_LAYER_DIR']
    #dir = '/home/parkinsonj/mushy-layer/'

    return mushy_layer_dir


def get_executable(base_name, dim=2):
    """ Given some program name, and assuming compilation options (mpiCC, gfortran, OPT, MPI etc.)
     return the correct executable for this system"""

    host = ''

    # If we're on gyre, append set host string to .GYRE
    hostname = socket.gethostname()
    if 'gyre' in hostname:
        host = '.GYRE'

    executable_name = '%s%dd.Linux.64.mpiCC.gfortran.OPT.MPI%s.ex' % (base_name, dim, host)

    return executable_name


def get_executable_name(exec_dir='', exec_name='mushyLayer2d'):
    """ Get the executable in MUSHY_LAYER_DIR/execSubcycle/ which contains mushyLayer2d,
     prioritising OPT over DEBUG compilations.
     Independent of the architecture (i.e. will find for different C++ compilers, fortran versions etc.) """

    if not exec_dir:
        mushy_layer_dir = get_mushy_layer_dir() # previously: os.path.dirname(os.path.dirname(__file__))
        exec_dir = os.path.join(mushy_layer_dir, 'execSubcycle')

    # Get files in exec dir starting with mushyLayer2d and ending with ex
    possible_exec_files = [f for f in os.listdir(exec_dir) if f[:len(exec_name)] == exec_name and f[-2:] == 'ex']

    if len(possible_exec_files) == 0:
        print('Searching %s' % exec_dir)
        print('Cannot find any executable files \'%s\' - have you compiled the code?' % exec_name)
        sys.exit(0)

    init_possible_exec_files = possible_exec_files

    # Remove/keep .GYRE in list as appropriate
    if 'gyre' in socket.gethostname():
        possible_exec_files = [f for f in possible_exec_files if 'GYRE' in f]
    else:
        possible_exec_files = [f for f in possible_exec_files if 'GYRE' not in f]

    if len(possible_exec_files) == 0:
        print('Cannot find any executable files that match the current system.')
        print('Executable founds: ' + '\n'.join(init_possible_exec_files))
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
    # exec = 'mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex'

    # Sanity check
    if not os.path.exists(os.path.join(exec_dir, exec_file)):
        print('Executable ' + exec_file +' not found in directory ' + exec_dir)
        sys.exit(0)

    return exec_file


def construct_run_name(params):
    # run_name = 'CR' + str(params['parameters.compositionRatio']) + 'RaC' + str(params['parameters.rayleighComp'])

    long_to_short_name = {'parameters.compositionRatio': 'CR',
                          'parameters.rayleighComp': 'RaC',
                          'parameters.lewis': 'Le',
                          'parameters.darcy': 'Da',
                          'parameters.nonDimReluctance': 'R',
                          'parameters.bottomEnthalpy': 'HBottom'}

    params_to_include = ['parameters.compositionRatio',
                         'parameters.rayleighComp',
                         'parameters.lewis',
                         'parameters.permeabilityFunction']

    param_format = {'parameters.compositionRatio': '%1.3f',
                    'parameters.lewis': '%d',
                    'parameters.darcy': '%.1e',
                    'parameters.nonDimReluctance': '%.1e',
                    'parameters.bottomEnthalpy': '%1.2f'}
                       
    if 'parameters.rayleighComp' in params:
        if float(params['parameters.rayleighComp']) > 1000:
            param_format['parameters.rayleighComp'] = '%.1e'
        else:
            param_format['parameters.rayleighComp'] = '%d'

    permeability_functions = ['UniformPermeability',
                             'ChiCubedPermeability',
                             'KozenyPermeability',
                             'LogPermeability',
                             'CustomPermeability']

    run_name = ''

    for p in params_to_include:
        if p in long_to_short_name.keys() and p in params:
            run_name = run_name + long_to_short_name[p] + param_format[p] % float(params[p])
        elif p == 'parameters.permeabilityFunction' and p in params:
            perm_func = int(params[p])
            run_name = run_name + permeability_functions[perm_func]
        else:
            # Don't know how to handle this. Just do nothing.
            pass
            
    if float(params['parameters.darcy']) > 0:
            run_name = run_name + long_to_short_name['parameters.darcy'] + param_format['parameters.darcy'] % float(params['parameters.darcy'])
            run_name = run_name + long_to_short_name['parameters.nonDimReluctance'] + param_format['parameters.nonDimReluctance'] % float(params['parameters.nonDimReluctance'])

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


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False
    
def isint(value):
    try:
        int(value)
        return True
    except ValueError:
        return False
    
def read_inputs(inputs_file):
    """ Load up an inputs file and parse it into a dictionary """

    params = {}

    # Read in the file
    file_data = None
    with open(inputs_file, 'r') as f:
        file_data = f.readlines()

    for line in file_data:
        # print(line)
        # Ignore lines that are commented out
        if line.startswith('#'):
            continue

        # Remove anything after a #
        line = re.sub('#[^\n]*[\n]', '', line)

        parts = re.findall('^([^=]*)=(.*)$', line)

        # parts = line.split('=')
        # if len(parts) > 1:
        if parts:
            match = parts[0]
            key = match[0].strip()
            val = match[1].strip()

            # Convert to float/int as appropriate
            if isint(val):
                val = int(val)
            elif isfloat(val):
                val = float(val)

            params[key] = val
            
    # print(params)
    return params


def write_inputs(location, params, ignore_list = None, doSort=True):
    output_file = ''

    key_list = params.keys()
    if doSort:
        key_list.sort()
    
    for key in key_list:
        if not ignore_list or key not in ignore_list:
            
            if isinstance(params[key], list):
                key_val = ' '.join([str(a) for a in params[key]])
                
            else:
                key_val = str(params[key])

            output_file = output_file + '\n' + \
                key + '=' + key_val

    with open(location, 'w') as f:
        f.write(output_file)


def has_reached_steady_state(folder):
    time_table_file = os.path.join(folder, 'time.table.0')

    if os.path.exists(time_table_file):
        return True
        

def time_since_folder_updated(directory):
    smallest_t = 1e100
    for filename in os.listdir(directory):
        this_time_diff = time_since_file_updated(os.path.join(directory, filename))
        smallest_t = min(smallest_t, this_time_diff)
        
    return smallest_t


def time_since_file_updated(filename):
    current_t = time.time()
    t = os.path.getmtime(filename)
        
    this_time_diff = abs(current_t - t)
    
    return this_time_diff
    
def get_restart_file(most_recent_path):

    restart_file = ''
    
    chk_files = [f for f in os.listdir(most_recent_path) if len(f) > 4 and f[:3] == 'chk']
    
    # print('Chk files in ' + most_recent_path + ': ')
    # print(chk_files)

    if len(chk_files) > 0:
        restart_file = chk_files[-1]
        
    print('Found ' + str(len(chk_files)) + ' checkpoint files, most recent: ' + restart_file)

    return restart_file


def get_final_plot_file(directory):
    files_dir = [f for f in os.listdir(directory)  if (os.path.isfile(os.path.join(directory, f))) ]
    plt_files = []
    # print(files_dir)
    for f in files_dir:

        # Match if it doesn't start with chk, and finished like 01234.2d.hdf5
        pattern = '^(?!chk.*$).*\d+\.\dd\.hdf5'
        if re.match(pattern, f):
        #if len(f) > 5 and not 'chk' in f and f[-5:] == '.hdf5':
            plt_files.append(f)

    plt_files = sorted(plt_files)

    if plt_files:
        return plt_files[-1]
    else:
        return None


def get_final_chk_file(directory):
    files_dir = [f for f in os.listdir(directory)  if (os.path.isfile(os.path.join(directory, f))) ]
    plt_files = []

    for f in files_dir:

        # Match if file starts with chk and ends in [frame number].2d.hdf5
        pattern = '^chk.*\d+\.\dd\.hdf5'

        if re.match(pattern, f):
        #if len(f) > 5 and 'chk' in f and f[-5:] == '.hdf5':
            plt_files.append(f)

    plt_files = sorted(plt_files)

    if plt_files:
        return plt_files[-1]
    else:
        return None


def is_power_of_two(n):
    test = math.log(n)/math.log(2)
    
    if round(test) == test:
        return True
    else:
        return False


def string_to_array(string, conversion = lambda x:int(x)):

    if isinstance(string, list):
        return string

    parts = string.split(' ')
    array = [conversion(i) for i in parts]
    return array


def array_to_string(array):
    str_array = [str(a) for a in array]
    string = ' '.join(str_array)
    return string



def check_exec_exists(exec_dir, exec_name):
    base_name = os.path.join(exec_dir, exec_name)
    exec_name = get_executable(base_name, 2)

    # If executable exists, we're all good
    if os.path.exists(exec_name):
        return True

    # If executable doesn't exist, warn the user and ask if they want us to try and build it
    print('Error - executable does not exist: %s' % exec_name)
    response = raw_input('Try and build this executable? (y/n) ')

    if response == 'y':
        cmd = 'cd %s; make all; cd -;' % exec_dir
        os.system(cmd)

        if os.path.exists(exec_name):
            return True

        else:
            print('Unable to build the executable, please try yourself manually before running the test problems.')
            return False



    print('Please build the executable yourself before running the test problems.')

    return False


