from timeout import timeout
import os, sys
import subprocess
from subprocess import Popen
import paramiko
import time
import re
from datetime import date, datetime

# Call this function from the command line like:
#  python restartRun.py restart directory
# where directory is the one to copy to the cloud, and time delay (minutes)
# determines how old folders must be to be copied. I.e. folders modified more
# recently than now - timeDelay will not be copied. This is to stop messing
# with runs that are currently in progress. 
#
# Example usage: 
#  python restartRun.py restart mushyLayerLowC-upperBranch-insulating/CR1.1RaC600Le200ChiCubedPermeabilitypts128-0

# Copy all folders in dir which haven't been modified
# since timeDelay minutes ago to the cloud
def restart(folder):
    # get full folder
    cwd = os.getcwd()
    full_dir = os.path.join(cwd, folder)
    restart_dir = os.path.join(full_dir, '/restart')
    inputs_file = os.path.join(full_dir, 'inputs')
    restart_inputs_file = os.path.join(full_dir, 'inputs-restart')

    # get number of processors:
    # get all pout files in the run directory, and count them - this is the number of processors
    poutfiles = [f for f in listdir(full_dir) if (isfile(join(full_dir, f)) and 'pout.' in f)]
    num_proc = len(onlyfiles)

    # create restart inputs and folder
    os.mkdir(restart)

    # get previous inputs 
    params = readInputs(inputs_file)

    # restart from previous final chk file
    chkfiles = [f for f in listdir(full_dir) if (isfile(join(full_dir, f)) and 'chk' in f)]
    finalChkFile = chkfiles[-1]
    params['main.restart_file'] = finalChkFile

    # write out new inputs
    output_str = ''
    param_keys = params.keys()
    sorted(param_keys)
    for key in param_keys:
        output_str = output_str + '\n' + \
                key + '=' + str(params[key])

    with open(restart_inputs_file, 'w') as f:
        f.write(output_str)


    # execute program
    machine_specific = MachineSpecific()
    mushyLayerExec = machine_specific.get_home_dir() + \ 
        '/convection-in-sea-ice/MushyLayer/execSubcycle/' + \
        machine_specific.get_program_name()
    command = 'mpirun -np ' + str(num_proc) + ' ' + mushyLayerExec + ' ' + restart_inputs_file
    

if __name__ == '__main__':
    import argparse
    function_map = { 
        'restart': restart,
        'restartRun': restart,
    }
    
    parser = argparse.ArgumentParser()
    parser.add_argument( 'command', nargs=1 )
    parser.add_argument( 'folder', nargs='+' )
    args= parser.parse_args()
    function = function_map[args.command[0]]
    function( args.folder[0] )
    
