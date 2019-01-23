# Function to run a convergence test
import os
import sys
from colorama import Fore, Style
import math
from subprocess import Popen

parDir = os.path.abspath(os.pardir)
pythonDir = os.path.join(parDir, 'python')
sys.path.append(pythonDir)

from MushyLayerRunSimple import MushyLayerRunSimple
from mushyLayerRunUtils import construct_run_name, read_inputs, write_inputs, get_restart_file, get_executable_name
from SlurmTask import SlurmTask
# Params is a vector of params dicts



def AMRConvergenceTest(params, full_output_dir, physicalProblem, Nzs, num_procs = [1], numRestarts=0, analysis_command=''):

    os.environ["CH_TIMER"] = "1"
    subcycled = True # subcycled code

    dependencies = []

    while len(num_procs) < len(Nzs):
        num_procs.append(num_procs[-1])

    param_i = -1

    for p in params:
        param_i = param_i + 1

        num_proc = p['num_proc']
                   
        
        # Don't repeat runs, unless we're restarting
        
        #cwd = os.getcwd()
        run_name = p['main.plot_prefix']
        test_name = run_name + '-0'
        full_path = os.path.join(full_output_dir, test_name)
        
        allowRestarts = False
        
        print(full_path)
        if os.path.exists(full_path):
        
            print(Fore.YELLOW + '    Run already done' + Fore.RESET)
            
            if numRestarts > 0:
            
                print('    ' + str(numRestarts) + ' restart(s) allowed \n ')
                # Find most recent folder and restart from there
                
                i = 0
                while os.path.exists(full_path):
                    i = i + 1
                    test_name = run_name + '-' + str(i)
                    full_path = os.path.join(full_output_dir, test_name)
                
                if i > numRestarts:
                    continue
                    
                most_recent_path = os.path.join(full_output_dir, run_name + '-' + str(i-1))  
                # Now get the restart file
                restartFile = get_restart_file(most_recent_path)
                if restartFile:
                    p['main.restart_file'] = os.path.join(most_recent_path, restartFile)
                    allowRestarts = True
                    
                    print('    Set restart file: ' + p['main.restart_file'])
                
                
                
            else:
            
                continue
                
        else:
            print(Fore.YELLOW + '    No previous runs' + Fore.RESET)
    
      
        print(Fore.GREEN + '    **Do run**' + Fore.RESET)

        # Don't need output dir or exec dir here, MushyLayerRun will fill these in
        s = SlurmTask('', p['concise_run_name'], '', num_proc)

        ml_run = MushyLayerRunSimple(full_output_dir, num_proc, p, s, allowRestarts, get_executable_name())
        ml_run.single_run(run_name)

        dependencies.append(s.jobID)

        print('==================')

        
    return dependencies




    
        

