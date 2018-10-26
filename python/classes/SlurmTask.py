import os
import math
import time
import subprocess
import re

class SlurmTask:

    MAX_TASKS_PER_NODE = 8
    
    def __init__(self, folder, jobname, execFile, num_proc=1, timeLimit = 7.0, memoryLimit=5000, numNodes=0, mpirun=True):

        self.jobname = jobname
        self.folder = folder
        self.num_proc = num_proc
        
        self.execFile = execFile
        self.timeLimit = timeLimit # Measured in days
        self.memoryLimit = memoryLimit # in MB
        self.mpirun=mpirun

        self.preprocessCommand = ''
        self.dependency = []
        self.jobID = -1
        self.customCommand = ''

        # Default setup: all tasks on one node
        # This will fail for large numbers of processors - need to implement 
        # some clever load balancing

        # Always fix this
        self.cpu_per_task = 1

        # Number of nodes to use
        if numNodes > 0:
            self.num_nodes = numNodes
        else:
            # Determine number of nodes based on number of processors and MAX_TASKS_PER_NODE
            self.num_nodes = math.ceil(float(self.num_proc)/float(self.MAX_TASKS_PER_NODE))

        # Distribute tasks onto nodes
        tasks_per_node, remainder = divmod(self.num_proc, self.num_nodes)
        if remainder == 0:
            self.tasks_per_node = tasks_per_node
        else:
            # Add an extra task per node to ensure we have enough processors

            self.tasks_per_node = tasks_per_node + 1
            print('SlurmTask - Warning - number of processors doesn\'t split evenly onto number of nodes. \n')

        # Each node usually has two sockets
        self.tasks_per_socket = math.ceil(self.tasks_per_node/2) #self.tasks_per_node #

    def setPreprocess(self, cmd):
        self.preprocessCommand = cmd

    def setDependency(self, dep):
        if isinstance(dep, list):
            self.dependency = [str(a) for a in dep]
        else:
            self.dependency = [str(dep)]


    def setCustomCommand(self, cmd):
        self.customCommand = cmd

    def writeSlurmFile(self, runFileName = 'run.sh', inputsFileName = 'inputs'):
    	
        #self.execFile = self.execFile.replace('OPT.', '')
        

    	fh = open(self.getRunFile(runFileName), 'w+')

        days,  remainder = divmod(self.timeLimit, 1)
        hours, remainder = divmod(remainder*24, 1)
        mins,  remainder = divmod(remainder*60, 1)
        secs = remainder*60

    	time_string = '%d-%d:%d:%d' % (days,hours,mins,secs)

        mpiStr = ''
        if self.mpirun:
            mpiStr = 'mpirun -np ' + str(self.num_proc)

        if len(self.dependency) > 0:
            dependencies = ','.join(self.dependency)
            depStr = '#SBATCH --dependency=afterok:' + dependencies + ' \n'
        else:
            depStr = ''

    	file_contents = ['#!/bin/bash \n',
		'# Set your minimum acceptable walltime, format: day-hours:minutes:seconds \n',
		'#SBATCH --time=' + time_string + '\n',
		'#SBATCH --partition=shared,priority-ocean' + '\n',
		'# Set name of job shown in squeue' + '\n',
		'#SBATCH --job-name ' + self.jobname + '\n',
		'# Request CPU resources' + '\n',
		'#SBATCH --ntasks=' + str(int(self.num_proc) ) + '                  # Number of MPI ranks' + '\n',
		'#SBATCH --cpus-per-task=' + str(int(self.cpu_per_task) ) + '            # Number of cores per MPI rank ' + '\n',
		'#SBATCH --nodes=' + str(int(self.num_nodes)) + '                    # Number of nodes' + '\n',
		'#SBATCH --ntasks-per-node=' + str(int(self.tasks_per_node) ) + '         # How many tasks on each node' + '\n',
		'#SBATCH --ntasks-per-socket=' + str(int(self.tasks_per_socket) ) + '        # How many tasks on each CPU or socket (not sure what this really means)' + '\n',
		'#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically on nodes and sockets' + '\n',
		'# Memory usage (MB)' + '\n',
		'#SBATCH --mem-per-cpu=' + str(self.memoryLimit) + '\n',
		'#SBATCH --output=' + os.path.join(self.folder, 'sbatch.out') + '   # Standard output and error log' + '\n',
        depStr, 
		'cd ' + self.folder +  '; \n \n ']

        if self.customCommand:
            file_contents = [file_contents, self.customCommand]

        else:
            file_contents = [file_contents, self.preprocessCommand + ' ' + mpiStr + ' ' + self.execFile + ' ' + os.path.join(self.folder, inputsFileName)]
		

    	fh.writelines(file_contents)

    	fh.close()


    def runTask(self, runFileName='run.sh'):
    	cmd = 'cd ' + self.folder + ' ; sbatch ' + self.getRunFile(runFileName)
    	
        #os.system(cmd)
        #result = subprocess.run(cmd, stdout=subprocess.PIPE)
        result = subprocess.check_output(    cmd,  stderr=subprocess.STDOUT,  shell=True)

        res_str = result

        print(result)

        # Submitted batch job 223086
        res = re.findall('Submitted batch job (\d+)', res_str)

        #print(res_str)
        #print(res)

        if res:
            self.jobID = int(res[0])

        # Make a file in the output folder stating the slurm job id
        F = open(os.path.join(self.folder, 'jobid'), 'w')
        F.write(str(self.jobID))
        F.close()

        # Pause briefly in case submitting lots of slurm jobs
        time.sleep(0.5)

    def getRunFile(self, runFileName):
        return os.path.join(self.folder, runFileName)

   