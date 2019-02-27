import os
import math
import time
import subprocess
import re


class BatchJob:
    '''
    This is setup for SLURM, but should be fairly straightforward to modify for other queuing systems
    '''

    MAX_TASKS_PER_NODE = 16
    MAX_MEMORY = 128000.0  # MB

    postprocess_command = ''

    def __init__(self, folder, jobname, exec_file, num_proc=1, time_limit=7.0, memory_limit=4000, num_nodes=0,
                 mpi_run=True):

        self.jobname = jobname
        self.folder = folder
        self.num_proc = num_proc

        self.exec_file = exec_file
        self.time_limit = time_limit  # Measured in days
        self.memory_limit = memory_limit  # in MB
        self.mpi_run = mpi_run

        self.preprocess_command = ''
        self.dependency = []
        self.job_id = -1
        self.custom_command = ''

        self.partitions = ['shared', 'priority-ocean']
        self.exclude = None

        # Default setup: all tasks on one node
        # This will fail for large numbers of processors - need to implement 
        # some clever load balancing

        # Always fix this
        self.cpu_per_task = 1

        # Number of nodes to use
        if num_nodes > 0:
            self.num_nodes = num_nodes
        else:
            # Determine number of nodes based on number of processors and MAX_TASKS_PER_NODE
            self.num_nodes = math.ceil(float(self.num_proc) / float(self.MAX_TASKS_PER_NODE))

        # Distribute tasks onto nodes
        tasks_per_node, remainder = divmod(self.num_proc, self.num_nodes)
        if remainder == 0:
            self.tasks_per_node = tasks_per_node
        else:
            # Add an extra task per node to ensure we have enough processors

            self.tasks_per_node = tasks_per_node + 1
            print('BatchJob - Warning - number of processors doesn\'t split evenly onto number of nodes. \n')

        # Each node usually has two sockets
        self.tasks_per_socket = math.ceil(self.tasks_per_node / 2)  # self.tasks_per_node #

        # Make sure tasks per node * memory < MAX_MEMORY
        if self.tasks_per_node * self.memory_limit > self.MAX_MEMORY:
            self.memory_limit = 0.9 * self.MAX_MEMORY / float(self.tasks_per_node)
            print('Reducing memory limit to %1.0f GB per node to avoid exceeding node memory limit. \n' % (
                    self.memory_limit / 1000))

    def set_preprocess(self, cmd):
        self.preprocess_command = cmd

    def set_post_process(self, cmd):
        self.postprocess_command = cmd

    def set_dependency(self, dep):
        if isinstance(dep, list):
            self.dependency = [str(a) for a in dep]
        else:
            self.dependency = [str(dep)]

    def set_custom_command(self, cmd):
        self.custom_command = cmd

    def set_partition_exclude(self, partitions=None, exclude=None):
        self.partitions = partitions
        self.exclude = exclude

    def set_exec_file(self, exec_file):
        self.exec_file = exec_file

    def write_batch_file(self, runFileName='run.sh', inputs_file_name='inputs'):

        fh = open(self.get_run_file(runFileName), 'w+')

        # Get time limit in the format days-hours:mins:secs
        days, remainder = divmod(self.time_limit, 1)
        hours, remainder = divmod(remainder * 24, 1)
        mins, remainder = divmod(remainder * 60, 1)
        secs = remainder * 60

        time_string = '%d-%d:%d:%d' % (days, hours, mins, secs)

        mpi_str = ''
        if self.mpi_run:
            mpi_str = 'mpirun -np ' + str(self.num_proc)

        # If you don't have SLURM, the dependency command will be different
        if len(self.dependency) > 0:
            dependencies = ':'.join(self.dependency)
            dependency_str = '#SBATCH --dependency=afterok:' + dependencies + ' \n'
        else:
            dependency_str = ''

        if self.partitions and len(self.partitions) > 0:
            partitions_str = '#SBATCH --partition=' + ','.join(self.partitions) + '\n'
        else:
            partitions_str = ''

        if self.exclude and len(self.exclude) > 0:
            exclude_str = '#SBATCH -x ' + ','.join(self.exclude) + '\n'
        else:
            exclude_str = ''

        # If you don't have slurm, this should be different
        slurm_header = ['# Set your minimum acceptable walltime, format: day-hours:minutes:seconds \n',
                         '#SBATCH --time=' + time_string + '\n', partitions_str, exclude_str,
                         '# Set name of job shown in squeue' + '\n', '#SBATCH --job-name ' + self.jobname + '\n',
                         '# Request CPU resources' + '\n', '#SBATCH --ntasks=' + str(int(self.num_proc)) + '                  # Number of MPI ranks' + '\n',
                         '#SBATCH --cpus-per-task=' + str(int(self.cpu_per_task)) + '            # Number of cores per MPI rank ' + '\n',
                         '#SBATCH --nodes=' + str(int(self.num_nodes)) + '                    # Number of nodes' + '\n',
                         '#SBATCH --ntasks-per-node=' + str(int(self.tasks_per_node)) + '         # How many tasks on each node' + '\n',
                         '#SBATCH --ntasks-per-socket=' + str(int(self.tasks_per_socket)) + '        # How many tasks on each CPU or socket (not sure what this really means)' + '\n',
                         '#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically on nodes and sockets' + '\n',
                         '# Memory usage (MB)' + '\n', '#SBATCH --mem-per-cpu=' + str(self.memory_limit) + '\n',
                         '#SBATCH --output=' + os.path.join(self.folder, 'sbatch.out') + '   # Standard output and error log' + '\n']


        # The commands to execute will look the same regardless of whether you're using SLURM or not
        commands_to_execute = [self.preprocess_command + '\n \n',
                         'cd ' + self.folder + '; \n \n']

        if self.custom_command:
            commands_to_execute.append(self.custom_command)

        else:
            print('Create batch command, exec file = %s' % self.exec_file)
            run_str = mpi_str + ' ' + self.exec_file + ' ' + os.path.join(self.folder, inputs_file_name)

            # If we may be running on legacy, need to possibly use a different executable
            if 'legacy' in self.partitions:

                if '.GYRE.ex' not in run_str:
                    legacy_run_str = run_str.replace('.ex', '.GYRE.ex')
                else:
                    legacy_run_str = run_str

                # exec_dir_parts = self.execFile.split('/')
                print('Exec file: %s ' % self.exec_file)
                exec_dir = os.path.dirname(self.exec_file)
                print('Exec dir: %s' % exec_dir)

                custom_cmd = 'hs=$HOSTNAME \n' \
                             'if [[ $hs == *"gyre"* ]]; then \n' \
                             ' buildmush; \n cd %s ; \n' \
                             ' %s \n' \
                             'else\n' \
                             ' %s\n' \
                             'fi\n' % (self.folder, legacy_run_str, run_str)

                commands_to_execute.append(custom_cmd)

            else:
                commands_to_execute.append(run_str)

        commands_to_execute.append('\n' + self.postprocess_command + '\n')

        file_contents = ['#!/bin/bash \n']
        file_contents.extend(slurm_header)
        file_contents.append(dependency_str)
        file_contents.extend(commands_to_execute)


        fh.writelines(file_contents)

        fh.close()

    def run_task(self, runFileName='run.sh'):
        '''
        This method will write out a batch file for the slurm queuing system, then run it.
        If you don't have slurm, you'll need to rewrite this method
        '''
        self.write_batch_file(runFileName)
        slurm_command = 'sbatch' # change this if you're not using slurm!

        cmd = 'cd ' + self.folder + '; ' + slurm_command + ' ' + self.get_run_file(runFileName)

        result = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)

        print(result)

        job_id_search = re.findall('Submitted batch job (\d+)', result)

        if job_id_search:
            self.job_id = int(job_id_search[0])

        # Make a file in the output folder stating the slurm job id
        F = open(os.path.join(self.folder, 'jobid'), 'w')
        F.write(str(self.job_id))
        F.close()

        # Pause briefly in case submitting lots of slurm jobs
        time.sleep(0.5)

    def get_run_file(self, runFileName):
        return os.path.join(self.folder, runFileName)
