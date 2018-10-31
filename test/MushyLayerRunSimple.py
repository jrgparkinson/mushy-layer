import sys
import os
import re
import math
from os import listdir, rename
from os.path import isfile, join, exists
import subprocess
import socket

from mushyLayerRunUtils import readInputs, writeParamsFile

# This class basically just makes the output directory and dumps an input file into it
# It also performs various sanity checks
class MushyLayerRunSimple:

    def __init__(self, base_output_dir, num_proc, parameters, slurmJob, allowMultipleDirs, programName):
        
        self.num_proc = num_proc
        self.parameters = parameters
        self.program_name = programName
            
        # Allow us to create output dirs -0, -1, -2 etc.
        # If set to false, and the output dir ending -0 already exists,
        # then we won't do the run
        self.allowMultipleOutputDirs = allowMultipleDirs

        self.makeSlurmJob = True
        self.slurmJob = slurmJob

        self.base_output_dir = base_output_dir

        mushyLayerDir = os.path.dirname(os.path.dirname(__file__))
        self.exec_dir = os.path.join(mushyLayerDir, 'execSubcycle')

        #self.base_command = "mpirun --mca mpi_warn_on_fork 0  -np " + \
        #    str(self.num_proc) + " " + self.exec_dir + \
        #    self.program_name + " "

        # Make outputDir if it doesn't exist
        
        if not (exists(self.base_output_dir)):
            os.makedirs(self.base_output_dir)

        # Choose max grid size sensibly
        max_level = 0
        if 'main.max_level' in parameters:
                max_level = int(parameters['main.max_level'])

        if 'main.gridfile' in parameters and max_level > 0:
            
            with open(parameters['main.gridfile'], 'r') as f:
                #all_lines = f.read()
                maxBoxLength = 0.0
                pattern = '\(\((\d+),(\d+)\) \((\d+),(\d+)\) \((\d+),(\d+)\)\)'
                
                for line in f.readlines():
                    m = re.search(pattern, line)
                    if m:
                        
                        # Determine size of box
                        x_lo = int(m.group(1))
                        y_lo = int(m.group(2))
                        x_hi = int(m.group(3))
                        y_hi = int(m.group(4))
                        z_lo = int(m.group(5))
                        z_hi = int(m.group(6))

                        maxBoxLength = max(maxBoxLength, y_hi+1-y_lo)
                        maxBoxLength = max(maxBoxLength, x_hi+1-x_lo)
                        maxBoxLength = max(maxBoxLength, z_hi+1-z_lo)

                # Make sure max_grid_size is at least as big as the largest
                # box in the specified grid
                parameters['main.max_grid_size'] = str(max(int(parameters['main.max_grid_size']), maxBoxLength))
        elif 'main.num_cells' in parameters and 'main.max_grid_size' not in parameters:
                # want to choose max_grid_size to do parallel stuff most effectively
                # Want to choose a box size to make the parallel computations most
                # efficiently
                ncells = parameters['main.num_cells'].split(' ')
                grid_res_x = int(ncells[0])
                grid_res_y = int(ncells[0])
                if math.log(grid_res_x, 2.0) == math.floor(math.log(grid_res_x, 2.0)):
                        min_box_size = math.sqrt(grid_res_x * grid_res_x / num_proc)

                        # Round this up to the nearest power of 2
                        min_box_size = int(2**(math.ceil(math.log(min_box_size, 2))))
                        parameters['main.max_grid_size'] = str(min_box_size)
                
                

    def single_run(self, run_name):

        test = self.make_next_dir(join(self.base_output_dir, run_name))
        if test == -1:
            print('Run already done, skipping \n \n')
            return -1

        this_run_directory = self.get_most_recent_directory(run_name)

        exitStatus = self.run_model(self.parameters['main.num_cells'],
                                    join(self.base_output_dir, this_run_directory),
                                    False)

        return exitStatus

    def run_model(self,  grid_res,  output_folder,  analytic):

        # Perform various checks on the inputs

        # Remove trailing slash from output_folder
        if output_folder[-1] == "/":
            output_folder = output_folder[:-1]

        self.parameters['main.output_folder'] = output_folder

        # If this isn't specified, don't output plots for each iteration
        # (in case an input file has this option left in)
        # don't want to generate loads of output unnecessarily 
        if 'main.iteration_plot_interval' not in self.parameters:
            self.parameters['main.iteration_plot_interval'] = '-1'

        # Unless a max step has been explicitly specified, set it to something
        # big so we hit a steady state
        if 'main.max_step' not in self.parameters:
            self.parameters['main.max_step'] = '999999'

        if 'main.max_time' not in self.parameters:
            self.parameters['main.maxTime'] = '9999'


        # Finished checking inputs

        # Now need to pass all our parameters to the code. The easy way is to
        # do it with arguments, the hard way is to create a new inputs file
        # and edit it. The latter means we have a record of all parameters
        # used, which is preferable.

        # This is the easy way (just append this to the command)
        # parameters_string = ' '.join([key+'='+value for key, value
        #                             in self.parameters.items()])

        # Hard way
        # 1) get base inputs file
        new_input_file = output_folder + '/inputs'

        # copyfile(inputs_file, new_input_file)
        # Read in the file
        PARAMS_KEY = 'main.params_file'

        #inputsParams = readInputs(inputs_file)
        inputsParams = {}
        
        extraParams = {}
        if PARAMS_KEY in self.parameters:
            extraParams = readInputs(self.parameters[PARAMS_KEY])
        elif PARAMS_KEY in inputsParams:
            extraParams = readInputs(inputsParams[PARAMS_KEY])

        # Merge all dictionary's
        # Important that it's done in this order, so that parameters
        # are overwritten in the correct order. Priority is:
        # 1) params specified in this python scrupt
        # 2) params in main.params_file 
        # 3) params in the inputs file

        for key in extraParams.keys():
            if key not in self.parameters:
                self.parameters[key] = extraParams[key]

        for key in inputsParams.keys():
            if key not in self.parameters:
                self.parameters[key] = inputsParams[key]


        # If we have Kozeny-Carman permeability, make sure we have hele-shaw cell to limit the max permeability
        if self.parameters['parameters.permeabilityFunction'] == 2:
            self.parameters['parameters.heleShaw'] = 'true'
            
            
        # Let's add another item to the parameters: the mercurial repository version. Get this from BASH
        pipe = subprocess.Popen(
        ["hg", "identify", "--num"],
        stdout=subprocess.PIPE
        )
        repoVersion = pipe.stdout.read()
        repoVersion = repoVersion.replace('\n', '') # clean up a bit
        self.parameters['mercurialRepository'] = repoVersion
        
        # Also add the current machine to the parameter list so we now where each run was done from
        self.parameters['machineCodeIsRunOn'] = socket.gethostname()

        # Write out final params file
        writeParamsFile(new_input_file, self.parameters)

        # Set remaining parameters on the slurm jobs

        self.slurmJob.folder = output_folder
        self.slurmJob.execFile = self.exec_dir + self.program_name

        self.slurmJob.writeSlurmFile()
        self.slurmJob.runTask()

        return 1

    # Take some basename, and make the directory basename-(n+1) where
    # directories basename-0 to basename-n already exist
    def make_next_dir(self,  basename):
        i = 0
        new_dir = basename + "-" + str(i) + "/"
        while(os.path.isdir(new_dir)):
            print('    Path exists: ' + new_dir + ', ' + str(os.path.isdir(new_dir)))
            onlyfiles = [f for f in listdir(new_dir) if isfile(join(new_dir, f))]
            print(str(onlyfiles))
            
            i = i + 1
            new_dir = basename + "-" + str(i) + "/"
            
        # If we already have directory -0, we may not want to create -1, -2 etc.
        if i > 0 and not self.allowMultipleOutputDirs:
            return -1
            
        # Also check for directories ending in -steady with the same name
        if not self.allowMultipleOutputDirs:
            steadyName = new_dir.replace('-' + str(i), '-steady')
            if os.path.isdir(steadyName):
                print('    Path exists: ' + steadyName)
                return -1

        os.makedirs(new_dir)

        #print("New output directory is " + new_dir)
        return new_dir

    def get_most_recent_directory(self, run_name):
        basename = join(self.base_output_dir, str(run_name))
        i = 0
        new_dir = basename + "-" + str(i) + "/"
        while(os.path.exists(new_dir)):
            i = i + 1
            new_dir = basename + "-" + str(i) + "/"

        most_recent_dir = basename + "-" + str(i - 1) + "/"

        return most_recent_dir



