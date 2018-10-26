import sys
import os
pythonDir = os.path.abspath(os.pardir)
sys.path.append(pythonDir)

from os import listdir, rename

from os.path import isfile, join, exists
from shutil import copyfile
import subprocess
from subprocess import Popen
from MachineSpecific import MachineSpecific
import re
from util.timeout import timeout
#import numpy
import math
from util.mushyLayerRunUtils import readInputs, writeParamsFile
from SlurmTask import SlurmTask

# For running processes in parallel
from multiprocessing import Process, Queue

# For checking disk space every 60 seconds
import time

# For sending emails
import smtplib
from email.mime.text import MIMEText

# For working out what machine we're on
import socket


class MushyLayerRun:
    PLOT_SUFFIX = '.2d.hdf5'
    NARROW = 'narrow'
    SQUARE = 'square'
    WIDE = 'wide'
    CUSTOM = 'custom'

    FIXED_GRID = 'fixed'
    ADAPTIVE_GRID = 'adaptive'
    UNIFORM_GRID = ''  # replace single level with uniform grid

    def __init__(self, base_output_dir, num_proc, parameters, subcycled=False, programName = ''):
        
        self.num_proc = num_proc
        self.parameters = parameters
        self.subcycled = subcycled
        self.prevFinalChkFile = ''
        self.finalFlux = float('nan')
        
        # Machine specific variables (where the executable is located etc.)
        self.machine_specific = MachineSpecific()
        
        if programName:
            self.program_name = programName
            
        else:
            self.program_name = self.machine_specific.get_program_name()
        
        # Allow us to create output dirs -0, -1, -2 etc.
        # If set to false, and the output dir ending -0 already exists,
        # then we won't do the run
        self.allowMultipleOutputDirs = True

        self.makeSlurmJob = False
        self.slurmJob = -1

        #self.base_output_dir = self.get_home_dir() + '/convection-in-sea-ice/MushyLayer/test'
        self.base_output_dir = base_output_dir
        if subcycled:
            self.exec_dir = self.get_home_dir() + \
                '/convection-in-sea-ice/MushyLayer/execSubcycle/'
            #self.parameters['main.use_subcycling'] = '1'
        else:
            self.exec_dir = self.get_home_dir() + \
                '/convection-in-sea-ice/MushyLayer/execNonSubcycle/'
            #self.parameters['main.use_subcycling'] = '0'

        self.base_command = "mpirun --mca mpi_warn_on_fork 0  -np " + \
            str(self.num_proc) + " " + self.exec_dir + \
            self.program_name + " "

        # Make outputDir if it doesn't exist
        
        if not (exists(self.base_output_dir)):
            os.makedirs(self.base_output_dir)

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
                
                
        #print('Max box size = ' + str(parameters['main.max_grid_size']) + '\n')

    def finalChkFile(self):
        return self.prevFinalChkFile

    def getFinalFlux(self):
        return self.finalFlux

    def grid(self, res):
        if self.grid_type == self.NARROW:
            height = res
            width = int(height / 4)
            # Can't have a width (or height) less than 4 (number of ghost cells
            # or some things won't work)
            width = max(width, 4)
            height = max(height, 4)

            return [width, height, width]

        elif self.grid_type == self.WIDE:
            height = res
            width = int(height * 4)
            return [width, height, width]

        # Square grid is default
        elif self.grid_type == self.CUSTOM:
            return self.custom_grid
        else:
            return [res, res, res]

    # def get_new_inputs_file(self, res):
    #    this_res_prefix = self.get_plot_prefix() + str(res) + "pts"
    #    return self.exec_dir + "inputs-" + this_res_prefix

    def run_model(self,  grid_res,  output_folder,  analytic):

        #if not "main.restart_file" in self.parameters:
        #    print("    n_cell for this run = " + self.parameters['main.num_cells'])

        #inputs_file = self.exec_dir + "inputs-standard"

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

        # We now have all parameters
        # Get an estimate for the max possible size of the run
        #ncells = self.parameters['main.num_cells'].split(' ')
        #numPlotFiles = int(float(self.parameters['main.max_step'])/float(self.parameters['main.plot_interval']))
        #if float(self.parameters['main.checkpoint_interval']) > 0.0:
        #    numChkFiles = int(float(self.parameters['main.max_step'])/float(self.parameters['main.checkpoint_interval']))
        #else:
        #    numChkFiles = 0
            
        #filesize = int(ncells[0])*int(ncells[1])*0.000366 # empirical number, in MB
        #totalsize = filesize*(numPlotFiles + numChkFiles) # plot file size is roughly same as chk file size
        #if totalsize > 10000:
        #    print("Warning, total size of run potentially ~" + str(totalsize/1000) + "GB" )

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


        # Write out final file
        # Sort params by key so there's a reasonably logical order
        param_keys = self.parameters.keys()
        param_keys = sorted(param_keys)

        # output_file = ''
        # for key in param_keys:
        #     if key == PARAMS_KEY:
        #         continue
        #     output_file = output_file + '\n' + \
        #         key + '=' + str(self.parameters[key])

        # with open(new_input_file, 'w') as f:
        #     f.write(output_file)
        writeParamsFile(new_input_file, self.parameters)
            
        if self.makeSlurmJob:
            # Set remaining parameters on the slurm jobs

            self.slurmJob.folder = output_folder
            self.slurmJob.execFile = self.exec_dir + self.program_name

            self.slurmJob.writeSlurmFile()
            self.slurmJob.runTask()

        

            return 1

        else:

            command = (
                'cd ' + output_folder + ' ; ' + self.base_command + ' ' + new_input_file)

            print(command + "\n")

            process = Popen(command, shell=True)
            exit_code = process.wait()

            
            self.prevFinalChkFile = output_folder + '/' + self.get_final_output(output_folder,  self.parameters['main.chk_prefix'])
            print('    Final checkpoint file: ' + self.prevFinalChkFile)
            print('    File name: ' + self.get_final_output(output_folder,  self.parameters['main.chk_prefix']))

            return exit_code

   

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

    def move_output_file(self,  filename,  this_res_directory):
        for n in range(0, self.num_proc):
            fname = filename + "." + str(n)
            if os.path.isfile(fname):
                rename(fname, this_res_directory + fname)

    # No longer have to do this - done at start of run
    # def copy_inputs(self,  this_res_directory, res):
    #    inputs_file_name = self.get_new_inputs_file(res)
    #    os.rename(inputs_file_name, this_res_directory + "inputs")

    def copy_final_output(self, grid_res, destination_dir=""):
        directory = self.get_most_recent_directory(grid_res)
        #prefix = self.get_res_prefix(grid_res)
        prefix = self.get_plot_prefix()

        if not os.path.exists(directory):
            return

        if not destination_dir:
            destination_dir = self.base_output_dir

        # If we get this far the directory must exist
        output_dir = join(self.base_output_dir, directory)
        old_name = self.get_final_output(output_dir, prefix)
        new_name = prefix + "-steady" + self.PLOT_SUFFIX

        old_location = join(directory, old_name)
        new_location = join(destination_dir, new_name)

        copy_cmd = 'cp ' + old_location + ' ' + new_location
        # print (copy_cmd)
        os.system(copy_cmd)

    def copy_final_output_dt(self, grid_res, dt, destination_dir=""):
        directory = self.get_most_recent_directory(dt)
        prefix = self.get_plot_prefix()
        new_prefix = self.get_dt_prefix(dt)

        #print(directory)
        #print(prefix)

        if not os.path.exists(directory):
            return

        if not destination_dir:
            destination_dir = self.base_output_dir

        # If we get this far the directory must exist
        output_dir = join(self.base_output_dir, directory)
        old_name = self.get_final_output(output_dir, prefix)
        new_name = new_prefix + "-steady" + self.PLOT_SUFFIX

        old_location = join(directory, old_name)
        new_location = join(destination_dir, new_name)

        copy_cmd = 'cp ' + old_location + ' ' + new_location
        #print (copy_cmd)
        os.system(copy_cmd)

    @timeout(5)
    def mkdirOnCloud(self, dir):
        os.system('ssh cloud mkdir -p ' + dir)

    # Copy across entire output directory
    def copy_all_output_to_cloud(self, res):
        directory = self.get_most_recent_directory(res)
        prefix = self.get_plot_prefix()

        cloud_base_dir = '/local/home/cloud/ice/Users/parkinson/'

        if not os.path.exists(directory):
            return

        # directory currently looks like bm2-Ra100/64-5
        # need to change the number on the end so we don't overwrite data on
        # cloud

        # first ensure the output folder exists on cloud
        self.mkdirOnCloud(cloud_base_dir + self.base_output_dir)

        # then get directory list from cloud
        dir_list_command = 'ls -d ' + cloud_base_dir + \
            self.base_output_dir + '/' + str(res) + '-*/'
        full_command = "ssh cloud " + dir_list_command
        print(full_command + '\n')
        ssh = subprocess.Popen(["ssh", "%s" % 'cloud', dir_list_command],
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        dirs = ssh.stdout.readlines()
        indexes = []

        # Take just the number after the resolution '128-i'
        for dir in dirs:
            parsed = dir.replace(cloud_base_dir + self.base_output_dir, '')
            parsed = parsed.replace('\n', '')
            parsed = parsed.replace('/' + str(res) + '-', '')
            parsed = parsed.replace('/', '')
            indexes.append(int(parsed))

        # Generate new directory name of the form 'res-i'
        indexes.sort()
        if indexes:
            new_i = indexes[-1] + 1
        else:
            new_i = 1

        new_full_directory = cloud_base_dir + self.base_output_dir + \
            '/' + str(res) + '-' + str(new_i) + '/'

        # scp_command = 'scp -r ' + directory + ' cloud:' + new_full_directory
        # print(scp_command)
        rsync_command = 'rsync -avz --remove-source-files -e ssh ' + \
            directory + ' cloud:' + new_full_directory

        exit_code = os.system(rsync_command)

        # Remove directory if files successfully copied across
        if exit_code == 0:
            # This will only reomve empty directories (which is good)
            os.rmdir(directory)

    # Requires absolute path to directory
    def get_final_output(self,  directory,  prefix):
        # Get all files in the directory with the given prefix
        allfiles = [
            f for f in listdir(directory) if isfile(join(directory, f))]
        #print(allfiles)
        #print(prefix)
        this_res_files = [f for f in allfiles if (prefix in f)]
        this_res_files.sort()

        # Find the final file, rename and copy it
        if (len(this_res_files) > 0):
            return this_res_files[-1]

        return ''

    def get_res_prefix(self, grid_res):
        return self.get_plot_prefix()
        #return self.get_plot_prefix() + str(grid_res) + "pts"

    def get_dt_prefix(self, dt):
        return self.get_plot_prefix() + str(dt) + "dt"

    def get_plot_prefix(self):
        return self.parameters['main.plot_prefix']
    
    def get_check_prefix(self):
        return self.parameters['main.chk_prefix']
    
    def diskSpaceCheck(self, q):
        starttime = time.time()
        sleepTime = 60.0  # seconds
        warningDiskSpace = 10000  # 10000 #MBs
        runFinished = q.get()

        while runFinished[0] is False:
            st = os.statvfs(self.machine_specific.local_dir())
            freeSpace = st.f_bavail * st.f_frsize / 1024 / 1024  # in MB

            #print('Free space remaining: ' + str(freeSpace) + 'MBs')

            # Send warning email if less than 10GB
            if (freeSpace < warningDiskSpace):
                msgText = "Warning - low disk space on " + self.machine_specific.get_machine_name(
                ) + ". Remaining space: " + str(freeSpace) + "MBs."
                print(msgText)
                msg = MIMEText(msgText)

                # me == the sender's email address
                # you == the recipient's email address
                email = 'james.parkinson@physics.ox.ac.uk'
                msg['Subject'] = 'Warning - low disk space on ' + \
                    self.machine_specific.get_machine_name()
                msg['From'] = email
                msg['To'] = email

                # Send the message via our own SMTP server, but don't include the
                # envelope header.
                #try:
                #    # send email
                #    s = smtplib.SMTP('localhost')
                #    s.sendmail(email, [email], msg.as_string())
                #    s.quit()
                #except socket.error:
                #    print("Couldn't send warning email")
               

            remainingSleepTime = sleepTime - ((time.time() - starttime) % 60.0)
            if remainingSleepTime < 0:
                remainingSleepTime = sleepTime
            # print("Remaining sleep time = " + str(remainingSleepTime) + "\n")

            time.sleep(remainingSleepTime)
            try:
                runFinished = q.get_nowait()
            except:
                runFinished = [False]
                pass

            # print("Run finished? " + str(runFinished) + ". Sleep time remaining = "+ str(remainingSleepTime))

    # Retained for backward compatibility
    def run(self, grid_resolutions, fixed_grid=False):

        exitStatuses = []

        for grid_res in grid_resolutions:

            self.parameters['main.num_cells'] = self.grid(grid_res)
            

            if fixed_grid:
                grids_file = (self.exec_dir + "grids-" + str(grid_res) + "-" +
                              str(2 * grid_res) + "-" + str(4 * grid_res))

                self.parameters['main.grids_file'] = grids_file
                self.parameters['main.max_box_size'] = '1024'
                self.parameters['main.regrid_interval'] = '-1'

            exitStatus = self.single_run(str(grid_res))
            # self.process_output(grid_res)

            exitStatuses.append(exitStatus)

        return exitStatuses

    def single_run(self, run_name):
        #q = Queue()
        #q.put([False])
        #p2 = Process(target=self.diskSpaceCheck, args=(q,))
        #p2.start()

        test = self.make_next_dir(join(self.base_output_dir, run_name))
        if test == -1:
            print('Run already done, skipping \n \n')
            return -1
            
        this_run_directory = self.get_most_recent_directory(run_name)

        exitStatus = self.run_model(self.parameters['main.num_cells'],
                                    join(self.base_output_dir, this_run_directory),
                                    False)

        # Still want to copy stuff even if program crashed
        #print("Execution finished, move files in directory: " +
        #      this_run_directory + " with prefix: " +
        #      self.get_plot_prefix() + "\n")

        #q.put([True])

        # Get the final flux
        #diagFile = os.path.join(this_run_directory, 'diagnostics.out')
        #print('Getting final flux from diagnostics file: ' + diagFile)
        #f = open(diagFile,'r') 
        #fContents = f.read() 
        #print(fContents)
        # \d+\.\d+
        #m = re.findall(r"Solute flux\: top \= ([^,]*),", fContents)
        #if m:
        #    self.finalFlux = float(m[-1])
        #else:
        #    print('Could not find final flux')

        #self.copy_final_output(run_name)

        return exitStatus

    def runForDt(self, grid_res, dt, fixed_grid=False):
        self.make_next_dir(join(self.base_output_dir, str(dt)))
        this_dt_directory = self.get_most_recent_directory(dt)
        this_dt_prefix = self.get_res_prefix(dt)

        self.parameters['main.fixed_dt'] = str(dt)

        if fixed_grid:
            grids_file = (self.exec_dir + "grids-" + str(grid_res) + "-" +
                          str(2 * grid_res) + "-" + str(4 * grid_res))

            self.parameters['main.grids_file'] = grids_file
            self.parameters['main.max_box_size'] = '1024'
            self.parameters['main.regrid_interval'] = '-1'

            # Ensure this is set to 2
            self.max_level = 2

        self.run_model(self.grid(grid_res),
                       join(self.base_output_dir, this_dt_directory), False)

        # Still want to copy stuff even if program crashed
        print("Execution finished, move files in directory: " +
              this_dt_directory + " with prefix: " + this_dt_prefix + "\n")
        # self.process_output_dt(grid_res, dt)

        self.copy_final_output_dt(grid_res, dt)

    # restartDir should look like 256-10
    # restartFile (if specified) should look like   bm1-chk-64pts0325.2d.hdf5
    def continue_run(self,  restart_dir, restart_file=""):
        print("continue previous model")
        full_restart_dir = join(self.base_output_dir, self.base_output_dir, restart_dir)

        if not restart_file:
            # get the latest chk file from this run
            restart_file = self.get_final_output(
                full_restart_dir,  self.get_check_prefix())
            if not restart_file:
                print(
                    "No checkpoint file found in directory " +
                    full_restart_dir)
                return

        print("Restarting from file: " + restart_file)

        restart_file_full_path = self.full_path(
            restart_dir + "/" + restart_file)
        # restart_command = " main.restart_file = " + restart_file_full_path
        self.parameters["main.restart_file"] = restart_file_full_path
        # print(restart_command)

        # get the grid resolution from the restart directory
        parts = restart_dir.split("-")
        height = int(parts[0])

        run_output_dir = self.make_next_dir(
            self.base_output_dir + "/" + restart_dir + "-cont")
        print("output_dir = " + run_output_dir)
        self.run_model(
            self.grid(height),  run_output_dir,  False)

        self.move_output_file('pout', run_output_dir)
        self.move_output_file('time.table',  run_output_dir)
        self.copy_inputs(run_output_dir, height)

    def calc_analytic_soln(self,  grid_resolutions):
        for grid_res in grid_resolutions:
            self.run_model(self.grid(grid_res),  join(
                self.base_output_dir, self.base_output_dir), True)

    # Append /..../test/output_dir/ to the start of the relative path
    def full_path(self,  relative_path):
        full_path = self.base_output_dir + "/" + self.base_output_dir + "/" + relative_path
        return full_path

    def get_home_dir(self):
        return self.machine_specific.get_home_dir()

