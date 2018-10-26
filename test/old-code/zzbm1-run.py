import socket
hostname = socket.gethostname()

#Get the home directory, which should contain /Chombo and /convection-in-sea-ice
if (hostname == 'atmlxlap005'):
	homeDir = '/home/parkinsonjl'
elif (hostname == 'gyre4'):
	homeDir = '/local/home/gyre4/ice/Users/some2962'
else:
	homeDir = '~' #pray this works
	
testDir = homeDir + '/convection-in-sea-ice/test'

from os import listdir, rename, chdir
from shutil import copyfile
from os.path import isfile, join
from subprocess import Popen, PIPE
import os

import sys
chomboCompareDir = homeDir + '/Chombo/lib/util/ChomboCompare'
sys.path.insert(0, chomboCompareDir)
import compare

'''
This script is designed to automate the following processes:

1) Run the sea ice model for the directional solidification problem at grid resolutions of 
	32, 64, 128, 256, 512, 1024
2) Rename the final (steady state) hdf5 files appropriately and move to /bm1/. Move the rest into the directory
	/bm1/[num-pts]pts/
along with the pout.n files
3) Calculate analytic solutions at each grid resolution and rename
4) Run chomboCompare and produce convergence plots for verification by eye that the code is working correctly

We assume that we are running in parallel, so output is sent to pout.[n] and the program is run with mpirun

'''

#### Vars to change
grid_resolutions = [256]
#grid_resolutions=[1024]
#grid_resolutions=[2048]
restartDirs = ['256-13/bm1-chk-256pts0013.2d',  '64-4', '128-3',  '256-3',  '512-2'] #can start each run from a checkpoint file if we want
num_proc = 2
execDir = homeDir + "/convection-in-sea-ice/MushyLayer/execNonSubcycle/"
base_command = "mpirun -np "+str(num_proc)+ " " + execDir + "mushyLayer2d.Linux.64.mpiCC.gfortran.DEBUG.MPI.ex"
plot_prefix = "bm1-"
check_prefix = "bm1-chk-"
plot_suffix = ".2d.hdf5"
outputDir = "bm1/" #this is relative to current directory
maxLevel = 0 # How many levels in the AMR hierarchy? (0=one level, 1=two levels etc.)
run_model = True
calc_analytic = True

plot_interval = 1 #20

run_from_restart = False
run_chomboCompare = False



### Don't change anything else ###########

# Loop over the different grid sizes
grid_i = 0
for grid_res in grid_resolutions:

	### 1) Run sea ice model 
	this_res_prefix = plot_prefix + str(grid_res) + "pts"
	this_res_check_prefix = check_prefix + str(grid_res) + "pts"
	width = grid_res / 8
	n_cell = str(width) + " " + str(grid_res) + " " + str(width)

	print("n_cell for this run = " + n_cell)

	if run_model:
		
		restartCommand = ""
		if run_from_restart:
			restartDir = restartDirs[grid_i]
			restartFile = testDir + "/" + outputDir + restartDir
			print("Restarting from file: " + restartFile)
			restartCommand = " main.restart_file = " + restartFile
		

        
		i = 0
        #Make sure the destination directory exists before we start copying
        #i.e. ensure we have bm1/64-i/ where i increases for each run number
		thisResDirectory = outputDir + str(grid_res) + "-" + str(i) + "/"
		while(os.path.exists(thisResDirectory)):
			i=i+1
			thisResDirectory = outputDir + str(grid_res) + "-" + str(i) + "/"
			
		os.makedirs(thisResDirectory)
		print("New output directory is " + thisResDirectory)
		
        
		print ("Running model with command: ")
		
		print(base_command)

         
		process = Popen(base_command + 
" ../MushyLayer/execNonSubcycle/inputs" + 
" main.num_cells=" + n_cell +
" main.plot_interval = " + str(plot_interval)+ 
" main.plot_prefix=" + this_res_prefix + 
" main.check_prefix=" + this_res_check_prefix + 
" main.maxLevel = " + str(maxLevel) + 
" main.output_folder = " + thisResDirectory +
restartCommand + 
" main.maxStep=5000", #make sure we're running to steady state
shell=True)
		(output, err) = process.communicate()
		exit_code = process.wait()
		print("Finished simulation, moving files")
	else:
		print("Not running model")
	
	### 2) Rename files and copy output to a new directory
	#First find the most recent directory of runs
	i = 0
	#Make sure the destination directory exists before we start copying
	#i.e. ensure we have bm1/64-i/ where i increases for each run number
	thisResDirectory = outputDir + str(grid_res) + "-" + str(i) + "/"
	while(os.path.exists(thisResDirectory)):
		i=i+1
		thisResDirectory = outputDir + str(grid_res) + "-" + str(i) + "/"
		
	thisResDirectory = outputDir + str(grid_res) + "-" + str(i-1) + "/"
	
	if(os.path.exists(thisResDirectory)):
		allfiles = [f for f in listdir(thisResDirectory) if isfile(join(thisResDirectory, f))]
		this_res_files = [f for f in allfiles if (this_res_prefix in f)]
		this_res_files.sort()

		# Find the final output, rename and copy it
		if (len(this_res_files) > 0 and run_model):
			old_name = this_res_files[-1]
			new_name = this_res_prefix + "-steady" + plot_suffix
			#copyfile(thisResDirectory + old_name, outputDir + new_name)
			os.system('cp '+ thisResDirectory+old_name + ' ' + outputDir + new_name)


	#Copy across the relevant pout files
	if run_model:
		#Copy across the inputs file
		copyfile(execDir + "inputs", thisResDirectory + "inputs")

		for n in range(0, num_proc):
			fname = "pout." + str(n)
			if os.path.isfile(fname):
				rename(fname, thisResDirectory + fname)
			
			fname = "time.table." + str(n)
			if os.path.isfile(fname):
				rename(fname, thisResDirectory + fname)

	
	print("Done moving files")
	
	### 3) calculate analytic solution
	if calc_analytic:
		print("Calculating analytic solution, running ")
		print(base_command)
	
		
		process = Popen(base_command + 
" ../MushyLayer/execNonSubcycle/inputs" + 
" main.num_cells=" + n_cell +
" main.plot_prefix=" + this_res_prefix + 
" main.check_prefix=" + this_res_check_prefix + 
" main.maxLevel = " + str(maxLevel) + 
" main.output_dir = " + outputDir +
" main.printAnalyticSoln=true", shell=True)
		(output, err) = process.communicate()
		exit_code = process.wait()
	
		#Get the analytic solution and move/rename it - it will be the plot file at time 0
		#analytic_file = outputDir + this_res_prefix + "-analytic" + plot_suffix
		#new_location = outputDir + this_res_prefix + "-analytic" + plot_suffix

#		print(analytic_file)
	#	rename(analytic_file, new_location)

		print("Calculated analytic solution, moving to next grid size")
	else:
		print("Moving to next grid size")

	print("======================================================")
	grid_i = grid_i + 1

print("All simulations complete!")
print("======================================================")

### 4) Run chomboCompare
if run_chomboCompare:
	print("Running chomboCompare")
	print("======================================================")
	chdir(chomboCompareDir)
	plotLog = True
	parallel = True
	compareToPts = [2048, 2048, 2048, 2048]


	compare.convergenceStudy(compare.THETA, grid_resolutions, compareToPts, compare.STEADY, compare.ANALYTIC, compare.L1ERR, plotLog, parallel, execDir + outputDir)

		


	


