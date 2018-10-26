from MushyLayerRun import MushyLayerRun
import os
#import numpy
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile, hasReachedSteadyState
from ChkFile import ChkFile
import sys, getopt
from os import listdir
from os.path import isfile, join

#thisDir = os.getcwd()
#sys.path.append('/path/to/application/app/folder')
currentDir = os.getcwd()
MushyLayerDir = currentDir.replace('python', '')

executeCommands = True
# This script will find all valid runs in the folder specified, re run with
# the additional options specified, and put the output in the new folder

# For local laptop:
originalBaseDir = '/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/optimalStates-highRes-newSoluteFlux/'
newBaseDir = '/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/optimalStates-highRes-WellsUnits/'
execLoc = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/mushyLayer2d.Linux.64.mpiCC.gfortran.DEBUG.MPI.ex'
#originalBaseDir = '/media/parkinsonjl/DATA/optimalStates-highRes/'
#newBaseDir = '/media/parkinsonjl/DATA/optimalStates-highRes-new/'


# For gyre:
originalBaseDir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/optimalStates-highRes-newSoluteFlux/'
newBaseDir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/optimalStates-highRes-WellsUnits/'
execLoc = MushyLayerDir + '/execSubcycle/mushyLayer2d.Linux.64.mpiCC.gfortran.OPTHIGH.MPI.ex'

if not os.path.isfile(execLoc):
	sys.exit("Error: executable not found")


# 4 processors should speed things up a little
num_proc = 1

foldersToProcess = []
# Get list of folders that end in -0

for f in next(os.walk(originalBaseDir))[1]:
	#print(f[-2:])
	
	# Only do folders which have reached steady state
	full_folder = os.path.join(originalBaseDir, f)
	steadyState = hasReachedSteadyState(full_folder)
	
	if f[-2:] == '-0' and steadyState:
		foldersToProcess.append(f)
		
		
# Just do a few as a test
foldersToProcess = foldersToProcess[:2]
#foldersToProcess = ['CR1.5RaC80Le200ChiCubedPermeabilitypts64-0']
print(foldersToProcess)

for folder in foldersToProcess:
	old_folder = os.path.join(originalBaseDir, folder)
	output_folder = os.path.join(newBaseDir, folder)
	
	if os.path.isdir(output_folder):
		continue
	
	if executeCommands:
		os.mkdir(output_folder)
	
	inputs = readInputs(os.path.join(originalBaseDir, folder, 'inputs'))
	inputsLoc = os.path.join(newBaseDir, folder, 'inputs')
	
	# Get the restart file
	files = [f for f in listdir(old_folder) if (isfile(join(old_folder, f)) and f[:3] == 'chk')]
	files = sorted(files)
	if len(files) > 0:
		restart_file = files[-1]
	else:
		print('no restart files!')
		print(files)
		continue
	
	print(restart_file)
	restart_file = os.path.join(old_folder, restart_file)
	# Modify inputs as required
	inputs['main.output_folder'] = output_folder
	inputs['main.restart_file'] = restart_file
	inputs['main.plot_interval'] = '100000'
	inputs['main.max_time'] = '1000'
	inputs['main.max_step']='99999999'
	inputs['main.steady_state']= '5e-3'
	
	# Wells units:
	inputs['parameters.referenceTemperature'] = '-3.5'
	inputs['parameters.referenceSalinity'] = '35'
	
	# Katz units:
	inputs['parameters.prevReferenceTemperature'] = '-23'
	inputs['parameters.prevReferenceSalinity'] = '230'
	
	# Convert boundary conditions:
	inputs['parameters.bottomEnthalpy'] = float(inputs['parameters.bottomEnthalpy']) - 1.0
	inputs['parameters.topEnthalpy'] = float(inputs['parameters.topEnthalpy']) - 1.0
	inputs['parameters.topBulkConc'] = float(inputs['parameters.topBulkConc']) + 1.0
	inputs['parameters.bottomBulkConc'] = float(inputs['parameters.bottomBulkConc']) + 1.0
	
	# Convert non dim vals:
	inputs['parameters.compositionRatio'] = float(inputs['parameters.compositionRatio']) - 1.0
	
	# For testing
	#inputs['main.plot_interval'] = '1'
	#inputs['
	
	if executeCommands:
		writeParamsFile(inputsLoc, inputs)
	
	cmd = 'cd ' + output_folder + '; mpirun -np ' + str(num_proc) + ' '  + execLoc + ' ' + inputsLoc 
	
	print(cmd)
	if executeCommands:
		# Execute command
		os.system(cmd)
