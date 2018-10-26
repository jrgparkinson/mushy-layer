import os
import sys
parDir = os.path.abspath(os.pardir)
pythonDir = os.path.join(parDir, 'python')
sys.path.append(pythonDir)

from SlurmTask import SlurmTask
from MachineSpecific import MachineSpecific
from mushyLayerRunUtils import readInputs, writeParamsFile

# May need to 
# chmod a+x ../execSubcycle/mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex
# if valgrind doesn't have permission to run

machineSpecific = MachineSpecific()
baseDir = os.path.join(machineSpecific.get_home_dir(), 'convection-in-sea-ice', 'MushyLayer')
memleakFolder = os.path.join(baseDir, 'memleakTest')

print('Home directory: %s' % machineSpecific.get_home_dir())
print('Base directory: %s' % baseDir)
print('Memleak folder: %s' % memleakFolder)


i = 0
outputFolder = os.path.join(memleakFolder, 'test-' + str(i))
while os.path.exists(outputFolder):
	i = i + 1
	outputFolder = os.path.join(memleakFolder, 'test-' + str(i))
	
print('Output folder for memleak test: %s \n' % outputFolder)
os.makedirs(outputFolder)     

inputsFile = outputFolder + '/inputs'
params = readInputs(memleakFolder + '/inputs')
params['main.output_folder'] = outputFolder
writeParamsFile(inputsFile, params)


program = os.path.join(baseDir, 'execSubcycle' , machineSpecific.get_program_name())

os.system('chmod a+x ' + program)

cmd = 'valgrind --leak-check=yes ' + program


num_proc = 1
timeLimit = 1.0
memoryLimit = 5000
numNodes = 1
mpirun = False

s = SlurmTask(outputFolder, 'memleakTest', cmd,  num_proc, timeLimit , memoryLimit, numNodes, mpirun)
s.writeSlurmFile()
s.runTask()