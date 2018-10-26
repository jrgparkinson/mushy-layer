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
dataDir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/'
brinicleFolder = os.path.join(dataDir, 'brinicle')

print('Home directory: %s' % machineSpecific.get_home_dir())
print('Base directory: %s' % baseDir)
print('Brinicle folder: %s' % brinicleFolder)


i = 0
outputFolder = os.path.join(brinicleFolder, 'run-' + str(i))
while os.path.exists(outputFolder):
	i = i + 1
	outputFolder = os.path.join(brinicleFolder, 'run-' + str(i))
	

print('Output folder for brinicle simulation: %s \n' % outputFolder)
os.makedirs(outputFolder)     

inputsFile = outputFolder + '/inputs'
params = readInputs(brinicleFolder + '/inputs')
params['main.output_folder'] = outputFolder

# Other changes to params go here
params['parameters.compositionRatio']='1.01'
params['parameters.darcy']='0.01'
params['parameters.rayleighComp']='1e5'
params['parameters.rayleighTemp']='1e3'
params['parameters.prandtl']='1.0'
params['parameters.stefan']='5'

params['parameters.plumeBounds']= '0.35 0.65' #'0.45 0.55'
params['parameters.enthalpyPlume'] = '5.1'
params['parameters.bulkConcPlume'] = '-0.01'
params['parameters.inflowVelocity'] = '-1.0'

params['main.nondimensionalisation']='1'
params['main.num_cells']='64 64 8'
params['main.initial_cfl']='0.0001'
params['main.plot_interval'] = '50'

params['main.cfl']='0.01'

params['bc.velHi']='0 7' 
params['bc.velLo']='0 0'

writeParamsFile(inputsFile, params)


program = os.path.join(baseDir, 'execSubcycle' , machineSpecific.get_program_name())

#os.system('chmod a+x ' + program)
execFile = program


num_proc = 8
timeLimit = 7.0
memoryLimit = 5000
numNodes = 1
mpirun = True

s = SlurmTask(outputFolder, 'brinicle', execFile,  num_proc, timeLimit , memoryLimit, numNodes, mpirun)
s.writeSlurmFile()
s.runTask()