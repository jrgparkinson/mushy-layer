# Run all convergence tests for the methods paper

from runAMRConvergenceTest import runTest
import os

from colorama import Fore

output_dir = '/home/parkinsonjl/mushy-layer/test/output/'
output_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/Test/'


if os.path.exists(output_dir):
	# check empty
	if os.listdir(output_dir) :
		print(Fore.RED + 'Warning - output directory not empty')
		print(Style.RESET_ALL)
				
else:
    os.makedirs(output_dir)

# 1) Diffusive solidification problem
physicalProblem = 'noFlow'
AMRSetup = [{'max_level': 0, 'ref_rat': 0},
	{'max_level': 1, 'ref_rat': 2},
    {'max_level': 1, 'ref_rat': 4},
    {'max_level': 2, 'ref_rat': 2}]

Nzs 	  = [16, 32, 64]
num_procs = [1] * len(Nzs)

dataFolder = os.path.join(output_dir, 'NoFlow')
figureName = os.path.join(output_dir, 'NoFlow', 'noFlow.pdf')

parDir = os.path.abspath(os.pardir)
matlabFolder = os.path.join(parDir, 'matlab', 'MethodsPaper')
analysis_command = 'cd ' + parDir + '; \n \n matlab -nodisplay -nosplash -nodesktop -r "noFlowSolution(\'' + dataFolder + '\', \'' +figureName+ '\'); exit;"' 

runTest(output_dir, physicalProblem, AMRSetup, Nzs, num_procs, analysis_command)

# 2) Convection in a fixed porous medium
physicalProblem = 'convectionDB'
#runTest(output_dir, physicalProblem, AMRSetup, Nzs, num_procs)


# 3) Convection in a fixed porous medium with variable porosity

# 4) Fully coupled porous hole

# 5) Fixed chill Hele-Shaw cell


