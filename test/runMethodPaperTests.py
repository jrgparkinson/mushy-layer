# Run all convergence tests for the methods paper

from runAMRConvergenceTest import runTest
import os

from colorama import Fore, Style

base_output_dir = '/home/parkinsonjl/mushy-layer/test/output/'
base_output_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/Test/'


if os.path.exists(base_output_dir):
	# check empty
	if os.listdir(base_output_dir) :
		print(Fore.RED + 'Warning - output directory not empty \n')
		print(Style.RESET_ALL)
				
else:
    os.makedirs(base_output_dir)

######################################
# 1) Diffusive solidification problem
#######################################
physicalProblem = 'noFlow'
AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform']}, 
{'max_level': 1, 'ref_rat': 2, 'run_types': ['amr']},
{'max_level': 2, 'ref_rat': 2, 'run_types': ['amr']}]
# While testing:
#AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform']}];

Nzs 	  = [16, 32, 64]
num_procs = [1] * len(Nzs)

# Setup up the post processing command
dataFolder = os.path.join(base_output_dir, 'AMRConvergenceTestNoFlow')
figureName = os.path.join(dataFolder, 'noFlow.pdf')
parDir = os.path.abspath(os.pardir)
matlabFolder = os.path.join(parDir, 'matlab', 'MethodsPaper')
analysis_command = 'cd ' + matlabFolder + '; \n \n matlab -nodisplay -nosplash -nodesktop -r "noFlowSolution(\'' + dataFolder + '\', \'' +figureName+ '\'); exit;"' 

# Run
extra_params = {'main.debug':'true'}
runTest(dataFolder, physicalProblem, AMRSetup, Nzs, num_procs, analysis_command)

######################################
# 2) Convection in a fixed porous medium
######################################
physicalProblem = 'convectionDB'
AMRSetup = [{'max_level': 0, 'ref_rat': 2, 'run_types': ['uniform']}, 
{'max_level': 1, 'ref_rat': 2, 'run_types': ['variable']}]

Nzs 	  = [128]
num_procs = [4]

chi  = 0.4
base_dataFolder = os.path.join(base_output_dir, 'AMRConvergenceTestConvectionDB')

Da_Ra_vals = [{'Da': 1e-6, 'RaT': [1e7, 1e8, 1e9]},
{'Da': 1e-2, 'RaT': [1e3, 1e4, 1e5, 5e5]}]

for Da_Ra in Da_Ra_vals:
	Da = Da_Ra['Da']

	for Ra in Ra_Ra['RaT']:

		extra_params['parameters.rayleighTemp'] = Ra
	 
		extra_params['parameters.darcy'] = Da*pow(1-chi,2)/pow(chi,3.0)
	                
		extra_params['bc.porosityHiVal']= str(chi) + ' ' + str(chi) # 0.4 0.4
		extra_params['bc.porosityLoVal']= str(chi) + ' ' + str(chi)

		output_dir = 'chi'+str(chi)+'-Da' + str(Da) + '-Ra' +  params['parameters.rayleighTemp']

		#extra_params = {}
		runTest(dataFolder, physicalProblem, AMRSetup, Nzs, num_procs, analysis_command, extra_params)


# 3) Convection in a fixed porous medium with variable porosity

# 4) Fully coupled porous hole

# 5) Fixed chill Hele-Shaw cell


