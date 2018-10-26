# Compare a series of steady state solutions with different numbers of
# grid points and capture the differences in the primary variables (theta, solid fraction, ThetaL)
from subprocess import Popen, PIPE
import numpy as np
import matplotlib.pyplot as plt

import socket
hostname = socket.gethostname()

#Get the home directory, which should contain /Chombo and /convection-in-sea-ice
if (hostname == 'atmlxlap005'):
	get_home_dir = '/home/parkinsonjl'
elif (hostname == 'gyre4'):
	get_home_dir = '/local/home/gyre4/ice/Users/some2962'
elif(hostname== 'jamie-Ubuntu'):
    get_home_dir =  "/home/jamie"
else:
	get_home_dir = '~' #pray this works

#Machine specific vars
parallelProgram = "/compare2d.Linux.64.mpiCC.gfortran.DEBUG.MPI.ex"
serialProgram = "/compare2d.Linux.64.g++.gfortran.DEBUG.ex"

#Constants - don't change! (unless I change some of the mushyLayer code)
ANALYTIC = "analytic"
STEADY = "steady"
L1ERR = "L1"
L2ERR = "L2"
MAX = "Max"
SUM = "Sum"
CONCENTRATION = 4
THETA = 5
POROSITY = 8
LIQUID_CONC = 10
ChomboCompareDir = get_home_dir + "/Chombo/lib/util/ChomboCompare"
testDir = '/convection-in-sea-ice/test'

'''
##Vars - feel free to change
#############################################
#Variable to compare
VARIABLE = THETA

#Type of run to compare from
compareFrom = STEADY

#Type of run to compare to (the supposedly exact solution)
compareTo = STEADY

#Choose to compare either L1 or L2 error
errType = L1ERR

#Do a log-log plotError?
plotLog = True

#Are we working in parallel?
parallel = True
#############################################
'''

#Switching to a function interface
'''
Inputs:
VARIABLE - a number giving the line in the chomboCompare outputOneLevel file of the variable
             we're interested in. THETA and POROSITY are already defined
compareFromPts_arr - array of the different resolutions we'll be comparing
compareFrom - STEADY or ANALYTIC. Type of run to compare from (the lower resolution file)
compareTo - STEADY or ANALYTIC. Type of run to compare to (usually the same as compareFrom)
errType - L1ERR, L2ERR, MAXERR, SUMERR. Type of error to study convergence of.
plotLog - whether or not to make a log-log plotError. Usually true.
parallel - whether or not chomboCompare is being run in parallel. Usually true.
data_dir - where the data files to analyse are located

Outputs:
1. error files for each comparison made, of the format
	err[low_res_points][low_res_type - steady/analytic]-[high_res_points][high_res_type - steady/analytic]
2. Produces a plotError of error convergence
'''


def convergenceStudy(VARIABLE, compareFromPts_arr, compareToPts_arr, compareFrom, compareTo, errType, plotLog, parallel, data_dir):
	# if (len(compareFromPts_arr) != len(compareToPts_arr)):
   	#  print("Warning: compareFromPts and compareToPts are different lengths!!")
   	 
	compareFromPts_log_arr = []
	L1 = []
	L2 = []
	Max = []
	Sum = []

	pts_i=0

	while pts_i < len(compareFromPts_arr):
		compareFromPts = compareFromPts_arr[pts_i]
		filename = get_home_dir + data_dir + "bm1-" + str(compareFromPts) + "pts-" + compareFrom + ".2d.hdf5"

		compareToPts = compareToPts_arr[pts_i]

		if (compareTo == ANALYTIC):
			exactFile = get_home_dir + data_dir + "bm1-"+str(compareToPts)+"pts-analytic.2d.hdf5"
		elif(compareTo == STEADY):
			exactFile = get_home_dir + data_dir + "bm1-"+str(compareToPts)+"pts-steady.2d.hdf5"

		print("File to compare against: " + exactFile)

		if (filename == exactFile):
			print("Skipping " + filename)
			continue


		compareFromPts_log_arr.append(np.log10(compareFromPts))
	

		print("Comparing " + filename)

		errFile = "err" + str(compareFromPts) + compareFrom + "-" + str(compareToPts) + compareTo + ".hdf5";
	
		if (parallel):
			program = ChomboCompareDir + parallelProgram
		else:
			program = ChomboCompareDir + serialProgram

        
		process = Popen(program + " " + get_home_dir + testDir + "/bm1inputs.compare" + 
	" compare.computedRoot=" + filename + 
	" compare.exactRoot=" + exactFile + 
	" compare.errorRoot=" + errFile, shell=True)
		(outputOneLevel, err) = process.communicate()
		exit_code = process.wait()

		#print("exit code = " + str(exit_code))

		#in parallel, the outputOneLevel is actually put in a file - pout.0
		if (parallel):
			lines = open('pout.0','r')
		else:
			lines = outputOneLevel.split("\n")

		#print(lines)

		#now deal with the outputOneLevel
		L1.append([])
		L2.append([])
		Max.append([])
		Sum.append([])	
		i = 0
		for line in lines:
			#print (str(i) + ": " + line)
			#print (line)

			# Don't parse lines that don't contain errors
			if (i < 5 or i > 27):
				i = i+1
				continue

			parts = line.split(":")
			varName = str(parts[0])
			err = parts[1].split(",")
			err_L1 = float(err[0])
			err_L2 = float(err[1])
			err_Max = float(err[2])
			err_Sum = float(err[3])

			L1[-1].append(err_L1)
			L2[-1].append(err_L2)
			Max[-1].append(err_Max)
			Sum[-1].append(err_Sum)
			
			i=i+1

		#lines.close()
		#print("===========================================")
		pts_i = pts_i + 1

	#print (L1)
	#Now we have the data, do some analysis

	#Let's look at convergence for the variable and error type specified
	var = VARIABLE
	theta_L1_err = []
	theta_L2_log_err = []
	log_err = []
	secondOrderCompL1 = []
	firstOrderCompL1 = []
	secondOrderCompL2 = []
	firstOrderCompL2 = []
	halfOrderCompL2 = []
	i = 0
	
	if errType == L1ERR:
		allErrors = L1
	elif errType == L2ERR:
		allErrors = L2
	elif errType == MAX:
		allErrors = Max
	elif errType == SUM:
		allErrors = Sum
		
	errArray = []
	log_errArray = []
	
	#Go through the array of this error type for all vars, and extract the var we're interested in
	for err in allErrors:
		print(err)
		thisVarErr = abs(err[var])
		errArray.append(thisVarErr)
		log_errArray.append(np.log10(thisVarErr))
		firstOrderCompL1.append(log_errArray[0] - (compareFromPts_log_arr[i] - compareFromPts_log_arr[0]))
		secondOrderCompL1.append(log_errArray[0] - 2*(compareFromPts_log_arr[i] - compareFromPts_log_arr[0]))
		i=i+1



	plt.plotError(compareFromPts_log_arr, log_errArray, label="convergence", marker='x')
	plt.plotError(compareFromPts_log_arr, firstOrderCompL1, label="1st order", marker='x')
	plt.plotError(compareFromPts_log_arr, secondOrderCompL1, label="2nd order", marker='x')
	plt.legend(loc = 3)
	plt.title(errType + " Error convergence of " + compareFrom + " solution to " + compareTo + " solution")
	plt.xlabel('log10(grid points)')
	plt.ylabel('log10(' + errType + ' error)')
	plt.show()


	#print(theta_L1_err)
