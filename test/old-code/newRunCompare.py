from convergenceStudy import convergenceStudy, plotError, plotErrorTime, plotTime
from modelOutput import plotNumIters
import ModelOutput.modelOutput
from MachineSpecific import  MachineSpecific
import numpy as np
import matplotlib.pyplot as plt

parallel = True


compareFromPts = [16, 32, 64]
compareToPts = [128]

machSpec = MachineSpecific()


compareFrom = convergenceStudy.CALCULATED
compareTo = convergenceStudy.ANALYTIC

compareTo = 'bm2-Ra100'
comapareFrom = 'bm2-fixed-Ra100'


cStudy1 = convergenceStudy(compareFromPts, compareToPts, compareFrom, compareTo, parallel,  '/home/parkinsonjl/Dropbox/DPhil/Data/')
#cStudy3 = convergenceStudy(compareFromPts, compareToPts, compareFrom, compareTo, parallel, '/convection-in-sea-ice/test/bm1-2level-2proc/')

cStudy2 = convergenceStudy(compareFromPts, compareToPts, compareFrom, compareTo, parallel, '/home/parkinsonjl/Dropbox/DPhil/Data/')
#cStudy4 =  convergenceStudy(compareFromPts, compareToPts, compareFrom, compareTo, parallel, '/convection-in-sea-ice/test/bm1-3level-16proc/')
#cStudy5 = convergenceStudy(compareFromPts, compareToPts, compareFrom, compareTo, parallel, '/convection-in-sea-ice/test/bm1-1level-16proc/')
# cStudy1.plotErrorResolution(convergenceStudy.THETA, convergenceStudy.L1ERR)
# cStudy1.plotErrorTime(convergenceStudy.THETA, convergenceStudy.L1ERR)
# cStudy1.plotTimeResolution()


cStudyArr = [cStudy1, cStudy2] #cStudy3, cStudy4, cStudy5, 
plotError(cStudyArr, convergenceStudy.THETA, convergenceStudy.L1ERR, False)
plotErrorTime(cStudyArr, convergenceStudy.THETA, convergenceStudy.L1ERR)
plotTime(cStudyArr)


# 
range = np.arange(5)
range = [0,10, 50, 100, 300]
#range = [10]
# 
# #Load just a single outputOneLevel file to analyse the iterations
outputOneLevel = ModelOutput(64, 1024, ModelOutput.CALCULATED, ModelOutput.ANALYTIC, machSpec.get_home_dir() + '/convection-in-sea-ice/test/bm1-1level-1proc/', parallel)
# outputOneLevel.plotConvergence(range, ModelOutput.CONC_RESIDUAL)

# outputOneLevel.plotNumIterations()
# outputOneLevel.plotConvergenceWithSign(range)
#  
# # 
# # 
outputAdaptive = ModelOutput(32, 1024, ModelOutput.CALCULATED, ModelOutput.ANALYTIC, machSpec.get_home_dir() + '/convection-in-sea-ice/test/bm1-2level-1proc/', parallel)
# outputAdaptive.plotConvergence(range, ModelOutput.CONC_RESIDUAL)

# outputAdaptive.plotNumIterations()
# outputAdaptive.plotConvergenceWithSign(range)


# plotNumIters([outputAdaptive, outputOneLevel])


# 
# outputAdaptive2 = ModelOutput(64, 1024, ModelOutput.CALCULATED, ModelOutput.ANALYTIC, machSpec.get_home_dir() + '/convection-in-sea-ice/test/bm1-2level-2proc/', parallel)
# outputAdaptive2.plotConvergence(range)
# # outputAdaptive.plotNumIterations()
# outputAdaptive2.plotConvergenceWithSign(range)
# 
# 




#Always have this!
plt.show()