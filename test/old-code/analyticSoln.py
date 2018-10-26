import numpy as np
import matplotlib.pyplot as plt

zVals = []
thetaVals = []
solidFracVals = []
porosityVals = []
enthalpyVals = []
v = 0.7
L = 1
stefan = 5.7

analyticSolnWorster = 'analyticSolnWorster.data'
analyticSolnKatz = 'analyticSolnKatz.data'
steadyStateSoln64 = 'pltFromZero64points.data'
steadyStateSoln128 = 'pltFromZero128points.data'
analyticSoln = 'analyticSoln.data'

filename = analyticSoln

with open(filename, 'r') as f:
	for line in f:
		data = line.split(',')
		z = float(data[0])
		theta = float(data[1])
		solidFrac = float(data[2])
		porosity = 1-solidFrac
		enthalpy = stefan*porosity + theta
	
		zVals.append(z)
		thetaVals.append(theta)
		solidFracVals.append(solidFrac)
		porosityVals.append(porosity)
		enthalpyVals.append(enthalpy)


dz = zVals[2]-zVals[1]
dz=dz*L
print("dz = " + str(dz))

enth = np.array(enthalpyVals, dtype=np.float)
gradEnth = np.gradient(enth, dz)
enthalpyAdvection = [v*x for x in gradEnth]

thet = np.array(thetaVals, dtype=np.float)
gradThet = np.gradient(thet, dz)
lapThet = np.gradient(gradThet, dz)

diff = enthalpyAdvection + lapThet

plotVars = False
if (plotVars):
	plt.plotError(zVals, thetaVals)
	plt.plotError(zVals, solidFracVals)
	plt.show()
	
plotGrad = True
if(plotGrad):
	plt.plotError(zVals,  gradThet,  label='dtheta/dz')
	plt.plotError(zVals,  lapThet,  label='d^2theta/dz^2')
	plt.legend()
	plt.show()


plotDiff = True
if (plotDiff):
	h1, = plt.plotError(zVals, diff, label='diff')
	h2, = plt.plotError(zVals, lapThet, label='lapalacian(theta)')
	h3, = plt.plotError(zVals, enthalpyAdvection, label='enthalpy advection')
	plt.legend()
	plt.show()
	
	
	
	
