from PltFile import PltFile
import numpy as np
#from matplotlib import mpl,pyplot

p = PltFile('channelSpacing-periodic-CR1.25RaC400Le200ChiCubedPermeabilitypts1024-024400.2d.hdf5')


[numChannels, relPos, channelSpacing] = p.channelProperties(False)

lev0Dx = p.levels[0][p.DX]
print('dx = ' + str(lev0Dx))
domainSize = p.domainSize
print('Domain extent: ' + str(domainSize))
avChanSpacing = np.mean(channelSpacing)
avChanSpacing = (p.domainSize[2] - p.domainSize[0])/float(numChannels)

print('Num channels: ' + str(numChannels) + ', relative positions (0 = left, 1 = right): ' + str(relPos))
print('Average channel spacing: ' + str(avChanSpacing))     
p.plotField('Porosity')

