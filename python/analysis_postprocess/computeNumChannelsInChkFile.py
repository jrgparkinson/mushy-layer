from ChkFile import ChkFile
#import numpy as np
#from matplotlib import mpl,pyplot

c = ChkFile('/home/parkinsonjl/mnt/sharedStorage/wettlaufer/CR1.4375RaC1000000.0Le200KozenyPermeabilityDa0.0001R1e-45pts512-S7.0-TB-20.0-domWidth0.5-0/chkfixedChill-periodic-CR1.4375RaC1000000.0Le200KozenyPermeabilityDa0.0001R1e-45pts512-S7.0-TB-20.0-domWidth0.5-014700.2d.hdf5')

[numChannels, relPos] = c.channelProperties(True)  

print('Num channels: ' + str(numChannels) + ', relative positions (0 = left, 1 = right): ' + str(relPos))     

