from PltFile import PltFile
import matplotlib.pyplot as plt

file_loc = '/home/parkinsonjl/mnt/sharedStorage/optimalStates-restructured/Le200/CR2.000/RaC100/' \
    'pts116-steady/CR2.0RaC100Le200ChiCubedPermeabilitypts116-033512.2d.hdf5'
p = PltFile(file_loc)

[num_channels, rel_pos, channel_spacing] = p.channel_properties(False)

lev0Dx = p.levels[0][p.DX]
print('dx = ' + str(lev0Dx))
domainSize = p.domainSize
print('Domain extent: ' + str(domainSize))
# av_chan_spacing = np.mean(channelSpacing)
av_chan_spacing = (p.domainSize[2] - p.domainSize[0]) / float(num_channels)

print('Num channels: ' + str(num_channels) + ', relative positions (0 = left, 1 = right): ' + str(rel_pos))
print('Average channel spacing: ' + str(av_chan_spacing))
p.plot_field('Porosity')

plt.title(r"$\mathscr{C} = 1.0$")

plt.savefig('tex_demo')
plt.show()
