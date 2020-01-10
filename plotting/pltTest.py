# import matplotlib
# matplotlib.use('Agg')

from PltFile import PltFile
import matplotlib.pyplot as plt

file_loc = '/home/parkinsonjl/mnt/sharedStorage/optimalStates-restructured/Le200/CR2.000/RaC100/' \
    'pts116-steady/CR2.0RaC100Le200ChiCubedPermeabilitypts116-033512.2d.hdf5'

file_loc = 'plt001142.2d.hdf5'

# file_loc = shared_storage.get_dir('SeaIceGrowth/T-10Outflow/chk135000.2d.hdf5')
p = PltFile(file_loc, load_data=True)
# p.compute_diagnostic_vars()

lev0Dx = p.levels[0][p.DX]
print('dx = ' + str(lev0Dx))
domainSize = p.domain_size
print('Domain extent: ' + str(domainSize))

fig, axes = plt.subplots(2, 2)
ax_list = axes.flatten()

# print('Num channels: ' + str(num_channels) + ', relative positions (0 = left, 1 = right): ' + str(rel_pos))
# print('Average channel spacing: ' + str(av_chan_spacing))
porosity = p.get_level_data('Porosity', 0)
x = porosity.coords['x']
y = porosity.coords['y']
img = ax_list[0].pcolormesh(x, y, porosity)

# v = p.get_level_data('xDarcy velocity', 0)
# img_v = ax_list[1].pcolormesh(x, y, v)

# plt.title('Concentration ratio = %s' % (p.inputs['parameters.compositionRatio'] - 1.0))

# plt.savefig('tex_demo')
plt.show()

