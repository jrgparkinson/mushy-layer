import matplotlib
# matplotlib.use('Agg')

from PltFile import PltFile
import matplotlib.pyplot as plt


# file_loc = '/home/parkinsonjl/mnt/sharedStorage/optimalStates-restructured/Le200/CR2.000/RaC100/' \
#     'pts116-steady/CR2.0RaC100Le200ChiCubedPermeabilitypts116-033512.2d.hdf5'

# file_loc = '/home/parkinsonjl/mnt/sharedStorage/TransientGrowth/HeleShaw-256x256-CR2.0-Rm100.0-Rel0.2-TBottom1.1-np2/chk024200.2d.hdf5'
file_loc = '/home/parkinsonjl/mnt/sharedStorage/SeaIceWarming/SeaIce-128x256-Ttop-10.0-Tbottom-2.8-d0.002-NewTtop-2laterRestart-np8-saltyTop-K01e-09-finished-restart-0/chk334800.2d.hdf5'
# file_loc = '/home/parkinsonjl/mnt/sharedStorage/SeaIceWarming/SeaIce-128x256-Ttop-10.0-Tbottom-2.8-d0.002-NewTtop-2laterRestart-np8-saltyTop-K01e-09-finished-restart-0/plt375860.2d.hdf5'
# file_loc = '/home/parkinsonjl/mnt/sharedStorage/TransientGrowth/HeleShaw-256x256-CR2.0-Rm100.0-Rel0.2-TBottom1.1-np2/plt029540.2d.hdf5'
p = PltFile(file_loc)

p.load_data()

p.compute_diagnostic_vars()

# [num_channels, rel_pos, channel_spacing] = p.channel_properties(False)

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

v = p.get_level_data('xDarcy velocity', 0)
img_v = ax_list[1].pcolormesh(x, y, v)


# plt.title('Concentration ratio = %s' % (p.inputs['parameters.compositionRatio'] - 1.0))

# plt.savefig('tex_demo')
plt.show()

