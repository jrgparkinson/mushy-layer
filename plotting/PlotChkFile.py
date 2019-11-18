import matplotlib.pyplot as plt
from util.plotUtils import latexify, format_axes
import os
import numpy as np
from ChkFile import ChkFile
from util.nondimensionalisation import dimensional_time
from util import shared_storage

plt.close('all')
latexify(fig_width=4.0, fig_height=3.0)

figure_output_directory = '/home/parkinsonjl/phd/figures/'

fig_name = 'chk088000'
# fig_name = 'chk277000'

base_folder = shared_storage.get_dir('SeaIceGrowth/T-20HigherPermeability/')
cf = ChkFile(os.path.join(base_folder, '%s.2d.hdf5' % fig_name))


fig = plt.figure()
ax = plt.gca()
# ax.set_position([0.05, 0.1, 0.7, 0.7])


cf.compute_diagnostic_vars()

# Viridis reversed (so water is blue)
cmap = plt.get_cmap('Blues')

plot_levels = range(0, cf.num_levels)
# plot_levels = range(1,2)

print('Plotting porosity')


porosity = cf.get_data('Porosity')

# X,Y =  cf.get_mesh_grid_xarray(porosity, True)
X,Y = cf.get_mesh_grid()
img = ax.pcolormesh(X, Y, porosity, cmap=cmap)


# make a color bar
cbar_porosity = plt.colorbar(img, cmap=cmap, ax=ax, shrink=0.75)
cbar_porosity.set_label(label='Porosity, $\chi$')

# Add streamlines
# print('Plotting streamlines')
# for level in plot_levels:
#     u = cf.get_data('xAdvection velocity')
#     v = cf.get_data('yAdvection velocity')
#
#     X,Y =  cf.get_mesh_grid_xarray(u, grow=False)
#     ax.streamplot(X, Y, u, v,  color='magenta')



# Temperature contours

T = cf.get_data('Temperature')
X, Y = cf.get_mesh_grid()
levels = np.linspace(-0.4, 1.4, 74)
delta_t = 20
t_ref = -23
T = T*delta_t + t_ref
levels = np.linspace(-5,5,21)
temperature_plt = ax.contour(X, Y, T, levels, linestyles='dashed', cmap='bwr')

cbar_temperature = plt.colorbar(temperature_plt, shrink=0.75)
cbar_temperature.set_label('$T$ ($^\circ$ C)')
cbar_temperature.set_ticks(levels[::2])
cbar_temperature.set_ticklabels(levels[::2])

ax.set_xlabel('$x$ (metres)')
ax.set_ylabel('$z$ (metres)')

ax = format_axes(ax)
ax.set_aspect('equal', adjustable='box')


#ax.set_xlim([0, 1.0])

time_days = dimensional_time(cf.time)
plt.title('$t = $ %d days' % time_days)

plt.tight_layout()

# Move axis left a bit
ax.set_position([-0.2, 0.15, 0.7, 0.75])



plt.draw()



figure_full_path = os.path.join(figure_output_directory, fig_name + '.jpg')
plt.savefig(figure_full_path,  format='jpg', dpi=800)
print('Saved figure to %s' % figure_full_path)

plt.show()



