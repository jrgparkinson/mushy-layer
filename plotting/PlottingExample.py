import matplotlib.pyplot as plt
import os
import numpy as np
from PltFile import PltFile
import sys

plt.close('all')
# latexify(fig_width=10.0, fig_height=10.0)

# 2D plotting example

figure_output_directory = os.path.join(os.environ['MUSHY_LAYER_DIR'], 'plotting')
figure_name = 'PlottingExample'

# data_file =  os.path.join(os.environ['MUSHY_LAYER_DIR'], 'plotting', 'plt001142.2d.hdf5')
# data_file =  os.path.join(os.environ['MUSHY_LAYER_DIR'], 'plotting', 'plt004312.3d.hdf5')
data_file = os.path.join(os.environ['MUSHY_LAYER_DIR'], 'matlab', 'ChomboMatlab', 'plt025000.2d.hdf5')

# Get data
pf = PltFile(data_file)
pf.load_data()

# Start making plot
fig = plt.figure()
ax = plt.gca()

# Viridis reversed (so water is blue)
cmap = plt.get_cmap('viridis_r')

plot_levels = range(0, pf.num_levels)
# plot_levels = range(1,2)

if not len(plot_levels):
    print('No data to plot')
    sys.exit(-1)

print('Plotting porosity')
img = None
for level in plot_levels:
    porosity = pf.get_level_data('Porosity', level)

    # Convert to python indexing
    # porosity = porosity.transpose('z', 'y', 'x')

    # porosity.coords

    # X,Y =  pf.get_mesh_grid_xarray(porosity, True)

    X = porosity.coords['x']
    Y = porosity.coords['y']

    # porosity = np.array(porosity).transpose()

    img = ax.pcolormesh(X, Y, porosity, cmap=cmap)

# plt.show()
# Draw level outlines
print('Plotting level outlines')
pf.plot_outlines(ax)

# make a color bar
cbar_porosity = plt.colorbar(img, cmap=cmap, ax=ax)
cbar_porosity.set_label(label=r'$\chi$', rotation=0, labelpad=10)

# Add streamlines
print('Plotting streamlines')
for level in plot_levels:
    u = pf.get_level_data('xAdvection velocity', level, valid_only=True)
    v = pf.get_level_data('yAdvection velocity', level, valid_only=True)

    # X,Y =  pf.get_mesh_grid_xarray(u, grow=False)

    x = np.array(u.coords['x'])
    y = np.array(u.coords['y'])

    # ax.streamplot(X, Y, u, v,  color='magenta')
    ax.streamplot(x, y, u, v, color='magenta')

# Temperature contours
for level in plot_levels:
    T = pf.get_level_data('Temperature', level, valid_only=True)
    # X, Y = pf.get_mesh_grid_xarray(T)

    X = T.coords['x']
    Y = T.coords['y']

    levels = np.linspace(-0.4, 1.4, 10)
    temperature_plt = ax.contour(X, Y, T, levels, linestyles='dashed', cmap='bwr')

ax.set_xlabel('$x$')
ax.set_ylabel('$z$', rotation=0, labelpad=10)

# ax = format_axes(ax)
ax.set_aspect('equal', adjustable='box')

# ax.set_xlim([0, 1.0])

plt.show()

figure_full_path = os.path.join(figure_output_directory, figure_name + '.eps')
plt.savefig(figure_full_path, format='eps')
print('Saved figure to %s' % figure_full_path)
