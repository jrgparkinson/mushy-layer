import matplotlib.pyplot as plt
from OptimalStates.plotUtils import latexify, format_axes
import os
import numpy as np
from PltFile import PltFile

plt.close('all')
latexify(fig_width=10.0, fig_height=10.0)

figure_output_directory = '/home/parkinsonjl/mushy-layer/matlab/ChomboMatlab/'

pf = PltFile('/home/parkinsonjl/mushy-layer/matlab/ChomboMatlab/plt025000.2d.hdf5')

fig = plt.figure()
ax = plt.gca()

pf.load_data()


# Viridis reversed (so water is blue)
cmap = plt.get_cmap('viridis_r')

plot_levels = range(0,pf.num_levels)
# plot_levels = range(1,2)

print('Plotting porosity')
for level in plot_levels:

    porosity = pf.get_level_data('Porosity', level)

    X,Y =  pf.get_mesh_grid_xarray(porosity, True)
    img = ax.pcolormesh(X, Y, porosity, cmap=cmap)

# Draw level outlines
print('Plotting level outlines')
pf.plot_outlines(ax)


# make a color bar
cbar_porosity = plt.colorbar(img, cmap=cmap, ax=ax)
cbar_porosity.set_label(label='$\chi$', rotation=0, labelpad=10)

# Add streamlines
print('Plotting streamlines')
for level in plot_levels:
    u = pf.get_level_data('xAdvection velocity', level, valid_only=True)
    v = pf.get_level_data('yAdvection velocity', level, valid_only=True)

    X,Y =  pf.get_mesh_grid_xarray(u, grow=False)
    ax.streamplot(X, Y, u, v,  color='magenta')



# Temperature contours
for level in plot_levels:
    T = pf.get_level_data('Temperature', level, valid_only=True)
    X, Y = pf.get_mesh_grid_xarray(T)
    levels = np.linspace(-0.4, 1.4, 10)
    temperature_plt = ax.contour(X, Y, T, levels, linestyles='dashed', cmap='bwr')


ax.set_xlabel('$x$')
ax.set_ylabel('$z$', rotation=0, labelpad=10)

ax = format_axes(ax)
ax.set_aspect('equal', adjustable='box')


ax.set_xlim([0, 1.0])

plt.show()


figName = 'PlottingExample'
figure_full_path = os.path.join(figure_output_directory, figName + '.eps')
plt.savefig(figure_full_path,  format='eps')
print('Saved figure to %s' % figure_full_path)



