import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from nondimensionalisation import dimensional_salinity, dimensional_temperature, dimensional_time
from plotting.MushyPltFile import latexify2
import matplotlib as mpl
import re
import math


def get_profiles(folder):
    nc_files = [f for f in os.listdir(folder) if re.match(r'^chk.*\.nc$', f)]

    full_paths = [os.path.join(folder, x) for x in nc_files]
    dfs = [xr.open_dataset(f) for f in full_paths]

    ds = xr.merge(dfs)

    return ds


def get_color(colors_arr, min_t, max_t, t):
    interp = (t - min_t) / (max_t - min_t)

    index = int(math.floor(len(colors_arr) * interp))

    return colors_arr[index]


def get_dimensional(profs, time, prof_name):
    prof = profs[prof_name].sel(t=time)
    if 'Salinity' in prof_name:
        prof = dimensional_salinity(prof)
    elif 'Temperature' in prof_name:
        prof = dimensional_temperature(prof)
    return prof


# Define these
control_folder = ''
warming_folder = ''

profiles_control = get_profiles(control_folder)
profiles_warming = get_profiles(warming_folder)

print(profiles_control)

rm_profile_control = profiles_control['Rm profile']
temperature_profile_control = profiles_control['Temperature']

c = rm_profile_control.coords
times_control = np.array(rm_profile_control.coords['t'])
z = rm_profile_control.coords['z']

time_warming = np.array(profiles_warming['Porosity'].coords['t'])

dim_times_control = dimensional_time(times_control)
dim_times_warming = dimensional_time(times_control)

valid_control_times = times_control[times_control > time_warming[0]]

latexify2(3, 3)
fig2, axes = plt.subplots(2, 1)
# num_steps = int(float(len(times))/float(plot_interval))
num_steps = len(times_control)
print('Num timesteps: %d' % num_steps)
colors = plt.cm.viridis(np.linspace(0, 1, num_steps + 1))

# profile_name = 'Harmonic permeability (min horiz) profile'
# profile_name = 'Porosity'
# profile_name= 'Temperature'
# profile_name = 'LiquidSalinity'
# profile_name = 'Porosity'
# profile_name = 'BulkSalinity'
# profile_name = 'Rm profile'
# profile_name = 'Harmonic permeability (min horiz) profile'
# profile_name = 'Rm profile no z dependence'
# profile_name = 'Permeability'
profile_name = 'Min Vertical Velocity'

comparison = False

min_time = time_warming[0]
max_time = time_warming[-1]

plot_interval = 10

plot_interval_time = 0.001

# last_plot_

if comparison:

    for i, time in enumerate(times_control):

        if time < min_time or time > max_time:
            continue

        if not i % plot_interval == 0:
            continue

        rm_time = rm_profile_control.sel(t=time)
        temperature_time = temperature_profile_control.sel(t=time)

        profile = get_dimensional(profiles_control, time, profile_name)
        this_color = get_color(colors, min_time, max_time, time)

        axes[1].plot(profile, z, linestyle=':', label='$t=%.5g$' % dim_times_control[i], color=this_color)

        i = i + 1

plot_interval = 10
for i in range(30, len(time_warming) - 1):

    if not i % plot_interval == 0:
        continue

    time = time_warming[i]

    if profile_name not in profiles_warming.keys():
        continue

    profile = get_dimensional(profiles_warming, time, profile_name)
    this_color = get_color(colors, min_time, max_time, time)

    axes[1].plot(profile, z, color=this_color)  # label='$t=%.5g$' % dim_times_warming[i],

    i = i + 1

axes[1].set_ylim([0.6, 1.0])
# axes[1].set_xlim([0.0, 31])

# axes[0].title(profile_name)
axes[1].set_ylabel('$z$')
axes[1].set_xlabel(profile_name)

norm = mpl.colors.Normalize(vmin=dimensional_time(min_time), vmax=dimensional_time(max_time))
cb1 = mpl.colorbar.ColorbarBase(axes[0], cmap=mpl.cm.viridis, orientation='horizontal', norm=norm)
cb1.set_label('Time (days)')
axes[0].set_position([0.15, 0.9, 0.8, 0.06])
axes[1].set_position([0.15, 0.14, 0.8, 0.6])

# plt.legend()


plt.draw()

figure_output_directory = '/home/parkinsonjl/phd/figures/'
fig_name = '%s_profile_warming' % profile_name

if comparison:
    fig_name = fig_name + '_comparison'

figure_full_path = os.path.join(figure_output_directory, fig_name + '.jpg')
plt.savefig(figure_full_path, format='jpg', dpi=800)
print('Saved figure to %s' % figure_full_path)

plt.show()
