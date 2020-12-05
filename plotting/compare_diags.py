import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

"""
This is a short script to quickly plot all the diagnostics from a simulation
"""


def get_regrid_times(pout):
    regrid_times = []
    with open(pout, 'r') as file:
        lines = file.read()

        # matches = re.findall(r"coarse time step\s+(\d+)\s+old time =\s+([\.\da-z-]+)\s+old.*\nAMR\:\:regrid", lines)
        matches = re.findall(r"New grids are different \(time = (.*)\)\? (\d)", lines)
        # print(matches)

        regrid_times = [float(match[0]) for match in matches if int(match[1]) == 1]

    return regrid_times


# Options to change:
base_path = '/home/parkinsonjl/mushy-layer/execSubcycle/freestream/'
folders = [os.path.join(base_path, f) for f in os.listdir(base_path) if
         os.path.exists(os.path.join(base_path, f, 'diagnostics.csv'))]

folders = [os.path.join(base_path, 'uniform'),
         os.path.join(base_path, 'amr-grid-buffer-4'),
         os.path.join(base_path, 'amr-newregrid'),
         os.path.join(base_path, 'amr-regrid_eta_small'),
         os.path.join(base_path, 'amr-newregrid-smallregrideta')]  # a

folders = [os.path.join(base_path, 'uniform'),
         os.path.join(base_path, 'amr-newregrid'),
         os.path.join(base_path, 'amr-newregrid-initpressure'),
         os.path.join(base_path, 'amr-oldregrid')]  # a
# folders = []


diag_x_axis = 'time'
diag_to_compare = 'Fs_vertical_av'
diag_to_compare = 'L0FsVertFluid'
# diag_to_compare = 'LambdaPostRegrid'
diag_to_compare = 'LambdaMax'
# diag_to_compare = 'maxFreestreamCorrection'
# Load diagnostics

fig = plt.figure()
ax = fig.gca()

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

for (folder, color) in zip(folders, colors[:len(folders)]):
    diags_file = os.path.join(folder, 'diagnostics.0.csv')
    nskip = 0
    with open(diags_file, 'r') as f:
        lines = f.readlines()
        if lines[0] == lines[1]:
            nskip = 1
    diags = pd.read_csv(diags_file, skiprows=nskip)

    label = diags_file.split('/')[-2]

    if diag_to_compare in diags.keys():
        diag_vals = np.array(diags[diag_to_compare])
        x_vals = np.array(diags[diag_x_axis])

        # remove nans
        x_vals = x_vals[~np.isnan(diag_vals)]
        diag_vals = diag_vals[~np.isnan(diag_vals)]

        ax.plot(x_vals, diag_vals, label=label, color=color)



lims = ax.get_ylim()
num_lines = len(folders)
lims_length = (lims[1]-lims[0])/num_lines

for i in range(0, len(folders)):
# for (folder, color) in zip(folders, colors[:len(folders)]):

    color = colors[i]
    regrid_y_lims = [lims[0] + i*lims_length,
                     lims[0] + (i+1) * lims_length]

    regrid_times = get_regrid_times(os.path.join(folders[i], 'pout.0'))

    for time in regrid_times:
        ax.plot([time, time],regrid_y_lims, linestyle='--', color=color)

ax.legend()
ax.set_title(diag_to_compare)

plt.savefig(diag_to_compare + '.pdf', format='pdf', dpi=800)
plt.show()
