import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

"""
This is a short script to quickly plot all the diagnostics from a simulation
"""

# Options to change:
diags_file = os.path.join('..', 'execSubcycle', '3DBug', 'diagnostics_backup.csv')
diag_x_axis = 'timestep'

# Load diagnostics
diags = pd.read_csv(diags_file)

# These diagnostic are reported at multiple locations in the domain - want to plot on the same set of axes
special_diags = {'horizontallyAveragedSalinity_': ['0', '20', '40', '60'], 'Fs_vert_': ['10pc', '20pc', '30pc', '40pc', '50pc']}
special_diags_names = list(special_diags.keys())
unique_diags = [d for d in diags.keys() if not (sum([d in s for s in special_diags_names]) or d == diag_x_axis) ]

# Create figure
rows = int(np.floor(np.sqrt(len(unique_diags))))
cols = int(np.ceil(len(unique_diags)/rows))
fig, axes= plt.subplots(rows, cols, sharex=True)
ax_list = axes.flatten()

# Plot each diagnostic
for d, ax in zip(unique_diags, ax_list):
    ax.plot(diags[diag_x_axis], diags[d])
    ax.set_title(d)

fig.subplots_adjust(left=0.06, right=0.95, bottom=0.05, top=0.95, hspace=0.4, wspace=0.5)

# Second figure for diagnostics sampled in multiple places
fig, axes= plt.subplots(1, len(special_diags_names), sharex=True)
ax_list = axes.flatten()
for d, ax in zip(special_diags_names, ax_list):
    for suffix in special_diags[d]:
        ax.plot(diags[diag_x_axis], diags[d + suffix], label=suffix)
        ax.set_title(d)

    ax.legend()

plt.show()
