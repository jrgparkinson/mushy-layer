"""

This script is for testing that the mushy-layer code correctly solves diffusion problems correctly
with all sorts of boundary conditions
Within this file there are a number of useful functions, and a 'DiffusiveSolution' class,
then at the bottom of we run everything which will:
  a) run an appropriate simulation using the mushy layer code to compute the temperature T(z, t)
  b) solve the steady-state problem using the solve_bvp function in scipy.integrate, to find T_{python}(z)
  c) plot the two solutions along with the error in the final mushy-layer simulation
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_bvp
import os
import matplotlib as mpl
import cycler
from pathlib import Path
from chombopy.inputs import read_inputs, write_inputs
import subprocess
import shutil
from plotting.MushyPltFile import MushyPltFile, latexify2
from mushyLayerRunUtils import get_mushy_layer_dir, get_executable_name


class DiffusiveSolution:
    a = 0
    b = 0
    T_ref = 0
    F = 0

    def __init__(self, min_temperature, max_temperature, method='Linear'):

        self.min_temperature = min_temperature
        self.max_temperature = max_temperature

        self.method = method

        self.V = np.nan
        self.cr = np.nan
        self.st = np.nan
        self.a = np.nan
        self.b = np.nan
        self.T_ref = np.nan
        self.F = np.nan

    def set_nonlinear_bcs(self, a, b, reference_temperature, flux, frame_advection_velocity, conc_ratio, stefan):
        self.a = a
        self.b = b
        self.T_ref = reference_temperature
        self.F = flux
        self.V = frame_advection_velocity
        self.cr = conc_ratio
        self.st = stefan

    # noinspection PyUnresolvedReferences
    def compute_solution(self, z_sol):

        if self.method == 'Linear':
            if self.V == 0:
                solution = self.max_temperature - z_sol
            else:
                print('Unable to compute solution for method="linear" and V!=0')
                sys.exit(-1)

        elif self.method == 'Mixed':
            z_solve = np.linspace(0, 1, 10)
            init_guess = np.zeros((2, z_solve.size))
            res_a = solve_bvp(self.diffusion_equation_function, self.diffusion_eq_bc, z_solve, init_guess)
            solution = res_a.sol(z_sol)[0]

        else:
            print('Unknown method "%s"' % self.method)
            sys.exit(-1)

        return solution


    def diffusion_equation_function(self, x, temperature):
        # For solving d^T/dz^2 - V*(dT/dz + dchi/dz) = 0

        # Define T_0 = T, T_1 = dT_0/dz
        # Then,
        # T_0' = T_1
        # T_1' = V*(T_1+dchi/dz)

        # Could in principle compute porosity in here if necessary
        # temperature = T[0]
        # porosity = self.cr/(temperature = self.cr)

        # Assuming either liquid or mushy

        chi = np.minimum(1.0, self.cr / (self.cr - temperature[0]))

        # if chi < 1:
        dchi_dz = temperature[1] * (self.cr / np.power(self.cr - temperature[0], 2))
        # else:
        dchi_dz[chi == 1.0] = 0.0

        # dchi_dz = 0

        return np.vstack((temperature[1],
                          self.V * (temperature[1] + self.st * dchi_dz)))

    def diffusion_eq_bc(self, temperature_left_boundary, temperature_right_boundary):
        # Current BCs:
        # T(z=0) = min_temperature
        # T(z=1) = max_temperature

        # First component of Ta, Tb is T
        # Second component if T'

        #

        # fixed_bc_res = np.array([Ta[0]-self.max_temperature,
        #                          Tb[0]-self.min_temperature])

        mixed_bc_res = np.array([temperature_left_boundary[0] - self.max_temperature,
                                 self.a * temperature_right_boundary[1] - self.F
                                 - self.b * (temperature_right_boundary[0] - self.T_ref)])

        return mixed_bc_res


def base_inputs():
    test_dir = os.path.join(get_mushy_layer_dir(), 'test/diffusiveGrowth')
    default_inputs = read_inputs(os.path.join(test_dir, 'inputs'))

    default_inputs['main.min_time'] = 0.5
    default_inputs['main.max_time'] = 3.0
    default_inputs['main.steady_state'] = 1e-8

    return default_inputs


if __name__ == "__main__":

    # Ensure we can find HDF5 libraries to run mushy layer code
    hdf5_library_path = '/home/parkinsonjl/soft/hdf5-1.8.21p-new/lib/'
    if 'LD_LIBRARY_PATH' in os.environ:
        os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + ':' + hdf5_library_path
    else:
        os.environ['LD_LIBRARY_PATH'] = hdf5_library_path

    # Chose some parameter options:

    # Options
    # T_min = 4.0
    # T_max = 5.0
    st = 1.0  # stefan number
    cr = 0.5

    # plotting_field = "Temperature"
    include_solution_diffusion = False

    plotting_field = "Bulk concentration"
    include_solution_diffusion = True

    # this is just a dirichlet BC, data matches well
    # opt = {'a': 0.0, 'b': 1.0, 'Tref': 4.0, 'F': 0.0, 'T_min': T_min, 'T_max': T_max, 'V': 1.0, 'CR': cr}

    # this is no flux
    # opt = {'a': 1.0, 'b': 0.0, 'Tref': 4.0, 'F': 0.0, 'T_min': T_min, 'T_max': T_max, 'V': 1.0, 'CR': cr}

    # this is a fixed flux (negative for cooling)
    # opt = {'a': 1.0, 'b': 0.0, 'Tref': 4.0, 'F': -1.0, 'T_min': T_min, 'T_max': T_max, 'V': 1.0, 'CR': cr}

    # this is a fixed flux and temperature
    # opt = {'a': 0.5, 'b': -1.0, 'Tref': 4.0, 'F': -0.5, 'T_min': T_min, 'T_max': T_max, 'V': 1.0, 'CR': cr}

    # this is a fixed flux and temperature
    # opt = {'a': 0.5, 'b': -1.0, 'Tref': 4.5, 'F': -1.6, 'T_min': T_min, 'T_max': T_max, 'V': 1.0, 'CR': cr}

    # Mushy layer solidification:
    T_min = -0.8
    T_max = 0.2

    # T_min = 2.0
    # T_max = 3.0
    opt = {'a': 0.5, 'b': -1.0, 'Tref': T_min, 'F': -1.0, 'T_min': T_min,
           'T_max': T_max, 'V': 1.0, 'CR': cr}  # this is a fixed flux and temperature
    # opt = {'a': 0.0, 'b': 1.1, 'Tref': T_min, 'F': 0.0, 'T_min': T_min,
    #        'T_max': T_max, 'V': 1.0, 'CR': cr}  # dirichlet

    # All simulations will be put in their own sub directory within this directory:
    base_dir = Path.cwd().parent.parent / 'execSubcycle' / 'ImperfectCooling'

    # a = 0.5
    # b = -1.0
    # Tref = 4.0
    # F = -0.5

    analytic_solution = DiffusiveSolution(T_min, T_max, method='Mixed')

    # this is just a dirichlet BC, data matches well
    analytic_solution.set_nonlinear_bcs(a=opt['a'], b=opt['b'], reference_temperature=opt['Tref'], flux=opt['F'],
                                        frame_advection_velocity=opt['V'], conc_ratio=cr, stefan=st)

    # Now let's run the mushy-layer simulation
    sim_name = 'ImperfectCooling_a%s_b%s_Tref%s_F%s_CR%s_st%s' % (opt['a'], opt['b'], opt['Tref'], opt['F'], cr, st)

    full_output_folder = os.path.join(base_dir, sim_name)
    figure_output_paths = [os.path.join(full_output_folder, 'errorComparison'),
                           os.path.join(base_dir, sim_name)]

    if os.path.exists(full_output_folder):
        print('Output directory already exists: %s' % full_output_folder)
        answer = input('Remove it? (y/n) ')
        if answer == 'y':
            shutil.rmtree(full_output_folder)
        else:
            print('Exiting')
            sys.exit(-1)

    os.makedirs(full_output_folder)

    inputs = base_inputs()
    inputs['parameters.stefan'] = st
    inputs['bc.temperatureHiVal'] = [0, opt['F']]
    inputs['bc.aHi'] = [0, opt['a']]
    inputs['bc.bHi'] = [0, opt['b']]
    inputs['bc.TRefHi'] = [0, opt['Tref']]
    inputs['bc.enthalpyLoVal'] = [0.0, T_max + st]
    inputs['bc.enthalpyHiVal'] = [0.0, T_max + st]
    inputs['parameters.nonDimVel'] = opt['V']
    inputs['parameters.compositionRatio'] = cr
    inputs['bc.HC_noflux'] = True

    if include_solution_diffusion:
        inputs['parameters.lewis'] = 10.0

    inputs['main.output_folder'] = full_output_folder

    new_inputs_loc = os.path.join(full_output_folder, 'inputs')
    write_inputs(new_inputs_loc, inputs)

    exec_loc = get_executable_name(return_full_path=True)

    cmd = 'cd %s; %s inputs' % (full_output_folder, exec_loc)

    print(cmd)
    print('Running... (usually takes ~ 30 seconds) ')

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(process.returncode)

    print('Finished')

    # Run comparison
    plt_files = [f for f in os.listdir(full_output_folder) if ('LinearDiffusion' in f and '.2d.hdf5' in f)]
    plt_files = sorted(plt_files)
    if not plt_files:
        print('No data to plot in folder %s' % full_output_folder)
        sys.exit(-1)

    # get some initial data
    init_pf = MushyPltFile(os.path.join(full_output_folder, plt_files[0]))
    init_pf.load_data()
    T_data = init_pf.get_level_data(plotting_field).mean('x').squeeze()
    z = T_data.coords['y']
    T_python = analytic_solution.compute_solution(z)

    # Set color cycle from colormap
    n = len(plt_files)
    color = plt.cm.viridis(np.linspace(0.1, 0.9, n))  # This returns RGBA; convert:
    mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)

    times = []

    latexify2(5.5, 3.5)

    fig = plt.figure()

    for plot_loc in plt_files:
        pf = MushyPltFile(os.path.join(full_output_folder, plot_loc))
        pf.load_data()
        T_data = pf.get_level_data(plotting_field).mean('x').squeeze()
        z = T_data.coords['y']

        plt.plot(z, T_data, label='')

        times.append(float('%.1g' % pf.time))

    if plotting_field == "Temperature":
        plt_python = plt.plot(z, T_python, label='Python', color='red', linestyle='--')

    ax = plt.gca()

    ax.set_xlabel('$z$')
    ax.set_ylabel(plotting_field)

    # Add colorbar for times
    cax_porosity = fig.add_axes([0.2, 0.86, 0.5, 0.03])
    cmap = plt.get_cmap('viridis')

    norm = mpl.colors.Normalize(vmin=min(times), vmax=max(times))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array(np.array([]))
    cbar = plt.colorbar(sm, cax=cax_porosity, ticks=[min(times), max(times)], orientation='horizontal')
    cax_porosity.xaxis.set_ticks_position('top')
    cax_porosity.xaxis.set_label_position('top')
    cbar.set_label('$t$', labelpad=-3)

    ax.set_xlim([0, 1])

    # Plot error
    if plotting_field == "Temperature":
        difference = T_data - T_python
        ax_err = ax.twinx()
        plt_err = ax_err.plot(z, difference, color='black', linestyle=':', label='Error (right axis)')
        ax_err.set_ylabel('Error')

        # Add legend
        lns = plt_python + plt_err
        labs = [line.get_label() for line in lns]
        ax.legend(lns, labs, loc=0)

        ax_err.set_frame_on(True)  # make sure there is any background
        ax_err.set_position([0.13, 0.15, 0.65, 0.7])

    # Sort axis positioning
    ax.set_zorder(1)  # make it on top
    ax.set_frame_on(False)  # make it transparent

    # Sort figure positioning
    ax.set_position([0.13, 0.15, 0.65, 0.7])

    # Save figure
    for figure_output_path in figure_output_paths:
        fig_name = figure_output_path + '.eps'
        print('Saving figure to %s' % fig_name)
        plt.savefig(fig_name, format='eps')

    plt.show()
