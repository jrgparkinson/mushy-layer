from PltFile import PltFile
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_bvp
import os
import matplotlib as mpl
import cycler
import mushyLayerRunUtils
import subprocess
import sys
import shutil


# This script is for testing that the mushy-layer code correctly solves diffusion problems correctly
# with all sorts of boundary conditions
# Within this file there are a number of useful functions, and a 'DiffusiveSolution' class,
# then at the bottom of we run everything which will:
#  a) run an appropriate simulation using the mushy layer code to compute the temperature T(z, t)
#  b) solve the steady-state problem using the solve_bvp function in scipy.integrate, to find T_{python}(z)
#  c) plot the two solutions along with the error in the final mushy-layer simulation

def latexify(fig_width=5.0, fig_height=4.0):

    font_size = 12

    params = {'backend': 'ps',
              'text.latex.preamble': ['\\usepackage{gensymb}', '\\usepackage{mathrsfs}'],
              'axes.labelsize': font_size,  # fontsize for x and y labels (was 10)
              'axes.titlesize': font_size,
              'legend.fontsize': font_size,  # was 10
              'xtick.labelsize': font_size,
              'ytick.labelsize': font_size,
              'lines.markersize': 3,
              'lines.linewidth': 1,
              'text.usetex': True,
              'figure.figsize': [fig_width, fig_height],
              'font.family': 'serif'
              }

    mpl.rcParams.update(params)


class DiffusiveSolution:

    a = 0
    b = 0
    T_ref = 0
    F = 0


    def __init__(self, min_temperature, max_temperature, method='Linear'):

        self.min_temperature = min_temperature
        self.max_temperature = max_temperature

        self.method = method

    def set_nonlinear_bcs(self, a, b, T_ref, F, V, cr, st):
        self.a = a
        self.b =b
        self.T_ref = T_ref
        self.F = F
        self.V = V
        self.cr = cr
        self.st = st

    def compute_solution(self, z):

        # z = np.linspace(self.min_z, self.max_z, 100)

        if self.method == 'Linear':

            if self.V == 0:

                solution = self.max_temperature - z



            # solution = z

        elif self.method == 'FixedFlux':

            pass

        elif self.method == 'Mixed':

            # init_guess = self.max_temperature - z

            # Ensure z solve goes to 0 and 1 so BCs are implemented properly
            z_solve = np.linspace(0, 1, 10)

            init_guess = np.zeros((2, z_solve.size))

            res_a = solve_bvp(self.diffusion_equation_function, self.diffusion_eq_bc, z_solve, init_guess)
            solution = res_a.sol(z)[0]

        return solution

    def diffusion_equation_function(self, x, T):
        # For solving d^T/dz^2 - V*(dT/dz + dchi/dz) = 0

        # Define T_0 = T, T_1 = dT_0/dz
        # Then,
        # T_0' = T_1
        # T_1' = V*(T_1+dchi/dz)

        # Could in principle compute porosity in here if necessary
        # temperature = T[0]
        # porosity = self.cr/(temperature = self.cr)

        # Assuming either liquid or mushy

        chi = np.minimum(1.0, self.cr/(self.cr-T[0]))

        #if chi < 1:
        dchi_dz = T[1]*(self.cr/np.power(self.cr-T[0], 2))
        #else:
        dchi_dz[chi==1.0] = 0.0

        # dchi_dz = 0

        return np.vstack((T[1],
                          self.V * (T[1] + self.st*dchi_dz)))

    def diffusion_eq_bc(self, Ta, Tb):
        # Current BCs:
        # T(z=0) = min_temperature
        # T(z=1) = max_temperature

        # First component of Ta, Tb is T
        # Second component if T'

        #

        fixed_bc_res = np.array([Ta[0]-self.max_temperature,
                                 Tb[0]-self.min_temperature])

        mixed_bc_res = np.array([Ta[0]-self.max_temperature,
                                 self.a*Tb[1] - self.F - self.b*(Tb[0]-self.T_ref)])

        return mixed_bc_res


def base_inputs():

    test_dir = os.path.join(mushyLayerRunUtils.get_mushy_layer_dir(), 'test/diffusiveGrowth')
    inputs = mushyLayerRunUtils.read_inputs(os.path.join(test_dir, 'inputs'))

    inputs['main.min_time'] = 0.5
    inputs['main.max_time'] = 3.0
    inputs['main.steady_state'] = 1e-8


    return inputs

if __name__ == "__main__":


    # Ensure we can find HDF5 libraries to run mushy layer code
    if 'LD_LIBRARY_PATH' in os.environ:
        os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + ':/home/parkinsonjl/soft/hdf5-1.8.21p-new/lib/'
    else:
        os.environ['LD_LIBRARY_PATH'] = '/home/parkinsonjl/soft/hdf5-1.8.21p-new/lib/'


    # Chose some parameter options:

    # Options
    T_min = 4.0
    T_max = 5.0
    st = 1.0 # stefan number
    cr = 0.5

    opt = {'a': 0.0, 'b': 1.0, 'Tref': 4.0, 'F': 0.0, 'T_min': T_min, 'T_max': T_max, 'V': 1.0, 'CR': cr}  # this is just a dirichlet BC, data matches well
    opt = {'a': 1.0, 'b': 0.0, 'Tref': 4.0, 'F': 0.0, 'T_min': T_min, 'T_max': T_max, 'V': 1.0, 'CR': cr}  # this is no flux
    opt = {'a': 1.0, 'b': 0.0, 'Tref': 4.0, 'F': -1.0, 'T_min': T_min, 'T_max': T_max, 'V': 1.0, 'CR': cr} # this is a fixed flux (negative for cooling)
    opt = {'a': 0.5, 'b': -1.0, 'Tref': 4.0, 'F': -0.5, 'T_min': T_min, 'T_max': T_max, 'V': 1.0, 'CR': cr} # this is a fixed flux and temperature
    opt = {'a': 0.5, 'b': -1.0, 'Tref': 4.5, 'F': -1.6, 'T_min': T_min, 'T_max': T_max, 'V': 1.0, 'CR': cr}  # this is a fixed flux and temperature


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
    base_dir = '/home/parkinsonjl/mushy-layer/execSubcycle/ImperfectCooling/'

    # a = 0.5
    # b = -1.0
    # Tref = 4.0
    # F = -0.5

    analytic_solution = DiffusiveSolution(T_min, T_max, method='Mixed')
    analytic_solution.set_nonlinear_bcs(a=opt['a'], b=opt['b'], T_ref=opt['Tref'], F=opt['F'], V=opt['V'], cr=cr, st=st)  # this is just a dirichlet BC, data matches well

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

    inputs['main.output_folder'] = full_output_folder

    new_inputs_loc = os.path.join(full_output_folder, 'inputs')
    mushyLayerRunUtils.write_inputs(new_inputs_loc, inputs)

    exec_loc = mushyLayerRunUtils.get_executable_name(return_full_path=True)

    cmd = 'cd %s; %s inputs' % (full_output_folder, exec_loc)

    print(cmd)
    print('Running... (usually takes ~ 30 seconds) ')

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print (process.returncode)

    print('Finished')

    # Run comparison
    plt_files = [f for f in os.listdir(full_output_folder) if ('LinearDiffusion' in f and '.2d.hdf5' in f)]

    plt_files = sorted(plt_files)

    # Set color cycle from colormap
    n = len(plt_files)
    color = plt.cm.viridis(np.linspace(0.1, 0.9, n))  # This returns RGBA; convert:
    mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)

    times = []

    latexify(5.5, 3.5)

    fig  = plt.figure()

    for plot_loc in plt_files:

        pf = PltFile(os.path.join(full_output_folder, plot_loc))
        pf.load_data()
        T_data = pf.get_level_data('Temperature').mean('x').squeeze()
        z = T_data.coords['y']

        T_python = analytic_solution.compute_solution(z)

        plt.plot(z, T_data, label='')

        times.append(float('%.1g' % pf.time))


    plt_python = plt.plot(z, T_python, label='Python', color='red', linestyle='--')

    ax = plt.gca()

    ax.set_xlabel('$z$')
    ax.set_ylabel('$T$')


    # Add colorbar for times
    cax_porosity = fig.add_axes([0.2, 0.86, 0.5, 0.03])
    cmap = plt.get_cmap('viridis')
    norm = mpl.colors.Normalize(vmin=min(times), vmax=max(times))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=cax_porosity, ticks=[min(times), max(times)], orientation='horizontal')
    cax_porosity.xaxis.set_ticks_position('top')
    cax_porosity.xaxis.set_label_position('top')
    cbar.set_label('$t$', labelpad=-3)


    ax.set_xlim([0, 1])

    # Plot error
    difference = T_data - T_python
    ax_err = ax.twinx()
    plt_err = ax_err.plot(z, difference, color='black', linestyle=':', label='Error (right axis)')
    # ax_err.set_yscale('symlog', linthreshx=1e-10)

    ax_err.set_ylabel('Error')

    # Add legend
    lns = plt_python + plt_err
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)

    # Sort axis positioning
    ax.set_zorder(1)  # make it on top
    ax.set_frame_on(False)  # make it transparent
    ax_err.set_frame_on(True)  # make sure there is any background

    # Sort figure positioning
    ax.set_position([0.13, 0.15, 0.65, 0.7])
    ax_err.set_position([0.13, 0.15, 0.65, 0.7])

    # Save figure
    for figure_output_path in figure_output_paths:
        fig_name = figure_output_path + '.eps'
        print('Saving figure to %s' % fig_name)
        plt.savefig(fig_name, format='eps')

    plt.show()

