from PltFile import PltFile, latexify2
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_bvp
import os
import mushyLayerRunUtils
import subprocess
import sys
import shutil


# noinspection PyUnresolvedReferences
class PoiseuilleSolution:

    def __init__(self, x_sol, porosity_sol, permeability_sol, da, body_force):
        """
        :param da: darcy number, scalar
        :param body_force: body_force, scalar
        :param permeability_sol: value of the permeability, scalar or vector
        :param porosity_sol: value of the porosity, scalar or vector
        :param x_sol: x values on which to compute the solution
        """

        self.da = da
        self.body_force = body_force
        self.permeability = permeability_sol
        self.porosity = porosity_sol
        self.x = x_sol

    def compute_solution(self):
        """
        :return: vertical velocity profile
        """

        # initial guess for v and v derivatives
        v_init = np.zeros((2, self.x.size))

        # Compute the result
        res = solve_bvp(self.ode, self.bc, self.x, v_init)

        # Get the vertical velocity out of the result
        v_analytic = res.sol(self.x)[0]
        return v_analytic

    def ode(self, x_ode, v_ode):
        """ Define ODE:
            v'' = (porosity/permeability)*(1/da)*v - porosity(x)*body_force
        which becomes
            v0' = v1
            v1'  = (porosity/permeability)*(1/da)v0 - porosity(x)*body_force

            :param x_ode: x positions in the domain
            :param v_ode: current vertical velocity
            :return: expression of ODE
            """
        chi = np.interp(x_ode, self.x, self.porosity)
        perm = np.interp(x_ode, self.x, self.permeability)
        return np.vstack((v_ode[1], (chi / perm) * (1 / self.da) * v_ode[0] - chi * self.body_force))

    @staticmethod
    def bc(va, vb):
        """
        Define boundary conditions (no flow at either boundary)
        :param va: vertical velocity on x=0 boundary
        :param vb: vertical velocity on x=L boundary
        :return: residual of the boundary condition, i.e. v-0 = v
        """
        return np.array([va[0], vb[0]])


def base_inputs():
    test_dir = os.path.join(mushyLayerRunUtils.get_mushy_layer_dir(), 'test/poiseuilleFlow')
    default_inputs = mushyLayerRunUtils.read_inputs(os.path.join(test_dir, 'inputs'))

    default_inputs['main.min_time'] = 0.2
    default_inputs['main.max_time'] = 3.0
    default_inputs['main.steady_state'] = 1e-5

    return default_inputs


if __name__ == "__main__":

    # Ensure we can find HDF5 libraries to run mushy layer code
    hdf5_library_path = '/home/parkinsonjl/soft/hdf5-1.8.21p-new/lib/'
    if 'LD_LIBRARY_PATH' in os.environ:
        os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + ':' + hdf5_library_path
    else:
        os.environ['LD_LIBRARY_PATH'] = hdf5_library_path

    # All simulations will be put in their own sub directory within this directory:
    base_dir = os.path.join(mushyLayerRunUtils.get_mushy_layer_dir(), 'test', 'poiseuilleFlow', 'output')
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    # Chose some parameter options:
    porosity_functions = ['Constant', 'Linear', 'Gaussian', 'GaussianSinusodial']
    permeability_functions = ['Pure fluid', 'Cubic', 'Kozeny', 'Log', 'Xsquare',
                              'Porosity'  # set permeability = porosity
                              ]

    opts = {'da': 0.01,  # darcy number
            'F': 100,  # magnitude of body force which drives flow
            'porosity_function': 2,  # see porosity_functions above
            'permeability_function': 1  # see permeability functions above
            }

    # Now let's run the mushy-layer simulation
    # First create files and directories
    sim_prefix = 'Poiseuille'
    sim_name = '%s_da%s_F%s_%sporosity_%spermeability' % (sim_prefix, opts['da'], opts['F'],
                                                          porosity_functions[opts['porosity_function']],
                                                          permeability_functions[opts['permeability_function']])

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

    # Create inputs file for simulation by taking the standard inputs and adding some extra options
    # which are specific to this run
    inputs = base_inputs()
    inputs['parameters.darcy'] = opts['da']
    inputs['parameters.body_force'] = opts['F']
    inputs['main.plot_prefix'] = sim_prefix
    inputs['main.porosity_function'] = opts['porosity_function']
    inputs['parameters.permeabilityFunction'] = opts['permeability_function']
    inputs['main.num_cells'] = [256, 32]
    inputs['main.max_grid_size'] = max(inputs['main.num_cells'])
    inputs['main.output_folder'] = full_output_folder

    new_inputs_loc = os.path.join(full_output_folder, 'inputs')
    mushyLayerRunUtils.write_inputs(new_inputs_loc, inputs)

    # Send the command to run the simulation
    exec_loc = mushyLayerRunUtils.get_executable_name(return_full_path=True)
    cmd = 'cd %s; %s inputs' % (full_output_folder, exec_loc)
    print(cmd)
    print('Running simulation... (usually takes ~ 10 seconds) ')

    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(process.returncode)

    print('Finished running simulation')

    # Now the simulation has finished, plot the results and compare with the semi-analytic solution
    plt_files = [f for f in os.listdir(full_output_folder) if (sim_prefix in f and '.2d.hdf5' in f)]
    if len(plt_files) == 0:
        print('No plot files found in %s with prefix %s' % (full_output_folder, sim_prefix))
        sys.exit(-1)

    plt_files = sorted(plt_files)

    # Create the figure window
    latexify2(7.5, 3.5)
    fig, axes = plt.subplots(1, 2)
    ax_list = axes.flatten()
    ax_left = ax_list[0]

    # Get the computed solution
    computed_solution = PltFile(os.path.join(full_output_folder, plt_files[-1]), load_data=True)
    v = computed_solution.get_level_data('yAdvection velocity').mean('y').squeeze()
    x = v.coords['x']
    porosity = computed_solution.get_level_data('Porosity').mean('y').squeeze()
    permeability = computed_solution.get_level_data('Permeability').mean('y').squeeze()

    plt_computed = ax_left.plot(x, v, label='Computed solution')

    # Compute and plot the analytic solution (scaled with max velocity)
    analytic_solution = PoiseuilleSolution(x, porosity, permeability, opts['da'], opts['F'])
    v_python = analytic_solution.compute_solution()
    plt_python = ax_left.plot(x, v_python, label='Semi-analytic solution', color='red', linestyle='--')

    ax_left.set_xlabel('$x$')
    ax_left.set_ylabel('$v$')

    ax_left.set_xlim([0, 1])

    # Plot error
    difference = abs(v - v_python)
    ax_err = ax_left.twinx()
    plt_err = ax_err.plot(x, difference, color='black', linestyle=':', label='Error (right axis)')
    ax_err.set_ylabel('Error')

    # Add legend
    lns = plt_computed + plt_python + plt_err
    labs = [line.get_label() for line in lns]
    ax_left.legend(lns, labs,
                   loc='lower center', bbox_to_anchor=(0.5, 1.0))  # , loc=0

    # Fix axis positioning
    ax_left.set_zorder(1)  # make it on top
    ax_left.set_frame_on(False)  # make it transparent
    ax_err.set_frame_on(True)  # make sure there is any background

    # Fix figure positioning
    ax_pos_left = [0.08, 0.15, 0.3, 0.55]
    ax_pos_right = [0.6, 0.15, 0.3, 0.55]
    ax_left.set_position(ax_pos_left)
    ax_err.set_position(ax_pos_left)

    # Plot the porosity and permeability on the right axis
    ax = ax_list[1]
    plt_porosity = ax.plot(x, porosity, label='Porosity', color='blue', linestyle='-')

    ax.set_xlabel('$x$')
    ax.set_ylabel(r'Porosity, $\chi$')

    ax_perm = ax.twinx()
    plt_permeability = ax_perm.plot(x, permeability, label='Permeability', color='red', linestyle='-')
    ax_perm.set_ylabel(r'Permeability, $\Pi$')
    ax.set_position(ax_pos_right)

    lns = plt_porosity + plt_permeability
    labs = [line.get_label() for line in lns]
    ax.legend(lns, labs, loc='lower center', bbox_to_anchor=(0.5, 1.0))

    # Save figure
    for figure_output_path in figure_output_paths:
        fig_name = figure_output_path + '.eps'
        print('Saving figure to %s' % fig_name)
        plt.savefig(fig_name, format='eps')

    # Display figure to the screen
    plt.show()
