import getopt
import math
import os
import sys
from colorama import Fore, Style
from AMRConvergenceTest import run_test
from makeFigures import porous_hole_command
from mushyLayerRunUtils import get_base_output_dir, get_matlab_base_command, read_inputs, get_mushy_layer_dir
from testFixedPorousHole import get_setup

##########################################################################
# 3) Convection in a mushy layer with an initial porous hole
##########################################################################


def porous_mushy_hole_resolution_specific_params(p):
    mushy_layer_base_dir = get_mushy_layer_dir()
    params_file = os.path.join(mushy_layer_base_dir, 'params/convergenceTest/porousMushyHole.parameters')

    params = read_inputs(params_file)

    nx_coarse = p.nz_coarse

    nx_grid = nx_coarse
    if p.ref_rat == 4:
        nx_grid = 2 * nx_coarse

    grid_file = mushy_layer_base_dir + '/grids/middleXSmall/' + str(nx_grid) + 'x' + str(nx_grid)

    # dx = 1 / float(nx_coarse)
    # dt = 5e-7 * 64.0 / float(nx_coarse)
    # Coarsest Nx is 8
    coarsest_dt = 1e-4
    coarsest_nx = 8
    dt = float(coarsest_dt) * float(coarsest_nx) / float(nx_coarse)
    params['main.fixed_dt'] = dt  # dt
    params['main.plot_period'] = float(params['main.max_time']) / 4  # produce same number of plot files per run (4)

    # Check this dt will go into max_time an integer number of times
    # If it doesn't, different simulations will be integrated for different lengths of time which can be an issue
    # when we come to comparing them (they won't be the same)
    num_steps = float(params['main.max_time']) / dt
    if not num_steps.is_integer():
        print(Fore.RED + 'WARNING: dt is not an integer division of max time')

    regrid_int = int(math.ceil(float(p.nz_coarse) / 16.0))
    regrid_int = max(regrid_int, 1)

    params['main.refine_thresh'] = float(params['main.radius']) * 10.0 / float(p.nz_coarse)

    bf = max(int(nx_coarse / 4), 4)

    max_grid_size = 1024  # bf * 2
    gbs = max(int(bf / 8), 2)

    params['main.block_factor'] = str(bf)
    params['main.grid_buffer_size'] = gbs
    params['main.tag_buffer_size'] = 0
    params['main.regrid_interval'] = str(regrid_int) + ' ' + str(regrid_int) + ' ' + str(regrid_int)
    params['main.max_grid_size'] = str(max_grid_size)  # Must be greater than block factor

    return nx_coarse, params, grid_file


def test_porous_mushy_hole(argv):
    # Defaults:
    max_time = 1.5e-4
    hole_radius = 0.04
    plot_interval = None

    try:
        opts, _ = getopt.getopt(argv, "t:h:p:")
    except getopt.GetoptError as err:
        print(str(err))
        print('testPorousMushyHole.py -t<max time> -h<hole radius=0.02>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in "-t":
            max_time = float(arg)
        elif opt in "-h":
            hole_radius = float(arg)
        elif opt in "-p":
            plot_interval = int(arg)

    base_output_dir = get_base_output_dir()
    matlab_command = get_matlab_base_command()

    print(Fore.GREEN + 'Setup tests for porous mushy hole' + Style.RESET_ALL)
    physical_problem = 'PorousMushyHole'
    min_concentration = -5.0  # default -5
    folder_name = 'PorousMushyHole-t%s-minC%s' % (max_time, min_concentration)
    data_folder = os.path.join(base_output_dir, folder_name)

    amr_setup, num_procs = get_setup()

    # Setup up the post processing command
    python_compare_file = os.path.join(get_mushy_layer_dir(), 'test', 'run_chombo_compare.py')
    chombo_compare_analyse = 'python %s -f %s -a -v Porosity -e L2 -r True -n 8 -d %s \n \n' % (python_compare_file,
                                                                                                data_folder,
                                                                                                base_output_dir)

    make_figures = porous_hole_command(folder_name=folder_name)
    analysis_command = '%s\n\n%s" %s; exit;"' % (chombo_compare_analyse, matlab_command, make_figures)

    # Run
    extra_params = {'main.max_time': max_time,
                    'main.radius': hole_radius,
                    'bc.bulkConcentrationHiVal': [-1, min_concentration]}

    if plot_interval:
        extra_params['main.plot_interval'] = plot_interval

    run_test(data_folder, physical_problem, porous_mushy_hole_resolution_specific_params, amr_setup, num_procs,
             analysis_command, extra_params)


if __name__ == "__main__":
    test_porous_mushy_hole(sys.argv[1:])
