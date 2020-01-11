# Run all convergence tests for the methods paper
import getopt
import os
import sys

from colorama import Fore, Style

from AMRConvergenceTest import run_test
from BatchJob import BatchJob
from MushyLayerRunSimple import MushyLayerRunSimple
from makeFigures import fixed_chill_command
from mushyLayerRunUtils import get_base_output_dir, get_matlab_base_command, read_inputs, get_executable_name, \
    get_mushy_layer_dir, add_params


##########################################################################
# 3) Convection in a mushy layer with an initial porous hole
##########################################################################

def hele_shaw_resolution_specific_params(p):
    """
    :param p: ConvergenceTestParams
    :return:
    """
    mushy_layer_base_dir = get_mushy_layer_dir()

    params_file = mushy_layer_base_dir + '/params/convergenceTest/FixedChill.parameters'

    params = read_inputs(params_file)
    nx_coarse = -1

    grid_file = ''  # mushyLayerBaseDir + '/grids/topMiddle/' + str(nx_coarse) + 'x' + str(nz_coarse)

    return nx_coarse, params, grid_file


def test_hele_shaw_fixed_chill(argv):
    mushy_layer_base_dir = os.path.abspath(os.pardir)
    # params_file = mushyLayerBaseDir + '/params/convergenceTest/FixedChill.parameters'
    params_file = mushy_layer_base_dir + '/params/convergenceTest/FixedChill.parameters'

    # Defaults
    default_params = read_inputs(params_file)
    extra_params = {'main.max_time': float(default_params['main.max_time']),
                    'parameters.rayleighComp': float(default_params['parameters.rayleighComp']),
                    'parameters.darcy': float(default_params['parameters.darcy']),
                    'parameters.compositionRatio': float(default_params['parameters.compositionRatio']),
                    'parameters.nonDimReluctance': float(default_params['parameters.nonDimReluctance']),
                    'parameters.prandtl': float(default_params['parameters.prandtl'])}

    do_amr = True

    # Pr = 10.0  # fix this for now
    periodic = False

    try:

        opts, _ = getopt.getopt(argv, "t:R:D:C:P:N:A")
    except getopt.GetoptError as err:
        print(str(err))
        print(
            'testHeleShawFixedChill.py -t<max time> -R<Compositional Rayleigh Number> -D<Darcy number> '
            '-C<Composition Ratio> -A<do amr> -P<true(1)/false(0) do plot files> -N<non dimensional reluctance>')

        sys.exit(2)

    for opt, arg in opts:
        if opt in "-t":
            extra_params['main.max_time'] = float(arg)
        elif opt in "-R":
            extra_params['parameters.rayleighComp'] = float(arg)
        elif opt in "-D":
            extra_params['parameters.darcy'] = float(arg)
        elif opt in "-C":
            extra_params['parameters.compositionRatio'] = float(arg)
        elif opt in "-N":
            extra_params['parameters.nonDimReluctance'] = float(arg)
        elif opt in "-P":
            do_plot_files = bool(int(arg))
            print('Do plot files: ' + str(do_plot_files) + ', arg = ' + str(arg))
            if not do_plot_files:
                extra_params['main.plot_interval'] = -1
        elif opt in "-A":
            do_amr = True

    base_output_dir = get_base_output_dir()

    print(Fore.GREEN + 'Setup tests for fixed chill in a Hele-Shaw cell' + Style.RESET_ALL)
    physical_problem = 'FixedChill'
    folder_name = "FixedChill-t%1.1e-Ra%.0e-Da%1.1e-C%1.2f-Rel%1.1e" % (extra_params['main.max_time'],
                                                                        extra_params['parameters.rayleighComp'],
                                                                        extra_params['parameters.darcy'],
                                                                        extra_params['parameters.compositionRatio'],
                                                                        extra_params['parameters.nonDimReluctance'])
    if periodic:
        folder_name = folder_name + '-periodic'
    data_folder = os.path.join(base_output_dir, folder_name)

    analysis_command = get_matlab_base_command() + ' "%s; exit;"' % fixed_chill_command(folder_name + '-0')

    single_run = True

    if single_run:

        num_proc = 1
        default_params['concise_run_name'] = folder_name

        default_params = add_params(default_params, extra_params)

        # Overwrite default params with any extra params we've specified
        default_params = add_params(default_params, extra_params)

        allow_restarts = False
        s = BatchJob('', default_params['concise_run_name'], '', num_proc)
        s.set_post_process(analysis_command)
        ml_run = MushyLayerRunSimple(base_output_dir, num_proc, default_params, s, allow_restarts,
                                     get_executable_name())
        ml_run.single_run(folder_name)

    else:

        # Nz_uniform = 256
        # Nz_amr_2 = int(float(Nz_uniform) / 2)
        # Nz_amr_4 = int(float(Nz_uniform) / 4)
        # amr_setup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': [Nz_uniform]},
        #             {'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [Nz_amr_2]},
        #             {'max_level': 2, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [Nz_amr_4]},
        #             {'max_level': 1, 'ref_rat': 4, 'run_types': ['amr'], 'Nzs': [Nz_amr_4]}]

        # While testing:
        nz_uniform = [64]
        # amr_setup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': Nz_uniform},
        #            {'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [64]}]
        amr_setup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': nz_uniform}]

        if do_amr:
            amr_setup.append({'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': nz_uniform})
            amr_setup.append({'max_level': 2, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': nz_uniform})

        num_procs = [1]  # Needs to be as long as the longest Nzs

        run_test(data_folder, physical_problem, hele_shaw_resolution_specific_params, amr_setup, num_procs,
                 analysis_command, extra_params, 0)


if __name__ == "__main__":
    test_hele_shaw_fixed_chill(sys.argv[1:])
