# Run all convergence tests for the methods paper
import os

from colorama import Fore, Style

from AMRConvergenceTest import run_test
from mushyLayerRunUtils import get_base_output_dir, get_matlab_base_command, read_inputs, get_mushy_layer_dir


######################################
# 1) Diffusive solidification problem
#######################################

def diffusive_solidification_resolution_specific_params(p):
    """
    :param p: ConvergenceTestParams
    :return:
    """
    mushy_layer_base_dir = get_mushy_layer_dir()
    params_file = os.path.join(mushy_layer_base_dir,'params/convergenceTest/noFlowConvTest.parameters')

    params = read_inputs(params_file)
    grid_file = ''

    nx_coarse = 4
    params['main.domain_width'] = str(4.0 * float(nx_coarse) / float(p.nz_coarse))
    linear_temperature_gradient = 1.0 / 4.0  # delta T / height
    params['main.refine_thresh'] = str(linear_temperature_gradient * 16.0 / float(p.nz_coarse))
    params['main.tag_buffer_size'] = str(int(float(p.nz_coarse) / 8))
    params['main.steady_state'] = '1e-8'  # str(1e-4 * pow(32.0/float(nz_coarse),2))
    params['main.max_grid_size'] = '32'
    params['main.max_step'] = 10000
    params['main.debug'] = 'true'  # make sure we output things like Temperature error


    return nx_coarse, params, grid_file


def test_diffusive_solidification():

    base_output_dir = get_base_output_dir()
    matlab_command = get_matlab_base_command()

    print(Fore.GREEN + 'Setup tests for solidification without flow' + Style.RESET_ALL)

    physical_problem = 'noFlow'
    amr_setup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': [8, 16, 32, 64, 128, 256]},
                {'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [8, 16, 32, 64, 128]},
                {'max_level': 2, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [8, 16, 32, 64]}]
    # While testing:
    # AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform']}];

    # Nzs 	  = [16, 32, 64]
    num_procs = [1] * 6  # Needs to be as long as the longest Nzs

    # Setup up the post processing command
    data_folder = os.path.join(base_output_dir, 'NoFlow')
    figure_name = os.path.join(base_output_dir, 'Fig4NoFlow.pdf')
    analysis_command = matlab_command + ' "Fig4NoFlow(\'' + data_folder + '\', \'' + figure_name + '\'); exit;"'

    # Run
    extra_params = {'main.debug': 'true'}
    run_test(data_folder, physical_problem, diffusive_solidification_resolution_specific_params, amr_setup, num_procs,
             analysis_command, extra_params)


if __name__ == "__main__":
    test_diffusive_solidification()