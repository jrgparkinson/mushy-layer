# Run all convergence tests for the methods paper
import os, sys
from colorama import Fore, Style

from runAMRConvergenceTest import runTest
from mushyLayerRunUtils import get_base_output_dir, get_matlab_base_command, read_inputs

##########################################################################
# 3) Convection in a fixed porous medium with variable porosity
##########################################################################

def fixed_porous_resolution_specific_params(nz_coarse, ref_rat, max_level, max_refinement):
    mushy_layer_base_dir = os.environ['MUSHY_LAYER_DIR']

    integration_time = 2e-4
    nx_coarse = nz_coarse

    if ref_rat == 2:
        grid_file = mushy_layer_base_dir + '/grids/middleXSmall/' + str(nx_coarse) + 'x' + str(nz_coarse)
    else:
        grid_file = mushy_layer_base_dir + '/grids/middleXSmall/' + str(nx_coarse * 2) + 'x' + str(nz_coarse * 2)
    # params = readInputs(mushyLayerBaseDir + )

    params_file = mushy_layer_base_dir + '/params/convergenceTest/DBVariablePorosityConvTest.parameters'

    params = read_inputs(params_file)

    # runTypes = ['uniform', 'amr', 'variable']
    # runTypes = ['uniform', 'amr']
    # runTypes = ['uniform']
    # runTypes = ['amr']

    # For a fixed dt:
    dt = 1e-5 * float(16.0 / float(nz_coarse))
    num_steps = float(integration_time) / dt

    params['main.plot_interval'] = str(int(num_steps / 5.0))

    params['main.max_step'] = '100000'

    params['main.min_time'] = integration_time
    params['main.max_time'] = integration_time

    # params['main.vel_refine_thresh'] = 5.0
    # params['main.stdev'] = '0.002'

    regrid_int = int(4.0 * float(nx_coarse) / 16.0)

    # Make sure we always have one grid
    max_grid_size = max(nx_coarse * max_refinement, 4)

    # bf = max(maxGridSize/2,4)
    bf = max(nx_coarse * max_refinement / 8, 4)

    grid_buffer = 0  # max(bf/4,1)
    tag_buffer = 0  # gridBuffer

    if max_level == 2:
        grid_buffer = max(bf / 4, 4)

    params['main.block_factor'] = str(bf)
    params['main.grid_buffer_size'] = str(grid_buffer)
    params['main.tag_buffer_size'] = str(tag_buffer)
    params['main.regrid_interval'] = str(regrid_int) + ' ' + str(regrid_int) + ' ' + str(regrid_int)
    params['main.max_grid_size'] = str(max_grid_size)  # Must be greater than block factor

    chi = '0.05'
    params['bc.porosityHiVal'] = chi + ' ' + chi  # 0.4 0.4
    params['bc.porosityLoVal'] = chi + ' ' + chi

    params['main.fixed_dt'] = str(dt)

    return nx_coarse, params, grid_file



def test_fixed_porous_hole():

    base_output_dir = get_base_output_dir()
    matlab_command = get_matlab_base_command()

    print(Fore.GREEN + 'Setup tests for (fixed) porous hole' + Style.RESET_ALL)
    physical_problem = 'DBVariablePorosity'
    data_folder = os.path.join(base_output_dir, 'FixedPorousHole-1proc')

    nz_uniform = [16, 32, 64, 128, 256, 512]
    amr_setup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': nz_uniform},
                {'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [16, 32, 64, 128]},
                {'max_level': 2, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [8, 16, 32, 64]},
                {'max_level': 1, 'ref_rat': 4, 'run_types': ['amr'], 'Nzs': [8, 16, 32, 64]}]

    # Nzs 	  = [16, 32, 64]
    num_procs = [1] * len(nz_uniform) #[1, 1, 1, 4, 4, 4]  # Needs to be as long as the longest Nzs

    # Setup up the post processing command

    uniform_prefix = 'Uniform-DBVariablePorosity-'

    python_compare_file = os.path.join(os.environ['MUSHY_LAYER_DIR'], 'test', 'run_chombo_compare.py')
    chombo_compare_analyse ='python %s -f %s -a -v \'xDarcy velocity\' -e L2 -r True -n 6 \n \n' % (python_compare_file, data_folder)

    # make_figures = 'Fig5FixedPorosityConvectionPlots(\'' + data_folder_nu + '\', \'' + data_folder_variable_porosity + '\', \'' + figure_directory + '\')'

    analysis_command = chombo_compare_analyse + '\n\n' + \
                       matlab_command + ' "analyseVariablePorosityTest(\'' + data_folder + '\', [' + ','.join([str(a) for a in nz_uniform]) + '], ' \
    'true, true, \'' + uniform_prefix + '\', \'xDarcyvelocity\', \'L2\'); exit;"'



    # Run
    extra_params = {}
    runTest(data_folder, physical_problem, fixed_porous_resolution_specific_params, amr_setup, num_procs, analysis_command, extra_params)

def main(argv):
    test_fixed_porous_hole()

if __name__ == "__main__":
    main(sys.argv[1:])

