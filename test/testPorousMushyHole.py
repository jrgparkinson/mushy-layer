# Run all convergence tests for the methods paper
import os, sys
from colorama import Fore, Style
import getopt
import math

from runAMRConvergenceTest import runTest
from mushyLayerRunUtils import get_base_output_dir, get_matlab_base_command, read_inputs

##########################################################################
# 3) Convection in a mushy layer with an initial porous hole
##########################################################################


def porous_mushy_hole_resolution_specific_params(nz_coarse, ref_rat, max_level, max_refinement):
    mushyLayerBaseDir = os.environ['MUSHY_LAYER_DIR']
    params_file = os.path.join(mushyLayerBaseDir,'params/convergenceTest/porousMushyHole.parameters')

    params = read_inputs(params_file)


    nx_coarse = nz_coarse

    Nx_grid = nx_coarse
    if ref_rat == 4:
        Nx_grid = 2 * nx_coarse

    gridFile = mushyLayerBaseDir + '/grids/middleXSmall/' + str(Nx_grid) + 'x' + str(Nx_grid)

    # dx = 1 / float(nx_coarse)
    # dt = 5e-7 * 64.0 / float(nx_coarse)
    # Coarsest Nx is 8
    coarsest_dt = 1e-4
    coarsest_nx = 8
    dt = float(coarsest_dt) * float(coarsest_nx) / float(nx_coarse)
    params['main.fixed_dt'] = dt  # dt
    params['main.plot_period'] = float(params['main.max_time'])/4  # produce 4 same number of plot files per run

    # Check this dt will go into max_time an integer number of times
    # If it doesn't, different simulations will be integrated for different lengths of time which can be an issue
    # when we come to comparing them (they won't be the same)
    num_steps = float(params['main.max_time'])/dt
    if not num_steps.is_integer():
        print(Fore.RED + 'WARNING: dt is not an integer division of max time')

    regrid_int = int(math.ceil(float(nz_coarse) / 16.0))
    regrid_int = max(regrid_int, 1)

    params['main.refine_thresh'] = float(params['main.radius']) * 1.0 / float(nz_coarse)  # 0.05*64.0/float(nz_coarse)

    bf = max(int(nx_coarse / 4), 4)

    maxGridSize = 1024 #bf * 2
    gbs = max(bf / 8, 2)

    params['main.block_factor'] = str(bf)
    params['main.grid_buffer_size'] = gbs
    params['main.tag_buffer_size'] = 0
    params['main.regrid_interval'] = str(regrid_int) + ' ' + str(regrid_int) + ' ' + str(regrid_int)
    params['main.max_grid_size'] = str(maxGridSize)  # Must be greater than block factor

    return nx_coarse, params, gridFile

def testPorousMushyHole(argv):

    # Defaults:
    max_time = 5.0e-4
    hole_radius = 0.04

    try:
        opts, args = getopt.getopt(argv, "t:h:")
    except getopt.GetoptError as err:
        print(str(err))
        print('testPorousMushyHole.py -t<max time> -h<hole radius=0.02>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in "-t":
            max_time = float(arg)
        elif opt in "-h":
            hole_radius = float(arg)

    base_output_dir = get_base_output_dir()
    matlab_command = get_matlab_base_command()

    print(Fore.GREEN + 'Setup tests for porous mushy hole' + Style.RESET_ALL)
    physicalProblem = 'PorousMushyHole'
    dataFolder = os.path.join(base_output_dir, 'PorousMushyHole-t' + str(max_time) + '-hole' + str(hole_radius) )

    Nz_uniform = [16, 32, 64, 128, 256, 512]
    AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': Nz_uniform},
                {'max_level': 1, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [16, 32, 64, 128]},
                {'max_level': 2, 'ref_rat': 2, 'run_types': ['amr'], 'Nzs': [8, 16, 32, 64]},
                {'max_level': 1, 'ref_rat': 4, 'run_types': ['amr'], 'Nzs': [8, 16, 32, 64]}]

    # While testing:
   # Nz_uniform = [16, 32, 64]
   # AMRSetup = [{'max_level': 0, 'ref_rat': 1, 'run_types': ['uniform'], 'Nzs': Nz_uniform}]

    # Nzs 	  = [16, 32, 64]
    #num_procs = [1, 1, 1, 4, 4, 4]  # Needs to be as long as the longest Nzs
    num_procs = [1] * len(Nz_uniform)

    # Setup up the post processing command

    # figureName = os.path.join(dataFolder, 'noFlow.pdf')
    uniform_prefix = 'Uniform-' + physicalProblem + '-'

    python_compare_file = os.path.join(os.environ['MUSHY_LAYER_DIR'], 'test', 'run_chombo_compare.py')
    chombo_compare_analyse ='python %s -f %s -a -v Porosity -e L2 -r True \n \n' % (python_compare_file, dataFolder)

    analysis_command = chombo_compare_analyse + '\n\n' + \
                       matlab_command + ' "analyseVariablePorosityTest(\'' + dataFolder + '\', [' + ','.join([str(a) for a in Nz_uniform]) + '],' \
     'true, true, \'' + uniform_prefix + '\', \'Porosity\', \'L2\'); exit;"'

    # Run
    extra_params = {'main.max_time':max_time,
                    'main.radius': hole_radius}

    runTest(dataFolder, physicalProblem, porous_mushy_hole_resolution_specific_params,
            AMRSetup, num_procs, analysis_command, extra_params)

if __name__ == "__main__":
    testPorousMushyHole(sys.argv[1:])



