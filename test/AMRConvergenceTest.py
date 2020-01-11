# This file contains a set of utilities for running a convergence test
import os
from colorama import Fore
from MushyLayerRunSimple import MushyLayerRunSimple
from mushyLayerRunUtils import get_restart_file, get_executable_name, \
    get_final_chk_file, get_mushy_layer_dir, is_power_of_two, string_to_array
from BatchJob import BatchJob


class ConvergenceTestParams:

    def __init__(self, nz_coarse, ref_rat, max_level, max_refinement):
        self.nz_coarse = nz_coarse
        self.ref_rat = ref_rat
        self.max_level = max_level
        self.max_refinement = max_refinement


def get_suffix(num):
    # return '-%d' % num
    return ''


def amr_convergence_test(params, full_output_dir, nzs, num_procs=1, num_restarts=0, restart_from_low_res=False):
    os.environ["CH_TIMER"] = "1"

    if not isinstance(num_procs, list):
        num_procs = [num_procs]

    dependencies = []
    prev_job_id = -1
    prev_dir = ''
    prev_num_cells = 0

    while len(num_procs) < len(nzs):
        num_procs.append(num_procs[-1])

    param_i = -1

    for p in params:
        param_i = param_i + 1

        num_proc = p['num_proc']

        # Don't repeat runs, unless we're restarting

        # cwd = os.getcwd()
        run_name = p['main.plot_prefix']
        test_name = run_name + get_suffix(0)
        full_path = os.path.join(full_output_dir, test_name)

        allow_restarts = False

        print(full_path)
        if os.path.exists(full_path):

            print(Fore.YELLOW + '    Run already done' + Fore.RESET)

            if num_restarts > 0:

                print('    ' + str(num_restarts) + ' restart(s) allowed \n ')
                # Find most recent folder and restart from there

                i = 0
                while os.path.exists(full_path):
                    i = i + 1
                    test_name = run_name + get_suffix(i)
                    full_path = os.path.join(full_output_dir, test_name)

                if i > num_restarts:
                    continue

                most_recent_path = os.path.join(full_output_dir, run_name + get_suffix(i - 1))
                # Now get the restart file
                restart_file = get_restart_file(most_recent_path)
                if restart_file:
                    p['main.restart_file'] = os.path.join(most_recent_path, restart_file)
                    allow_restarts = True

                    print('    Set restart file: ' + p['main.restart_file'])
            else:

                continue

        else:
            print(Fore.YELLOW + '    No previous runs' + Fore.RESET)

        print(Fore.GREEN + '    **Do run**' + Fore.RESET)

        # Don't need output dir or exec dir here, MushyLayerRun will fill these in
        s = BatchJob('', p['concise_run_name'], '', num_proc)

        if restart_from_low_res:
            python_file = os.path.join(get_mushy_layer_dir(), 'test', 'create_refined_restart.py')
            new_dir = full_path

            if prev_job_id > -1:

                refinement = int(float(p['main.num_cells'][0]) / float(prev_num_cells[0]))

                preprocess_cmd = 'python %s -p %s -n %s -r %d \n' % (python_file, prev_dir, new_dir, refinement)
                s.set_dependency(prev_job_id)
                s.set_preprocess(preprocess_cmd)
                p['main.restart_file'] = os.path.join(new_dir, 'restart.2d.hdf5')
            else:
                # Search for an old file we can restart from, just in case
                refinement = 2
                this_nx = int(p['main.num_cells'][0])
                coarse_nx = int(float(this_nx) / 2.0)
                coarser_dir = new_dir.replace(str(this_nx), str(coarse_nx))

                if os.path.exists(coarser_dir) and get_final_chk_file(coarser_dir):
                    preprocess_cmd = 'python %s -p %s -n %s -r %d \n' % (python_file, coarser_dir, new_dir, refinement)
                    s.set_preprocess(preprocess_cmd)
                    p['main.restart_file'] = os.path.join(new_dir, 'restart.2d.hdf5')

        ml_run = MushyLayerRunSimple(full_output_dir, num_proc, p, s, allow_restarts, get_executable_name())
        ml_run.single_run(run_name)

        prev_job_id = s.job_id
        prev_dir = full_path
        prev_num_cells = p['main.num_cells']
        dependencies.append(s.job_id)

        print('==================')

    return dependencies


def run_test(base_dir, physical_problem, resolution_specific_params, amr_setup, num_procs, analysis_command='',
             extra_params=None, num_restarts=0, restart_from_low_res=False):
    """ Driver for the AMRConvergenceTest class """

    # base_dir should be e.g.
    # '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/Test/AMRConvergenceTestNoFlow'

    mushy_layer_base_dir = os.path.abspath(os.pardir)

    all_params = []

    job_ids = []

    for setup in amr_setup:
        max_level = setup['max_level']
        ref_rat = max(setup['ref_rat'], 1)
        nzs = setup['Nzs']

        run_types = ['uniform']
        if 'run_types' in setup:
            run_types = setup['run_types']

        # finestRefinement = pow(ref_rat, max_level)
        max_refinement = ref_rat ** max_level
        # output_dir = ''
        nz_i = -1

        # Construct the params files
        for nz_coarse in nzs:

            nz_i = nz_i + 1

            res_params = ConvergenceTestParams(nz_coarse, ref_rat, max_level, max_refinement)
            nx_coarse, params, grid_file = resolution_specific_params(res_params)

            # Default options
            if nx_coarse == -1:
                print('Trying to convert to an array: ' + str(params['main.num_cells']))
                grid_pts = string_to_array(params['main.num_cells'])

                aspect_ratio = grid_pts[0] / grid_pts[1]
                nx_coarse = aspect_ratio * nz_coarse

            if not grid_file:
                grid_file = mushy_layer_base_dir + '/grids/middle/' + str(nz_coarse) + 'x' + str(nz_coarse)

            # Turn slope limiting off unless we've requested it
            if 'main.use_limiting' not in params:
                params['main.use_limiting'] = 'false'

            # also turn off debug to save disk space, unless we've explicitly asked for it
            if 'main.debug' not in params:
                params['main.debug'] = 'false'

            # Any extra params we may have
            if extra_params:
                for k, v in extra_params.items():
                    params[k] = v

            # numCellsAMR = str(nx_coarse) + ' ' + str(nz_coarse) + '  8'
            num_cells_amr = [int(nx_coarse), int(nz_coarse), 8]

            # gridFileQuarter = mushyLayerBaseDir + '/grids/rightQuarter/' + str(nx_coarse) + 'x' + str(nz_coarse)
            grid_file_three_levels = mushy_layer_base_dir + '/grids/rightHalfThreeLevel/' + str(nx_coarse) + 'x' + str(
                nz_coarse)
            # numCellsUniform = str(nx_coarse * finestRefinement) + ' ' + str(nz_coarse * finestRefinement) + '  8'
            # numCellsUniform = str(nx_coarse) + ' ' + str(nz_coarse) + '  8'
            num_cells_uniform = [int(nx_coarse), int(nz_coarse), 8]

            # params['main.gridfile'] = gridFile
            params['main.num_cells'] = num_cells_amr
            params['main.max_level'] = str(max_level)
            params['main.ref_ratio'] = str(ref_rat) + ' ' + str(ref_rat) + ' ' + str(ref_rat)

            num_proc = num_procs[nz_i]
            params['num_proc'] = num_proc

            if num_proc > 1:
                optimal_grid = float(nx_coarse) / (float(num_proc) / 2.0)

                # print('Initial optimal grid guess: %d' % optimalGrid)

                # increase to next power of 2
                while not is_power_of_two(optimal_grid):
                    optimal_grid = optimal_grid + 1

                # print('Final optimal grid guess: %d' % optimalGrid)

                max_grid = max(16, optimal_grid)
                # grid size must be greater than the blocking factor
                max_grid = max(max_grid, 2 * int(params['main.block_factor']))
                params['main.max_grid_size'] = str(int(max_grid))

            # Now construct vector different param sets
            param_sets = []

            # Run 1: uniform mesh
            if 'uniform' in run_types:
                p0 = dict(params)
                p0['main.max_level'] = '0'
                p0['main.num_cells'] = num_cells_uniform
                p0['run_name'] = 'Uniform'
                p0['concise_run_name'] = 'Uniform'
                p0['projection.eta'] = 0.0

                param_sets.append(p0)

            # This should do the job for AMR simulations
            if 'amr' in run_types:
                p1 = dict(params)
                p1['main.reflux_scalar'] = '1'
                p1['main.use_subcycling'] = '1'
                p1['main.refluxType'] = '2'
                p1['projection.eta'] = '0.99'

                p1['run_name'] = 'AMR-Subcycle-Reflux-Freestream' + str(p1['projection.eta']) + '-MaxLevel' + str(
                    max_level)
                p1['concise_run_name'] = 'AMR'
                # p13['main.initialize_VD_corr'] = 'true' # true by default

                param_sets.append(p1)

            # Variable mesh
            if 'variable' in run_types and max_level == 1:
                p2 = dict(params)
                p2['main.reflux_scalar'] = '1'
                p2['main.use_subcycling'] = '1'
                p2['main.refluxType'] = '2'
                p2['projection.eta'] = '0.95'

                p2['main.gridfile'] = grid_file
                p2['run_name'] = 'VM-Subcycle-Reflux-Freestream' + str(p2['projection.eta']) + '-MaxLevel' + str(
                    max_level)
                p2['concise_run_name'] = 'VM'

                param_sets.append(p2)

            if 'variable' in run_types and max_level == 2:
                p3 = dict(params)
                p3['main.reflux_scalar'] = '1'
                p3['main.use_subcycling'] = '1'
                p3['main.refluxType'] = '2'
                p3['projection.eta'] = '0.95'

                p3['main.gridfile'] = grid_file_three_levels
                p3['run_name'] = 'VM-Subcycle-Reflux-Freestream' + str(p3['projection.eta']) + '-MaxLevel' + str(
                    max_level)
                p3['concise_run_name'] = 'VM'

                param_sets.append(p3)

            # Do this for all the param sets just made:
            for p in param_sets:

                # if num_proc > 1:
                #    p['run_name'] = 'p' + p['run_name']

                run_name = p['run_name']
                run_name = run_name
                if int(p['main.max_level']) > 0:
                    run_name = run_name + '-ref' + str(ref_rat)

                run_name = run_name + '-' + physical_problem + '-' + str(nz_coarse)

                p['concise_run_name'] = p['concise_run_name'] + str(nz_coarse) + physical_problem
                p['main.plot_prefix'] = run_name + '-'
                p['main.chk_prefix'] = 'chk' + p['main.plot_prefix']

                # Now add to the full list of param sets
                all_params.append(p)

        # Actually run the convergence test
        # full_output_dir = os.path.join(base_dir, output_dir)
        full_output_dir = base_dir

        these_job_ids = amr_convergence_test(all_params, full_output_dir, nzs, num_procs, num_restarts,
                                             restart_from_low_res)

        # Concatenate lists
        job_ids = job_ids + these_job_ids

    # Once all these runs have been submitted, submit the analysis job
    print('analysis command: %s' % analysis_command)
    if analysis_command:
        run_analysis_name = 'runAnalysis.sh'

        # Don't redo analysis - we may be waiting on runs to finish
        if os.path.exists(os.path.join(base_dir, run_analysis_name)):
            print(Fore.YELLOW + 'Analysis job already submitted \n' + Fore.RESET)
        else:
            job_name = physical_problem + '-analysis'

            s = BatchJob(base_dir, job_name, '', 4)

            # If queuing system can't be found, try and run commands manually (will be slow)
            s.allow_manual_run = True

            s.set_dependency(job_ids)
            s.set_custom_command(analysis_command)

            s.run_task(run_analysis_name)
            print(Fore.GREEN + 'Submitted analysis job \n' + Fore.RESET)

    return job_ids
