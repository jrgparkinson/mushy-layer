# Script to test convergence of uniform and AMR solutions
import os
import sys
import math
import getopt
from colorama import Fore, Style
from mushyLayerRunUtils import read_inputs, is_power_of_two, string_to_array
from BatchJob import BatchJob
from AMRConvergenceTest import amr_convergence_test


#TODO - should allow this function to set which servers to run on

def runTest(base_dir, physical_problem, resolution_specific_params, AMRSetup, num_procs,
            analysis_command='', extra_params={},
            numRestarts=0, params_file='', restart_from_low_res=False):

    # base_dir should be e.g.
    # '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/Test/AMRConvergenceTestNoFlow'

    mushyLayerBaseDir = os.path.abspath(os.pardir)

    all_params = []

    job_ids = []

    for setup in AMRSetup:
        max_level = setup['max_level']
        ref_rat = max(setup['ref_rat'], 1)
        Nzs = setup['Nzs']

        runTypes = ['uniform']
        if 'run_types' in setup:
            runTypes = setup['run_types']

        # finestRefinement = pow(ref_rat, max_level)
        maxRefinement = ref_rat ** max_level
        # output_dir = ''
        Nz_i = -1

        # Construct the params files
        for nz_coarse in Nzs:

            Nz_i = Nz_i + 1

            # Some default options

            # Use same aspect ratio as already defined
            # nx_coarse = -1  # if this isn't changed, we'll eventually just use the predetermined aspect ratio
            # gridFile = ''

            # defaultParamsFile = os.path.join(mushyLayerBaseDir, '/params/convergenceTest/' + physical_problem + '.parameters')
            # if os.path.exists(defaultParamsFile):
            #    params = read_inputs(defaultParamsFile)

            output_dir = ''

            nx_coarse, params, gridFile = resolution_specific_params(nz_coarse, ref_rat, max_level, maxRefinement)



            # Default options
            if nx_coarse == -1:
                print('Trying to convert to an array: ' + str(params['main.num_cells']))
                gridPts = string_to_array(params['main.num_cells'])

                aspectRatio = gridPts[0] / gridPts[1]
                nx_coarse = aspectRatio * nz_coarse

            if not gridFile:
                gridFile = mushyLayerBaseDir + '/grids/middle/' + str(nz_coarse) + 'x' + str(nz_coarse)

            # Always turn off slope limiting for convergence tests
            params['main.use_limiting'] = 'false'
            params['main.debug'] = 'false'  # also turn off debug to save disk space

            # Any extra params we may have
            for k, v in extra_params.iteritems():
                params[k] = v

            #numCellsAMR = str(nx_coarse) + ' ' + str(nz_coarse) + '  8'    
            numCellsAMR = [int(nx_coarse), int(nz_coarse), 8]

            gridFileQuarter = mushyLayerBaseDir + '/grids/rightQuarter/' + str(nx_coarse) + 'x' + str(nz_coarse)
            gridFileThreeLevels = mushyLayerBaseDir + '/grids/rightHalfThreeLevel/' + str(nx_coarse) + 'x' + str(nz_coarse)
            # numCellsUniform = str(nx_coarse * finestRefinement) + ' ' + str(nz_coarse * finestRefinement) + '  8'
            # numCellsUniform = str(nx_coarse) + ' ' + str(nz_coarse) + '  8'
            numCellsUniform = [int(nx_coarse), int(nz_coarse), 8]

            # params['main.gridfile'] = gridFile
            params['main.num_cells'] = numCellsAMR
            params['main.max_level'] = str(max_level)
            params['main.ref_ratio'] = str(ref_rat) + ' ' + str(ref_rat) + ' ' + str(ref_rat) 
            
            num_proc = num_procs[Nz_i]
            params['num_proc'] = num_proc

            if num_proc > 1:
                optimalGrid = float(nx_coarse) / (float(num_proc) / 2.0)

                # print('Initial optimal grid guess: %d' % optimalGrid)
                
                # increase to next power of 2
                while not is_power_of_two(optimalGrid):
                    optimalGrid = optimalGrid + 1

                # print('Final optimal grid guess: %d' % optimalGrid)
                
                maxGrid = max(16, optimalGrid)
                # grid size must be greater than the blocking factor 
                maxGrid = max(maxGrid, 2 * int(params['main.block_factor'])) 
                params['main.max_grid_size'] = str(int(maxGrid))

            # Now construct vector different param sets
            param_sets = []
            
            # Run 1: uniform mesh
            if 'uniform' in runTypes:
                p0 = dict(params)
                p0['main.max_level'] = '0'
                p0['main.num_cells'] = numCellsUniform
                p0['run_name'] = 'Uniform'
                p0['concise_run_name'] = 'Uniform'
                p0['projection.eta'] = 0.0
            
                param_sets.append(p0)
            
            # This should do the job for AMR simulations
            if 'amr' in runTypes:
                p1 = dict(params)
                p1['main.reflux_scalar'] = '1'
                p1['main.use_subcycling'] = '1'
                p1['main.refluxType'] = '2'
                p1['projection.eta'] = '0.99'
                
                p1['run_name'] = 'AMR-Subcycle-Reflux-Freestream' + str(p1['projection.eta']) + '-MaxLevel' + str(max_level) 
                p1['concise_run_name'] = 'AMR'
                # p13['main.initialize_VD_corr'] = 'true' # true by default
            
                param_sets.append(p1)
            
            # Variable mesh
            if 'variable' in runTypes and max_level == 1:    
                p2 = dict(params)
                p2['main.reflux_scalar'] = '1'
                p2['main.use_subcycling'] = '1'
                p2['main.refluxType'] = '2'
                p2['projection.eta'] = '0.95'

                p2['main.gridfile'] = gridFile
                p2['run_name'] = 'VM-Subcycle-Reflux-Freestream' + str(p2['projection.eta']) + '-MaxLevel' + str(max_level)
                p2['concise_run_name'] = 'VM'
            
                param_sets.append(p2)
            
            if 'variable' in runTypes and max_level == 2: 
                p3 = dict(params)
                p3['main.reflux_scalar'] = '1'
                p3['main.use_subcycling'] = '1'
                p3['main.refluxType'] = '2'
                p3['projection.eta'] = '0.95'

                p3['main.gridfile'] = gridFileThreeLevels
                p3['run_name'] = 'VM-Subcycle-Reflux-Freestream' + str(p3['projection.eta']) + '-MaxLevel' + str(max_level)
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
        #full_output_dir = os.path.join(base_dir, output_dir)
        full_output_dir = base_dir
        
        these_job_ids = amr_convergence_test(all_params, full_output_dir, 
                                             physical_problem, Nzs, num_procs, 
                                             numRestarts, analysis_command,
                                             restart_from_low_res)

        # Concatenate lists
        job_ids = job_ids + these_job_ids

    # Once all these runs have been submitted, submit the analysis job
    if analysis_command:
        runAnalysisName = 'runAnalysis.sh'

        # Don't redo analysis - we may be waiting on runs to finish
        if os.path.exists(os.path.join(base_dir, runAnalysisName)):
            print(Fore.YELLOW + 'Analysis job already submitted \n' + Fore.RESET)
        else:
            jobName = physical_problem + '-analysis'

            s = BatchJob(base_dir, jobName, '', 4)

            s.set_dependency(job_ids)
            s.set_custom_command(analysis_command)

            # s.write_slurm_file(runAnalysisName)
            s.run_task(runAnalysisName)
            print(Fore.GREEN + 'Submitted analysis job \n' + Fore.RESET)

    return job_ids


def main(argv):
    
    # Defaults:
    num_proc = 4
    
    try:
       opts, args = getopt.getopt(argv, "n:f:H:")
    except getopt.GetoptError:
       print('runAMRConvergenceTest.py -n <num processors>')
       sys.exit(2)
    for opt, arg in opts:
        if opt in ("-n"):
            num_proc = int(arg)
        
    print('Num processors: ' + str(num_proc))

    #################################
    # These shouldn't change
    #################################
    os.environ["CH_TIMER"] = "1"

    base_dir = os.getcwd()  # for local runs
    base_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/'  # for writing to shared storage

    # Change things below here
    # num_proc = 1
    Nzs = [16, 32, 64, 128, 256, 512]
    # Nzs = [128]

    # num_procs = num_proc+0*Nzs
    num_procs = [1, 1, 2, 2, 4, 4, 8, 8, 8, 8, 8, 8]
    # num_procs = [1,1,4,4,4,4,4,4,4]
    # num_procs = [1,1,1,1,1,1,1,1,1,1,1]

    # num_procs = [16]

    # Nzs = [512,256,128,64,32] # more efficient to do it this way around
    
    # By default, don't try and restart from existing files
    # numRestarts = 0

    # max_level = 0
    # ref_rat = 2

    AMRSetup = [{'max_level': 1, 'ref_rat': 2},
    {'max_level': 1, 'ref_rat': 4},
    {'max_level': 2, 'ref_rat': 2}]

    # AMRSetup = [{'max_level': 1, 'ref_rat': 2}]

    # AMRSetup = [{'max_level': 0, 'ref_rat': 2}]

    # mushQuick, MushyDB, DBVariablePorosity, DirectionalDarcy, convectionDarcyBrinkman, DirectionalDB
    # MushyConvectionAnalyticVel, MushyConvection, MushyConvectionLiquid
    physicalProblem = 'MushyConvectionLiquid2' 

    base_dir = os.path.join(base_dir, 'AMRConvergenceTestMushyConvection')

    runTest(base_dir, physicalProblem, AMRSetup, num_procs)


if __name__ == "__main__":
    print('Running script')
    main(sys.argv[1:])

