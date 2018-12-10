# Script to test convergence of uniform and AMR solutions
import os
import sys
import math
import getopt
from colorama import Fore, Style


from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile, isPowerOfTwo, str2arr
from SlurmTask import SlurmTask
from AMRConvergenceTest import AMRConvergenceTest


def runTest(base_dir, physicalProblem, AMRSetup, num_procs, analysis_command = '', extra_params={}, numRestarts=0, params_file = ''):

    # base_dir should be e.g. 
    #'/network/group/aopp/oceans/AW002_PARKINSON_MUSH/Test/AMRConvergenceTestNoFlow'

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
    
        finestRefinement = pow(ref_rat, max_level)
        maxRefinement = ref_rat**max_level

        Nz_i = -1

        # Construct the params files 
        for Nz_coarse in Nzs:

            Nz_i = Nz_i + 1

            #num_proc = pow(Nz_coarse/32, 2)
            #num_proc = max(num_proc, 1)
            #print('Nz coarse = %d, initial num proc = %d \n' % (Nz_coarse, num_proc))
            #num_proc = min(num_proc, 16)
            #print('Final num_proc = %d' % num_proc)

            # Some default options

            # Use same aspect ratio as already defined
            Nx_coarse = -1 # if this isn't changed, we'll eventually just use the predetermined aspect ratio
            gridFile = ''

            defaultParamsFile = os.path.join(mushyLayerBaseDir, '/params/convergenceTest/'+physicalProblem+'.parameters')
            if os.path.exists(defaultParamsFile):
                params = readInputs(defaultParamsFile)

            output_dir = '' #AMRConvergenceTest/' + physicalProblem + '/'
            
            if physicalProblem == 'noFlow':
                #output_dir = 'AMRConvergenceTestNoFlow'
                #params = readInputs(mushyLayerBaseDir + '/params/convergenceTest/.parameters')
                if not params_file:
                    params_file = mushyLayerBaseDir + '/params/convergenceTest/noFlowConvTest.parameters'

                params = readInputs(params_file)
                
                
                Nx_coarse = 4
                params['main.domain_width'] = str(4.0*float(Nx_coarse)/float(Nz_coarse))
                linearTgradient = 1.0/4.0 # delta T / height
                params['main.refine_thresh'] = str(linearTgradient*16.0/float(Nz_coarse))
                params['main.tag_buffer_size'] = str(int(float(Nz_coarse)/8))
                params['main.steady_state'] = '1e-8' #str(1e-4 * pow(32.0/float(Nz_coarse),2))
                params['main.max_grid_size'] = '32'
                params['main.max_step'] = 10000
                params['main.debug'] = 'true' # make sure we output things like Temperature error


            elif physicalProblem == 'convectionDB':

                Nx_coarse = Nz_coarse

                gridFile = mushyLayerBaseDir + '/grids/leftRight/' + str(Nx_coarse) + 'x' + str(Nz_coarse)
                #params = readInputs(mushyLayerBaseDir + '/params/convergenceTest/.parameters')
                if not params_file:
                    params_file = mushyLayerBaseDir + '/params/convergenceTest/convectionDarcyBrinkmanConvTest.parameters'

                params = readInputs(params_file)
                

                params['main.refine_thresh'] = str(3.0/float(Nz_coarse))
                params['main.tag_buffer_size']=str(max(2,int(float(Nz_coarse)/8)))
                

                params['main.steady_state'] = '1e-8'
                params['main.max_time'] = '50000000.0'
                params['main.max_dt'] = '1e-1'

                
                params['main.plot_interval'] = '10000' #str(int(100.0*float(Nz_coarse)/16.0))
                params['main.checkpoint_interval'] = params['main.plot_interval']

                params['main.time_integration_order'] = '2'
                    
     
            elif physicalProblem == 'DBVariablePorosity':
                integrationTime = 2e-4
                
                Nx_coarse = Nz_coarse
                
                if ref_rat == 2:
                    gridFile = mushyLayerBaseDir + '/grids/middleXSmall/' + str(Nx_coarse) + 'x' + str(Nz_coarse)           
            
                else:
                    gridFile = mushyLayerBaseDir + '/grids/middleXSmall/' + str(Nx_coarse*2) + 'x' + str(Nz_coarse*2)   
                            
                    
                #params = readInputs(mushyLayerBaseDir + )
                if not params_file:
                    params_file = mushyLayerBaseDir + '/params/convergenceTest/DBVariablePorosityConvTest.parameters'

                params = readInputs(params_file)
                
                #runTypes = ['uniform', 'amr', 'variable']
                #runTypes = ['uniform', 'amr']
                #runTypes = ['uniform']
                #runTypes = ['amr']
                
                # For a fixed dt:
                dt = (1e-5)*float(16.0/float(Nz_coarse))
                numSteps = float(integrationTime)/dt
                
                params['main.plot_interval'] = str(int(numSteps/5.0))
                
                params['main.max_step'] = '100000'
                
                params['main.min_time'] = integrationTime
                params['main.max_time'] = integrationTime

                params['main.vel_refine_thresh'] = 5.0
                params['main.stdev'] = '0.002'
                
                regrid_int = int(4.0*float(Nx_coarse)/16.0)

                # Make sure we always have one grid
                maxGridSize = max(Nx_coarse*maxRefinement,4)

                #bf = max(maxGridSize/2,4)
                bf = max(Nx_coarse*maxRefinement/8,4)

                gridBuffer = 0 #max(bf/4,1)
                tagBuffer = 0 #gridBuffer

                if max_level == 2:
                    gridBuffer = max(bf/4, 4)

                #if ref_rat == 4:
                #    bf = bf/2
                #   gridBuffer = bf
                #   tagBuffer = gridBuffer


                params['main.block_factor'] = str(bf)
                params['main.grid_buffer_size'] = str(gridBuffer)
                params['main.tag_buffer_size'] = str(tagBuffer)
                params['main.regrid_interval']= str(regrid_int) + ' ' + str(regrid_int) + ' ' + str(regrid_int)
                params['main.max_grid_size'] = str(maxGridSize) # Must be greater than block factor

                # To pretend AMR is VM:
                #params['regrid.tagCenterOnly'] = '0.375'
                #params['main.tag_buffer_size'] = '0'
                #params['main.grid_buffer_size'] = '0'

                chi = '0.05'
                params['bc.porosityHiVal'] = chi + ' ' + chi # 0.4 0.4
                params['bc.porosityLoVal'] = chi + ' ' + chi

                params['main.fixed_dt'] = str(dt)



            elif physicalProblem == 'FixedChill':

                if not params_file:
                    params_file = mushyLayerBaseDir + '/params/convergenceTest/FixedChill.parameters'

                params = readInputs(params_file)

                output_dir = ''

                #Nx_coarse = Nz_coarse*2
                gridFile = mushyLayerBaseDir + '/grids/topMiddle/' + str(Nx_coarse) + 'x' + str(Nz_coarse)

                #pltInt = int(100.0*float(Nz_coarse)/64.0)
                #params['main.plot_interval'] = str(pltInt)
                #params['main.chk_interval'] = params['main.plot_interval']


            elif physicalProblem == 'PorousMushyHole' or physicalProblem == 'MushyConvectionLiquid' or physicalProblem == 'MushyConvectionLiquid2':

                if not params_file:
                    params_file = mushyLayerBaseDir + '/params/convergenceTest/porousMushyHole.parameters'

                params = readInputs(params_file)

                sameDt = 5e-7
                dt64 = sameDt

                # Scale dt with grid spacing
                sameDt=-1

                output_dir = '' # physicalProblem+'-t' + str(maxTime) + '/'
                #output_dir = 'AMRConvergenceTest/'+physicalProblem+'-t' + str(maxTime) + '-1proc/'

                Nx_coarse = Nz_coarse

                Nx_grid = Nx_coarse
                if ref_rat == 4:
                    Nx_grid = 2*Nx_coarse

                gridFile = mushyLayerBaseDir + '/grids/middleXSmall/' + str(Nx_grid) + 'x' + str(Nx_grid)

                dt = dt64*64.0/float(Nx_coarse)
                if sameDt > 0:
                    params['main.fixed_dt'] = sameDt #dt
                else:
                    params['main.fixed_dt'] = dt #dt


                pltInt = int(math.ceil(4*float(Nz_coarse)/32.0))
                pltInt = max(pltInt, 1)
                regrid_int = max(pltInt/4, 1)

                pltInt = 800 # large interval, save disk space
                params['main.plot_interval'] = str(pltInt)
                params['main.checkpoint_interval'] = params['main.plot_interval']


                # AMR options
                if 'main.vel_refine_thresh' in params:
                    params.pop('main.vel_refine_thresh')

                params['main.taggingVar'] = 3 # porosity
                params['main.refine_thresh'] = float(params['main.radius'])*1.0/float(Nz_coarse) #0.05*64.0/float(Nz_coarse)

                bf = max(int(Nx_coarse/4), 4)
                
                maxGridSize = bf*2
                gbs = max(bf/8, 2)

                params['main.block_factor'] = str(bf)
                params['main.grid_buffer_size'] = gbs
                params['main.tag_buffer_size'] = 0
                params['main.regrid_interval']= str(regrid_int) + ' ' + str(regrid_int) + ' ' + str(regrid_int)
                params['main.max_grid_size'] = str(maxGridSize) # Must be greater than block factor

            else:
                print(Fore.RED  + 'Unknown convergence test option' + Style.RESET_ALL)
                sys.exit(2)



            # Default options
            if Nx_coarse == -1:
                print('Trying to convert to an array: ' + str(params['main.num_cells']))
                gridPts = str2arr(params['main.num_cells'])

                aspectRatio = gridPts[0]/gridPts[1]
                Nx_coarse = aspectRatio*Nz_coarse

            if not gridFile:
                gridFile = mushyLayerBaseDir + '/grids/middle/' + str(Nz_coarse) + 'x' + str(Nz_coarse)

            # Always turn off slope limiting for convergence tests
            params['main.use_limiting'] = 'false'
            params['main.debug'] = 'false' # also turn off debug to save disk space

            # Any extra params we may have
            for k,v in extra_params.iteritems():
                params[k] = v

            numCellsAMR = str(Nx_coarse) + ' ' + str(Nz_coarse) + '  8'    
            

            gridFileQuarter = mushyLayerBaseDir + '/grids/rightQuarter/' + str(Nx_coarse) + 'x' + str(Nz_coarse)
            gridFileThreeLevels = mushyLayerBaseDir + '/grids/rightHalfThreeLevel/' + str(Nx_coarse) + 'x' + str(Nz_coarse)
            numCellsUniform = str(Nx_coarse*finestRefinement) + ' ' + str(Nz_coarse*finestRefinement) + '  8'
            numCellsUniform = str(Nx_coarse) + ' ' + str(Nz_coarse) + '  8'

            #params['main.gridfile'] = gridFile
            params['main.num_cells'] = numCellsAMR
            params['main.max_level'] = str(max_level)
            params['main.ref_ratio']= str(ref_rat) + ' ' + str(ref_rat) + ' ' + str(ref_rat) 
            
            num_proc = num_procs[Nz_i]
            params['num_proc'] = num_proc

            if num_proc > 1:
                optimalGrid = float(Nx_coarse)/(float(num_proc)/2.0)

                #print('Initial optimal grid guess: %d' % optimalGrid)
                
                # increase to next power of 2
                while not isPowerOfTwo(optimalGrid):
                    optimalGrid = optimalGrid + 1

                #print('Final optimal grid guess: %d' % optimalGrid)
                
                maxGrid = max(16, optimalGrid)
                # grid size must be greater than the blocking factor 
                maxGrid = max(maxGrid, 2*int(params['main.block_factor'])) 
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
                
                p1['run_name'] = 'AMR-Subcycle-Reflux-Freestream'+str(p1['projection.eta'])+'-MaxLevel' + str(max_level) 
                p1['concise_run_name'] = 'AMR'
                #p13['main.initialize_VD_corr'] = 'true' # true by default
            
                param_sets.append(p1)
            
            
            # Variable mesh
            if 'variable' in runTypes and max_level==1:    
                p2 = dict(params)
                p2['main.reflux_scalar'] = '1'
                p2['main.use_subcycling'] = '1'
                p2['main.refluxType'] = '2'
                p2['projection.eta'] = '0.95'

                p2['main.gridfile'] = gridFile
                p2['run_name'] = 'VariableMesh2SubcycleRefluxFreestream'+str(p2['projection.eta']) 
                p2['concise_run_name'] = 'VM'
            
                param_sets.append(p2)
            
            if 'variable' in runTypes and max_level==2: 
                p3 = dict(params)
                p3['main.reflux_scalar'] = '1'
                p3['main.use_subcycling'] = '1'
                p3['main.refluxType'] = '2'
                p3['projection.eta'] = '0.95'

                p3['main.gridfile'] = gridFileThreeLevels
                p3['run_name'] = 'VM3LevelsSubcycleRefluxFreestream' +str(p3['projection.eta']) 
                p3['concise_run_name'] = 'VM'
                
                param_sets.append(p3)
            
            
            
            # Do this for all the param sets just made:
            for p in param_sets:
                
                
                #if num_proc > 1:
                #    p['run_name'] = 'p' + p['run_name']
                
                run_name = p['run_name']
                run_name = run_name
                if int(p['main.max_level']) > 0:
                     run_name = run_name + '-ref' + str(ref_rat)
                     
                run_name = run_name  + '-' + physicalProblem + '-' + str(Nz_coarse)
                
                p['concise_run_name'] = p['concise_run_name'] + str(Nz_coarse) + physicalProblem
                p['main.plot_prefix'] = run_name + '-'
                p['main.chk_prefix'] =  'chk' + p['main.plot_prefix']
                
                # Now add to the full list of param sets
                all_params.append(p)
            

        # Actually run the convergence test
        full_output_dir = os.path.join(base_dir, output_dir)
        #print('Data dir: ' + base_dir + '\n')
        #print('Output dir: ' + output_dir + '\n')
        #print('Full output dir: ' + full_output_dir + '\n')
        
        
        these_job_ids = AMRConvergenceTest(all_params, full_output_dir, physicalProblem, Nzs, num_procs, numRestarts, analysis_command)

        # Concatenate lists
        job_ids = job_ids + these_job_ids


    # Once all these runs have been submitted, submit the analysis job
    if analysis_command:
        runAnalysisName = 'runAnalysis.sh'

        # Don't redo analysis - we may be waiting on runs to finish
        if os.path.exists(os.path.join(base_dir, runAnalysisName)):
            print(Fore.YELLOW + 'Analysis job already submitted \n' + Fore.RESET)
        else:
            jobName = physicalProblem + '-analysis'

            s = SlurmTask(base_dir, jobName, '', 4)

            s.setDependency(job_ids)
            s.setCustomCommand(analysis_command)

            s.writeSlurmFile(runAnalysisName)
            s.runTask(runAnalysisName)
            print(Fore.GREEN + 'Submitted analysis job \n' + Fore.RESET)


    return job_ids


def main(argv):
    
    # Defaults:
    num_proc = 4
    
    try:
       opts, args = getopt.getopt(argv,"n:f:H:")
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

    base_dir = os.getcwd() # for local runs
    base_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/' # for writing to shared storage

    # Change things below here
    #num_proc = 1
    Nzs = [16, 32,64,128,256,512]
    #Nzs = [128]

    #num_procs = num_proc+0*Nzs
    num_procs = [1,1,2,2,4,4,8,8,8,8,8,8]
    #num_procs = [1,1,4,4,4,4,4,4,4]
    #num_procs = [1,1,1,1,1,1,1,1,1,1,1]

    #num_procs = [16]

    #Nzs = [512,256,128,64,32] # more efficient to do it this way around
    
    # By default, don't try and restart from existing files
    #numRestarts = 0


    #max_level = 0
    #ref_rat = 2

    AMRSetup = [{'max_level': 1, 'ref_rat': 2},
    {'max_level': 1, 'ref_rat': 4},
    {'max_level': 2, 'ref_rat': 2}]

    #AMRSetup = [{'max_level': 1, 'ref_rat': 2}]

    #AMRSetup = [{'max_level': 0, 'ref_rat': 2}]

    #mushQuick, MushyDB, DBVariablePorosity, DirectionalDarcy, convectionDarcyBrinkman, DirectionalDB
    # MushyConvectionAnalyticVel, MushyConvection, MushyConvectionLiquid
    physicalProblem = 'MushyConvectionLiquid2' 

    base_dir = os.path.join(base_dir, 'AMRConvergenceTestMushyConvection')

    runTest(base_dir, physicalProblem, AMRSetup, Nzs, num_procs)

            

if __name__ == "__main__":
    print('Running script')
    main(sys.argv[1:])
                 


