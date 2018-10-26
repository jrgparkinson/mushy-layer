# Script to run mushy layer simulations for a variety of 
# grid setups to test how reflux performs/helps

import os
import sys
parDir = os.path.abspath(os.pardir)
pythonDir = os.path.join(parDir, 'python')
sys.path.append(pythonDir)

from MushyLayerRun import MushyLayerRun
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile

#################################
# These shouldn't change
#################################
#cwd = '/home/parkinsonj/convection-in-sea-ice/MushyLayer'
cwd = os.getcwd()
mushyLayerBaseDir = cwd.replace('/run', '')

os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
subcycled = True # subcycled code

  #output_dir = 'mushyLayerLowC-lowerBranch'
    
#output_dir = 'refluxTest'
#output_dir = 'refluxTestMushyLayer'
data_dir = os.getcwd() # for local runs
data_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/' # for writing to shared storage
output_dir = 'AMRConvergenceTestMushyLayer'
output_dir = 'AMRConvergenceTestHRL'
output_dir = 'AMRConvergenceTestNoFlow'
    
# Parallel stuff
num_proc = 1

Nzs = [8, 16] # for HRL etc
Nzs = [8, 16, 32, 64, 128, 256] # for mushy layer

# Grid stuff
for Nz_coarse in Nzs:
    #Nz_coarse = 32
    
    max_level = 2
    ref_rat = 2
    finestRefinement = pow(ref_rat, max_level)

    # Want to do several runs with slightly different parameters

    # Common parameters:
    
#     physicalProblem = 'MushyLayerDarcy'
#     params = readInputs(mushyLayerBaseDir + '/params/mushyLayerDarcy.parameters')
#     params['main.multiCompSolve'] = 'true'
#     params['picardSolver.max_iter'] = '1'
#     params['main.periodic_bc'] = '0 0 1'
#     params['main.steady_state']='1e-6'
#     params['main.checkpoint_interval'] = '2000'
#     params['main.plot_interval'] = '500'
#     params['main.max_grid_size'] = '256'
#     params['main.domain_width'] = '0.1'
#     params['main.enforceAnalyticSoln'] = 'true'
#     params['main.ignoreSolverFails'] = 'true'
#     params['amrmultigrid.verbosity'] = '0' # 2 for some logging
#     params['amrmultigrid.tolerance']   = '1.0e-10'
#     params['amrmultigrid.norm_thresh'] = '1.0e-10'
#     params['amrmultigrid.hang_eps']    = '1.0e-10'  
#     params['main.initialPerturbation'] = '0.0' 
#     params['main.fixed_dt'] = '0.00001' 
#     params['main.max_time'] = '0.01'
#     params['projector.verbosity'] = '0'
#       
#     params['parameters.rayleighComp'] = '50.0'
#     params['parameters.referenceTemperature'] = '-23'
#     params['parameters.referenceSalinity'] = '230'
#     params['parameters.topEnthalpy'] = '0.0'
#     params['parameters.bottomEnthalpy'] = '6.1'           
#     params['parameters.topBulkConc'] = '-1.0'
#     params['parameters.bottomBulkConc'] = '-1.0'    
#     params['parameters.stefan'] = '5'
#     params['parameters.compositionRatio'] = '5.0'
#     params['parameters.K'] = '1'
#     params['parameters.specificHeatRatio'] = '1'
#     params['parameters.lewis'] = '100'
#     params['parameters.heatConductivityRatio'] = '1'
#     Nx_coarse = Nz_coarse;
    
    physicalProblem = 'diffusiveSolidification'
    params = readInputs(mushyLayerBaseDir + '/params/solidificationNoFlow.parameters')
    params['main.multiCompSolve'] = 'true'
    params['picardSolver.max_iter'] = '1'
    params['main.periodic_bc'] = '1 0 1'
    params['main.steady_state']='1e-5'
    params['main.checkpoint_interval'] = '2000'
    params['main.plot_interval'] = '100'
    params['main.max_grid_size'] = '256'
    params['main.domain_width'] = '0.125'
    params['main.enforceAnalyticSoln'] = 'true'
    params['main.cfl'] = '0.05'
    params['main.ignoreSolverFails'] = 'true'
    params['amrmultigrid.verbosity'] = '0' # 2 for some logging
    params['amrmultigrid.tolerance']   = '1.0e-10'
    params['amrmultigrid.norm_thresh'] = '1.0e-10'
    params['amrmultigrid.hang_eps']    = '1.0e-10'  
    params['main.initialPerturbation'] = '0.0' #str(0.001*float(32)/float(Nx_coarse)) #'0.001'
    params['main.refine_thresh'] = str(0.3*16.0/float(Nz_coarse))
    Nx_coarse = int(Nz_coarse/4);

    
#     physicalProblem = 'HRL'
#     params = readInputs(mushyLayerBaseDir + '/params/hrl.parameters')
#     params['main.multiCompSolve'] = 'true'
#     params['picardSolver.max_iter'] = '1'
#     params['main.periodic_bc'] = '0 0 1'
#     params['main.steady_state']='1e-6'
#     params['main.checkpoint_interval'] = '2000'
#     params['main.plot_interval'] = '50'
#     params['main.max_grid_size'] = '256'
#     params['main.domain_width'] = '1.0'
#     params['main.enforceAnalyticSoln'] = 'true'
#     params['main.ignoreSolverFails'] = 'true'
#     params['amrmultigrid.verbosity'] = '0' # 2 for some logging
#     params['amrmultigrid.tolerance']   = '1.0e-15'
#     params['amrmultigrid.norm_thresh'] = '1.0e-15'
#     params['amrmultigrid.hang_eps']    = '1.0e-15'  
#     params['main.initialPerturbation'] = '0.0' 
#     params['main.fixed_dt'] = '0.00001' 
#     #params['main.fixed_dt'] = '0.0001' 
#     #params['main.max_time'] = '0.0005' #short
#     params['main.max_time'] = '0.1'#long
#     params['parameters.rayleighTemp'] = '500'
#     params['projector.verbosity'] = '0'
#     Nx_coarse = Nz_coarse;
#     
    
    numCellsAMR = str(Nx_coarse) + ' ' + str(Nz_coarse) + '  8'    
    gridFile = mushyLayerBaseDir + '/grids/topQuarter/' + str(Nx_coarse) + 'x' + str(Nz_coarse)
    numCellsUniform = str(Nx_coarse*finestRefinement) + ' ' + str(Nz_coarse*finestRefinement) + '  8'
    
    gridFile = ''
    #params['main.gridfile'] = gridFile
    params['main.num_cells'] = numCellsAMR
    params['main.max_level'] = str(max_level)
    

    # Now construct vector different param sets
    param_sets = []

    
    # Run 1: uniform mesh
    p0 = dict(params)
    p0['main.max_level'] = '0'
    p0['main.num_cells'] = numCellsUniform
    p0['run_name'] = 'Uniform'
    param_sets.append(p0)

    # Run 2: AMR/no subcycling/no reflux
    p1 = dict(params)
    p1['main.reflux_scalar'] = '0'
    p1['main.use_subcycling'] = '0'
    p1['projection.eta'] = '0.0'
    p1['run_name'] = 'AMRNoSubcycleNoReflux'
    #param_sets.append(p1)

    # Run 3: AMR/no subcycling/reflux
    p2 = dict(params)
    p2['main.reflux_scalar'] = '1'
    p2['main.use_subcycling'] = '0'
    p2['main.refluxType'] = '2'
    p2['projection.eta'] = '0.0'
    p2['run_name'] = 'AMRNoSubcycleReflux'
    #param_sets.append(p2)

    # Run 4: AMR/subcycling/no reflux
    p3 = dict(params)
    p3['main.reflux_scalar'] = '0'
    p3['main.use_subcycling'] = '1'
    p3['projection.eta'] = '0.0'
    p3['run_name'] = 'AMRSubcycleNoReflux'
    #param_sets.append(p3)

    # Run 5: AMR/subcycling/reflux
    p4 = dict(params)
    p4['main.reflux_scalar'] = '1'
    p4['main.use_subcycling'] = '1'
    p4['main.refluxType'] = '2'
    p4['projection.eta'] = '0.0'
    p4['run_name'] = 'AMRSubcycleReflux'
    param_sets.append(p4)
    
        # Run 3: AMR/no subcycling/reflux/freestream
    p10 = dict(params)
    p10['main.reflux_scalar'] = '1'
    p10['main.use_subcycling'] = '0'
    p10['main.refluxType'] = '2'
    p10['run_name'] = 'AMRNoSubcycleRefluxFreestream0.5'
    p10['projection.eta'] = '0.5'
    #param_sets.append(p10)
    
#    p11 = dict(params)
#    p11['main.max_level'] = str(max_level)
#    p11['main.gridfile'] = mushyLayerBaseDir + '/grids/topHalf/' + str(Nx_coarse) + 'x' + str(Nz_coarse)
#    p11['main.num_cells'] = numCellsAMR
#    p11['main.reflux_scalar'] = '1'
#    p11['main.use_subcycling'] = '1'
#    p11['main.refluxType'] = '2'
#    p11['run_name'] = 'AMRSubcycleRefluxFreestream'
#    p11['projection.eta'] = '0.55'
#    param_sets.append(p11)
    
    
    p13 = dict(p4)
    p13['run_name'] = 'AMRSubcycleRefluxFreestreamNoInit0.5'
    p13['projection.eta'] = '0.5'
    p13['main.initialize_VD_corr'] = 'false' # true by default
    #param_sets.append(p13)
    
    p13 = dict(p4)
    p13['run_name'] = 'AMRSubcycleRefluxFreestream0.5'
    p13['projection.eta'] = '0.5'
    #param_sets.append(p13)
    
    p14 = dict(p4)
    p14['run_name'] = 'AMRSubcycleRefluxFreestream0.49'
    p14['projection.eta'] = '0.49'
    #param_sets.append(p14)
    
    p15 = dict(p4)
    p15['run_name'] = 'AMRSubcycleRefluxFreestream0.51'
    p15['projection.eta'] = '0.51'
    #param_sets.append(p15)
    
    #p14 = dict(p10)
    #p14['run_name'] = 'AMRNoSubcycleRefluxFreestream0.52'
    #p14['projection.eta'] = '0.52'
    #param_sets.append(p14)
    
    #p14 = dict(p11)
    #p14['run_name'] = 'AMRSubcycleRefluxFreestream0.6'
    #p14['projection.eta'] = '0.6'
    #param_sets.append(p14)
    
    # Run 6: AMR/no subcycling/linear reflux
    p5 = dict(params)
    p5['main.reflux_scalar'] = '1'
    p5['main.use_subcycling'] = '0'
    p5['main.refluxType'] = '1'
    p5['run_name'] = 'AMRNoSubcycleLinearReflux'
    #param_sets.append(p5)
    
    # Run 7: AMR/subcycling/linear reflux
    p6 = dict(params)
    p6['main.reflux_scalar'] = '1'
    p6['main.use_subcycling'] = '1'
    p6['main.refluxType'] = '1'
    p6['run_name'] = 'AMRSubcycleLinearReflux'
    #param_sets.append(p6)
    

    
  
    

    for i in range(0, len(param_sets)):
        p = param_sets[i]
        run_name = p['run_name']
        run_name = run_name + '-' + physicalProblem + '-' + str(Nz_coarse)
        p['main.plot_prefix'] = output_dir + '-' + run_name + '-'
        p['main.chk_prefix'] =  'chk' + p['main.plot_prefix']
    
        # Don't repeat runs
        
        #cwd = os.getcwd()
        test_name = run_name + '-0'
        full_path = os.path.join(data_dir, output_dir, test_name)
        print(full_path + '\n')
        if os.path.exists(full_path):
            print('Run already done \n')
            continue
    
        # Do run
        full_output_dir =  os.path.join(data_dir, output_dir)

        ml_run = MushyLayerRun(full_output_dir, num_proc, p, subcycled)
        ml_run.single_run(run_name)



    
        

