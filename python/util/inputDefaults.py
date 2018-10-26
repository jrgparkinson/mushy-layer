# This is just a list of the default param options
# We will then use this to remove default options from other input files, 
# tidying them up

# Also list 'required' parameters, along with typical values
import re
#parDir = os.path.abspath(os.pardir)
#pythonDir = os.path.join(parDir, 'python')
#sys.path.append(pythonDir)

#from MushyLayerRun import MushyLayerRun
#import math
#from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile

required = {}
default = {}

required['main.cfl'] = 0.8
required['main.use_limiting'] = True
required['main.max_level'] = 0
required['main.regrid_interval'] = [4, 4, 4]
required['main.max_step'] = 0
required['main.tag_buffer_size'] = 4
required['main.refine_thresh'] = 0.3
required['main.block_factor'] = 1
required['main.max_grid_size'] = 32
required['main.fill_ratio'] = 0.75
required['main.checkpoint_interval'] = 0
required['main.plot_interval'] = 0
required['main.max_dt_growth'] = 1.1
required['main.fixed_dt'] = -1
required['main.dt_tolerance_factor'] = 1.1
required['main.grid_buffer_size'] = 1
required['main.verbosity'] = 0
required['main.max_time'] = 0.0
required['main.ref_ratio'] = [2, 2, 2]
required['main.num_cells'] = [64, 64, 64]
required['main.periodic_bc'] = [0, 0, 0]

required['amrmultigrid.num_smooth'] = 4
required['amrmultigrid.num_mg'] = 1
required['amrmultigrid.hang_eps'] = 1e-10
required['amrmultigrid.norm_thresh'] = 1e-10
required['amrmultigrid.tolerance'] = 1e-10
required['amrmultigrid.max_iter'] = 10
required['amrmultigrid.verbosity'] = 0

required['parameters.problem_type'] = 0 # 0 is the mushy layer problem
required['parameters.eutecticTemp'] = -23
required['parameters.eutecticComposition'] = 230
required['parameters.initialComposition'] = 30
required['parameters.liquidusSlope'] = -0.1
required['parameters.waterDistributionCoeff'] = 1e-5

required['bc.scalarLo'] = [0, 0, 0]
required['bc.scalarHi'] = [0, 0, 0]
required['bc.velLo'] = [0, 0, 0]
required['bc.velHi'] = [0, 0, 0]
required['bc.enthalpyLoVal'] = [0.0, 10.0, 0.0]
required['bc.enthalpyHiVal'] = [0.0, -1.0, 0.0]
required['bc.bulkConcentrationLoVal'] = [-1.0, -1.0, -1.0]
required['bc.bulkConcentrationHiVal'] = [-1.0, -1.0, -1.0]


# Now the optional parameters (but which really should be specified)
# along with their default values






default['parameters.fixed_porosity'] = ''
default['parameters.rampBuoyancy'] = 0.0
default['parameters.maxRaT'] = 0.0
default['parameters.maxRaC'] = 0.0
default['parameters.initRaT'] = 0.0
default['parameters.initRaC'] = 0.0
default['parameters.forcing_timescale'] = 0.5
default['parameters.referenceTemperature'] = -23
default['parameters.referenceSalinity'] = 230
default['parameters.prevReferenceTemperature'] = -23
default['parameters.prevReferenceSalinity'] = 230
default['parameters.restart_dtReduction'] = 1.0
default['parameters.'] = ''
default['parameters.'] = ''
default['parameters.'] = ''
default['parameters.'] = ''
default['parameters.'] = ''
default['parameters.'] = ''
default['parameters.'] = ''
default['parameters.prandtl'] = 0.0

default['advSrc.darcy'] = True
default['advSrc.pressure'] = False
default['advSrc.viscous'] = True
default['advSrc.buoyancy'] = True

default['amrmultigrid.MGtype'] = -1
default['amrmultigrid.relaxMode'] = 1 # GSRB
default['init.initVel'] = 0
default['main.initial_cfl'] = 0.1
default['main.analyticSoln'] = -1
default['main.perturbationSin'] = False
default['main.initSummerProfile'] = -1
default['main.summerProfileMushHeight'] = 0.0
default['main.initLinear'] = False
default['main.fixed_porosity'] = -1.0
default['main.porosity_function'] = 0
default['main.maxRestartWavenumber'] = 50
default['main.num_init_passes'] = 2
default['main.meltPondDepth'] = 0
default['main.ignoreBulkConcentrationSteadyState'] = False
default['main.min_time'] = 0.0
default['main.spongeHeight'] = 0.0
default['main.adv_vel_centering'] = 0.5
default['main.usePhiForImplicitAdvectionSolve'] = True
default['main.advPorosityLimit'] = 1e-4
default['main.analyticVel'] = False
default['main.analyticVelProject'] = False 
default['main.freestreamBeforeProjection'] = False
default['main.resetLambda'] = False
default['main.ccvel_porosity_cap'] = 1e-4
default['main.uDeluMethod'] = 0
default['main.uDelu_porosity'] = 1e-15
default['main.maxProjections'] = 1
default['main.maxChi'] = 1.05
default['main.stdev'] = 0.005
default['main.innerRadius'] = 0.2
default['main.vel_src_centering'] = 0.5
default['main.ignoreAdvectionViscousSrc'] = False
default['main.skipHCUpdate'] = False
default['main.analyticSourceTerm'] = False
default['main.diffusiveSrcForAdvection'] = True
default['main.printAccelDt'] = False
default['main.useAccelDt'] = False
default['main.taggingVar'] = 2 #temperature
default['main.refineMethod'] = 0 # undivided gradient
default['main.taggingVectorVar'] = -1 # none
default['main.vel_refine_thresh'] = 0.1
default['main.min_dt'] = 1e-7
default['main.refluxType'] = 2 # nonlinear correction
default['main.domain_length'] = 1.0
default['main.domain_width'] = 1.0
default['main.initAnalyticVel'] = False

default['main.restart_newTime'] = -1.0
default['main.'] = ''
default['main.'] = ''
default['main.'] = ''

default['main.addSubtractGradP'] = True
default['main.advectionInterpOrder'] = 1
default['main.chk_prefix'] = 'chk'
default['main.computeDiagnostics'] = True
default['main.debug'] = False
default['main.doAdvectionSolve'] = True
default['main.doAutomaticRestart'] = False
default['main.doEuler'] = True
default['main.doProjection'] = True
default['main.doScalarAdvectionDiffusion'] = True
default['main.doSyncOperations'] = True
default['main.enforceAnalyticSoln'] = False
default['main.explicitDarcyTerm'] = False
default['main.fixed_dt'] = -1
default['main.gridfile'] = 'grids'
default['main.ignoreSolverFails'] = True
default['main.implicitAdvectionSolve'] = False
default['main.initialize_pressure'] = True
default['main.initialize_VD_corr'] = True
default['main.initialPerturbation'] = 0.0
default['main.iter_plot_interval'] = -1
default['main.m_initial_dt_multiplier'] = 0.1
default['main.nondimensionalisation'] = 0
default['main.output_folder'] = 'output'
default['main.perturbationWavenumber'] = -1
default['main.plot_prefix'] = 'plt'
default['main.porosityInAdvectionOperator'] = True
default['main.project_initial_vel'] = True
default['main.reflux_momentum'] = True
default['main.reflux_normal_momentum'] = True
default['main.reflux_scalar'] = True
default['main.regrid_smoothing_coeff'] = 0.05
default['main.restart_file'] = 'restart.2d.hdf5'
default['main.restartPerturbation'] = 0.0
default['main.scalePwithChi'] = True
default['main.solverFailRestartMethod'] = 0 # restart timestep and halve dt
default['main.splitDarcySolve'] = False
default['main.steady_state'] = 1e-4
default['main.steadyStateNormType'] = 1
default['main.symmetric_domain'] = False
default['main.time_integration_order'] = 1
default['main.use_subcycling'] = False
default['main.viscous_num_smooth_down'] = 8
default['main.viscous_num_smooth_up'] = 8
default['main.viscous_solver_tol'] = 1e-10
default['main.viscousBCs'] = float(default['parameters.prandtl']) > 0
default['projection.usePiAdvectionBCs'] = True
default['projection.verbosity'] = 0
default['projection.scaleMACBCWithChi'] = False
default['projection.MACbcScale'] = 1.0

default['regrid.tagMushLiquidBoundary'] = False
default['regrid.tagDomainBoundary'] = False
default['regrid.marginalPorosityLimit'] = 0.9
default['regrid.tagCenterOnly'] = 0
default['regrid.testRegridCoarsening'] = False

# Test we haven't missed anything
files = ['/home/parkinsonjl/convection-in-sea-ice/MushyLayer/srcSubcycle/MushyLayerSubcycleUtils.cpp',
'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/srcSubcycle/AMRLevelMushyLayerSync.cpp',
'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/srcSubcycle/AMRLevelMushyLayerRegrid.cpp',
'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/srcSubcycle/AMRLevelMushyLayerIO.cpp',
'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/srcSubcycle/AMRLevelMushyLayerInit.cpp',
'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/srcSubcycle/AMRLevelMushyLayerFactory.H',
'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/srcSubcycle/AMRLevelMushyLayerFactory.cpp',
'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/srcSubcycle/AMRLevelMushyLayer.H',
'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/srcSubcycle/AMRLevelMushyLayer.cpp']

for file in files:

	print('Comparing with ' + file)
	#regexStr = '(pp[a-zA-Z]?+)\.[a-zA-Z]+\(\"([a-zA-Z\-\_]+)\".*'
	regexStr = '(p[^\.\n]+)\.[a-zA-Z]+\(\"([a-zA-Z\-\_]+)\".*'
	f = open(file, 'r')
	lines = f.readlines()
	for line in lines:
		m = re.findall(regexStr, line)

		for match in m:
			#print(match)

			#if match[0] == 'pp':
			#	prefix = 'main'
			#else:
			#	prefix = match[0][2:]

			#fullKey = prefix # + '.' + match[1]
			#print(fullKey)

			existsInKeys = False
			for k in default.keys():
				if match[1] in k:
					existsInKeys = True

			for k in required.keys():
				if match[1] in k:
					existsInKeys = True

					#print(k)

			if existsInKeys == False:
				print(match)
				print("Not found!")

# Test by making inputs from just required options
loc = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/inputs'
writeParamsFile(loc, required)

