
#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif



#include "MushyLayerSubcycleUtils.H"
#include "ParmParse.H"
#include "CoarseAverage.H"
#include "computeNorm.H"
#include "mushyLayerOpt.h"

#include "NamespaceHeader.H"

void
getAMRFactory(RefCountedPtr<AMRLevelMushyLayerFactory>&  a_fact)
{
  CH_TIME("getAMRFactory");

  // First, check for deprecated options
  // These are no longer read in, but used to be, so the user might be expecting some behaviour
  // which they won't get - we better warn them about this.

  ParmParse pp;
  Vector<string> old_options;
  old_options.push_back(string("main.explicitDarcyTerm"));
  old_options.push_back(string("main.lapVelNumSmooth"));
  old_options.push_back(string("main.lapVelSmoothScale"));
  old_options.push_back(string("advSrc.allow_lagged_lap_vel"));
  old_options.push_back(string("ccSrc.darcy"));

  old_options.push_back(string("main.postTraceSmoothing"));
  old_options.push_back(string("main.skipNewLevelScalars"));

  old_options.push_back(string("main.variable_eta_factor"));
  old_options.push_back(string("main.single_level_lambda_corr"));
  old_options.push_back(string("amrmultigrid.relaxMode"));
  old_options.push_back(string("main.implicitAdvectionSolve"));

  old_options.push_back(string("main.usePhiForImplicitAdvectionSolve"));
  old_options.push_back(string("main.noMultigrid"));
  old_options.push_back(string("main.noMultigridIter"));
  old_options.push_back(string("main.init_increase_dt"));
  old_options.push_back(string("main.analyticSourceTerm"));
  old_options.push_back(string("parameters.forcing_timescale"));

  old_options.push_back(string("projection.pre_smoothing"));

  old_options.push_back(string("main.skipTrickySourceTime"));
  old_options.push_back(string("main.iter_plot_interval"));


  for (int i=0; i< old_options.size(); i++)
  {
    string option = old_options[i];
    if (pp.contains(option))
    {
      string warning_msg = "Warning, using deprecated option: " + option;
      MayDay::Warning(warning_msg.c_str());
    }
  }

  ParmParse ppMain("main");

  ParmParse ppParams("parameters");
  ParmParse ppAMRMultigrid("amrmultigrid");
  ParmParse ppProjection("projection");
  ParmParse ppInit("init");
  ParmParse ppAdvsrc("advSrc");
  ParmParse ppCCSrc("ccSrc");
  ParmParse ppRegrid("regrid");
  ParmParse ppPatchGodunov("patchGodunov");
  ParmParse ppVelMultigrid("VelocityMultigrid");
  ParmParse HCMultigrid("HCMultigrid");

  MushyLayerOptions opt;

  MushyLayerParams params;
  params.getParameters();

  opt.cfl = 0.5;
  ppMain.query("cfl",opt.cfl);

  opt.min_time = 0;
  ppMain.query("min_time", opt.min_time);

  opt.output_dir = "";
  ppMain.query("output_folder", opt.output_dir);

  opt.plotfile_prefix = "plt";
  ppMain.query("plot_prefix", opt.plotfile_prefix);

  opt.minimalOutput = false;
  opt.debug = false;
  ppMain.query("debug", opt.debug);
  if (!opt.debug)
  {
    ppMain.query("minimalOutput", opt.minimalOutput);
  }

  opt.steadyStateCondition = 1e-3;
  ppMain.query("steady_state", opt.steadyStateCondition);

  opt.ignoreVelocitySteadyState = !params.isDarcyBrinkman();
  ppMain.query("ignoreVelocitySteadyState", opt.ignoreVelocitySteadyState);

  opt.ignoreBulkConcSteadyState=false;
  ppMain.query("ignoreBulkConcentrationSteadyState", opt.ignoreBulkConcSteadyState);

  // This is really the domain width, not length,
  // but changing it in the inputs files would be a right pain
  // at this point
  opt.domainWidth = -1;
  ppMain.query("domain_length", opt.domainWidth); // retained for backward compatability
  ppMain.query("domain_width", opt.domainWidth);

  std::vector<int> num_cells; // (num_read_levels,1);
  ppMain.getarr("num_cells",num_cells,0,SpaceDim);

  if (ppMain.contains("domain_height"))
  {
    ppMain.query("domain_height", opt.domainHeight);
    opt.domainWidth = opt.domainHeight*num_cells[0]/num_cells[SpaceDim-1];

  }
  else
  {
    // compute domainHeight from domainWidth
    opt.domainHeight = opt.domainWidth*num_cells[SpaceDim-1]/num_cells[0];
  }

  if (opt.domainWidth <= 0)
  {
    MayDay::Error("No domain width specified, or domain width is invalid");
  }


  opt.verbosity = 1;
  ppMain.query("verbosity", opt.verbosity);

  // 1 for volume averaged, 0 for max
  opt.steadyStateNormType = 1;
  ppMain.query("steadyStateNormType", opt.steadyStateNormType);

  // For use when restarting.
  opt.load_advVel = false;
  ppMain.query("load_advVel", opt.load_advVel);

  /**
   * Timestepping
   */
  opt.fixedDt = -1;
  ppMain.query("fixed_dt", opt.fixedDt);

  opt.max_dt_growth = 1.1;
  ppMain.query("max_dt_growth", opt.max_dt_growth);

  opt.init_dt_scale = 0.1;
  ppMain.query("init_dt_scale", opt.init_dt_scale);

  opt.ignoreSolveFails = true;
  ppMain.query("ignoreSolverFails", opt.ignoreSolveFails);

  opt.solverFailRestartMethod = 0;
  ppMain.query("solverFailRestartMethod", opt.solverFailRestartMethod);

  opt.adv_vel_centering_growth = 1.01;
  ppMain.query("adv_vel_centering_growth", opt.adv_vel_centering_growth);

  opt.initAdvVelCentering = 0.1;
  ppMain.query("init_advVel_centering", opt.initAdvVelCentering);

  opt.initial_dt_multiplier = 0.1;
  ppMain.query("initial_cfl", opt.initial_dt_multiplier);

  opt.computeDiagnostics = true;
  ppMain.query("computeDiagnostics", opt.computeDiagnostics);

  // Only relevant for darcy brinkman
  opt.doEulerPart = true;
  ppMain.query("doEuler", opt.doEulerPart);

  opt.doScalarAdvectionDiffusion = true;
  ppMain.query("doScalarAdvectionDiffusion", opt.doScalarAdvectionDiffusion);

  int advectionMethod =  velocityAdvectionTypes::m_porosityOutsideAdvection;
  ppMain.query("advectionMethod", advectionMethod);
  opt.advectionMethod = velocityAdvectionTypes(advectionMethod);

//  opt.skipTrickySourceTerm = -1;
//  ppMain.query("skipTrickySourceTime", opt.skipTrickySourceTerm);

  opt.allowMulticompAdvection = true;
  ppMain.query("allowMulticompAdvection", opt.allowMulticompAdvection);


  opt.scaleP_MAC = true;
  ppMain.query("scalePwithChi", opt.scaleP_MAC);
  // Makes more sense to store this parameter under projection
  ppProjection.query("scalePressureWithPorosity", opt.scaleP_MAC);

  opt.scaleP_CC = opt.scaleP_MAC;
  ppProjection.query("scaleCCPressure", opt.scaleP_CC);

  opt.projection_verbosity = 0;
  ppProjection.query("verbosity", opt.projection_verbosity);

  opt.usePiAdvectionBCs = true;
  ppProjection.query("usePiAdvectionBCs", opt.usePiAdvectionBCs);

//  opt.explicitDarcyTerm = false;
//  ppMain.query("explicitDarcyTerm", opt.explicitDarcyTerm);

  opt.solidPorosity = 0.05;
  ppMain.query("solidPorosity", opt.solidPorosity);

  opt.advPorosityLimit = opt.solidPorosity;
  ppMain.query("advPorosityLimit", opt.advPorosityLimit);

  opt.chiLimit = 0.0;
  ppMain.query("porousAdvVelLimit", opt.chiLimit);

  opt.ccVelPorosityLimit = opt.solidPorosity;
  ppMain.query("ccvel_porosity_cap", opt.ccVelPorosityLimit);

  opt.advVelsrcChiLimit = opt.solidPorosity; //1e-10
  ppMain.query("advVelSrcChiLimit", opt.advVelsrcChiLimit);

  opt.advVelChiLimit = min(pow(10,5)*opt.lowerPorosityLimit, pow(10,-10)) ; //was 1e-10
  ppMain.query("advPorosityLimit", opt.advVelChiLimit);



  // by default make this tiny (so essentially turned off)
  opt.uDelU_porosityLimit = 10*opt.lowerPorosityLimit; //1e-15;
  ppMain.query("uDelu_porosity", opt.uDelU_porosityLimit);

  opt.lapVelBCOrder = 0;
  ppMain.query("lap_vel_bc_order", opt.lapVelBCOrder);

  // Only do multiple projections on the base level
  opt.maxProjBaseLevel = 1;
  ppMain.query("maxProjections", opt.maxProjBaseLevel);

  opt.maxNumMACProj = 1;
  ppMain.query("max_num_MAC_proj", opt.maxNumMACProj);

  opt.enforceAnalyticVel = false;
  ppMain.query("analyticVel", opt.enforceAnalyticVel);

  opt.projectAnalyticVel = false;
  ppMain.query("correctAnalyticVel", opt.projectAnalyticVel);

  opt.analyticVelType = params.physicalProblem;
  ppMain.query ("analyticVelType", opt.analyticVelType);

  opt.initAnalyticVel=false;
  ppMain.query("initAnalyticVel", opt.initAnalyticVel);




  opt.useOldAdvVel = false;
  ppMain.query("useOldAdvVelForTracing", opt.useOldAdvVel);

//  opt.lapVelNumSmooth = 0;
//  ppMain.query("lapVelNumSmooth", opt.lapVelNumSmooth);
//
//  opt.lapVelSmoothScale = 0.0;
//  ppMain.query("lapVelSmoothScale", opt.lapVelSmoothScale);

  opt.CCVelSrcTermCentering = 0.5;
  ppMain.query("vel_src_centering", opt.CCVelSrcTermCentering);

//  opt.advSrcAllowLaggedLapVel = false;
//  ppAdvsrc.query("allow_lagged_lap_vel", opt.advSrcAllowLaggedLapVel);


  opt.advVelPressureSrc = false;
  opt.advVelDarcySrc = true;
  opt.advVelViscousSrc = true;
  opt.advVelBuoyancySrc = true;

  ppAdvsrc.query("pressure", opt.advVelPressureSrc);
  ppAdvsrc.query("darcy", opt.advVelDarcySrc);
  ppAdvsrc.query("viscous", opt.advVelViscousSrc);
  ppAdvsrc.query("buoyancy", opt.advVelBuoyancySrc);


  opt.uDeluMethod = 0;
  ppMain.query("uDeluMethod", opt.uDeluMethod);

  opt.uDelU_grow = 1;
  ppMain.query("uDelU_grow", opt.uDelU_grow);

  opt.uDelUConservativeForm = false;
  ppMain.query("uDelU_conservativeForm", opt.uDelUConservativeForm);

  // Cell-centred velocity source term

  opt.CCAdvSrc = true;
  ppCCSrc.query("advection", opt.CCAdvSrc);

  // This line must come after quering explicit Darcy Term
//  opt.CCDarcySrc = opt.explicitDarcyTerm;
//  ppCCSrc.query("darcy", opt.CCDarcySrc);

  opt.CCPressureSrcOverride = false;
  opt.CCPressureSrc = true;
  if (ppCCSrc.contains("pressure"))
  {
    opt.CCPressureSrcOverride = true;
    ppCCSrc.get("pressure", opt.CCPressureSrc);
  }


  opt.CCBuoyancySrc = true;
  ppCCSrc.query("buoyancy", opt.CCBuoyancySrc);

  opt.spongeHeight = 0.0;
  ppMain.query("spongeHeight", opt.spongeHeight);

//  opt.postTraceSmoothing = 0.0; // default is no smoothing
//  ppMain.query("postTraceSmoothing", opt.postTraceSmoothing);

  opt.rampBuoyancy = 0.0;
  ppParams.query("rampBuoyancy", opt.rampBuoyancy);

  opt.maxRaC = 0.0;
  opt.maxRaT = 0.0;

   ppParams.query("maxRaT", opt.maxRaT);
   ppParams.query("maxRaC", opt.maxRaC);

   opt.initRaC = 0;
   opt.initRaT = 0;
   ppParams.query("initRaT", opt.initRaT);
   ppParams.query("initRaC", opt.initRaC);

//   opt.skipNewLevelScalars = false;
//   ppMain.query("skipNewLevelScalars", opt.skipNewLevelScalars);

   opt.skipSaltUpdate = false;
   ppMain.query("skipSaltUpdate", opt.skipSaltUpdate);

   // Option here to not update enthalpy and salinity
   // (useful for debugging)
   opt.skipHCUpdate = false;
   ppMain.query("skipHCUpdate", opt.skipHCUpdate);

   opt.doDiffusionSrc = true;
   ppMain.query("diffusiveSrcForAdvection", opt.doDiffusionSrc);

   ppMain.query("consider_u_chi_dt", opt.forceUseUChiForCFL);

   opt.skipUnsafeCFL=false;
   ppMain.query("skip_unsafe_cfl", opt.skipUnsafeCFL);

   /**
    * Physics related options
   */

//  opt.iter_plot_interval = -1;
//  ppMain.query("iter_plot_interval", opt.iter_plot_interval);

  opt.lowerPorosityLimit = 1e-15;
  ppMain.query("lowPorosityLimit", opt.lowerPorosityLimit);

  /**
   * AMR related options
   */
  opt.useSubcycling = true;
  ppMain.query("use_subcycling", opt.useSubcycling);

  // 1st/2nd order interpolation for advection
  opt.CFinterpOrder_advection = 1;
  ppMain.query("advectionInterpOrder", opt.CFinterpOrder_advection);


  opt.refineThresh = 0.3;
  ppMain.query("refine_thresh", opt.refineThresh);

  opt.tagBufferSize = 4;
  ppMain.query("tag_buffer_size", opt.tagBufferSize);

  opt.doProjection=true;
  ppMain.query("doProjection", opt.doProjection);

  opt.useIncrementalPressure = false;
  ppProjection.query("useIncrementalPressure", opt.useIncrementalPressure);

//  opt.phiScale = 1;
//  ppProjection.query("phiScale", opt.phiScale);

//  opt.scaleMACBCWithChi = false;
//  opt.MACBCscale = 1.0;
//  ppProjection.query("scaleMACBCWithChi", opt.scaleMACBCWithChi);
//  ppProjection.query("MACbcScale", opt.MACBCscale);

  opt.doSyncOperations = true;
  ppMain.query("doSyncOperations", opt.doSyncOperations);

  opt.do_postRegrid_smoothing = false;
  ppMain.query("do_postRegrid_smoothing", opt.do_postRegrid_smoothing);

  opt.reflux_momentum = true;
  ppMain.query("reflux_momentum", opt.reflux_momentum);

  opt.reflux_normal_momentum = true;
  ppMain.query("reflux_normal_momentum", opt.reflux_normal_momentum);

  opt.reflux_enthalpy = true;
  opt.reflux_concentration = true;
  opt.reflux_lambda = true;
  ppMain.query("reflux_enthalpy", opt.reflux_enthalpy);
  ppMain.query("reflux_concentration", opt.reflux_concentration);
  ppMain.query("reflux_lambda", opt.reflux_lambda);

  opt.refluxAverageDown = true;
  ppMain.query("reflux_average_down", opt.refluxAverageDown);

  int refluxType = 2;
  ppMain.query("refluxType", refluxType);
  opt.refluxMethod = RefluxMethod(refluxType);

  opt.refluxBetaSign = -1;
  opt.refluxCorrSign = 1;
  ppMain.query("refluxBetaSign", opt.refluxBetaSign);
  ppMain.query("refluxCorrSign", opt.refluxCorrSign);

//  opt.variable_eta_factor = 1.0;
//  ppMain.query("variable_eta_factor", opt.variable_eta_factor);
//  CH_assert(opt.variable_eta_factor >= 1); // This must be >= 1, else eta will increase when it should be decreasing (and vice-versa)

//  opt.minEta = 0.99;
//  ppProjection.query("eta", opt.minEta);

//  opt.computeFreestreamCorrectionSingleLevel = false;
//  ppMain.query("single_level_lambda_corr", opt.computeFreestreamCorrectionSingleLevel);

  opt.compute_initial_VD_corr = true;
  ppMain.query("initialize_VD_corr", opt.compute_initial_VD_corr);

  opt.taggingVar = -1; // m_temperature
  opt.taggingVectorVar = -1;
  int tagMethod = 0;

  ppMain.query("taggingVar", opt.taggingVar);
  ppMain.query("refineMethod", tagMethod);
  ppMain.query("taggingVectorVar", opt.taggingVectorVar);

  opt.taggingMethod = TaggingMethod(tagMethod);

  opt.tagMLboundary = false;
  ppRegrid.query("tagMushLiquidBoundary", opt.tagMLboundary);

  opt.tagDomainBoundary = false;
  ppRegrid.query("tagDomainBoundary", opt.tagDomainBoundary);

  opt.fixed_grid_time = -1.0;
  ppRegrid.query("fixed_grid_time", opt.fixed_grid_time);

  opt.min_regrid_time = -1;
  ppRegrid.query("min_regrid_time", opt.min_regrid_time);

  opt.tag_velocity=false;
  if (ppMain.contains("vel_refine_thresh"))
  {
    opt.tag_velocity=true;
    ppMain.get("vel_refine_thresh", opt.vel_thresh);
  }

  opt.tag_plume_mush = false;
  if ((ppRegrid.contains("plume_vel") || ppRegrid.contains("plume_salinity")))
  {
    opt.tag_plume_mush = true;
    opt.plumeSalinityThreshold = -1.0 + log10(params.compositionRatio); // rough guess of salinity in channels
    opt.plumeVelThreshold = params.m_buoyancySCoeff/params.m_darcyCoeff; // rough guess of the velocity in the channels
    ppRegrid.query("plume_vel", opt.plumeVelThreshold); // think  this scales like da^3*ra
    ppRegrid.query("plume_salinity", opt.plumeSalinityThreshold);

  }
  
  opt.compositeChannelTagging=false;
  ppRegrid.query("compositeChannelTagging", opt.compositeChannelTagging);

  opt.taggingMarginalPorosityLimit = 1.0;
  ppRegrid.query("marginalPorosityLimit", opt.taggingMarginalPorosityLimit);

  opt.tagCenterBoxRegridTime = 0.0;
  ppRegrid.query("initTime", opt.tagCenterBoxRegridTime);

  opt.tagCenterBoxSize = 0;
  ppRegrid.query("tagCenterOnly", opt.tagCenterBoxSize);

  opt.testRegridCoarsening = false;
  ppRegrid.query("testRegridCoarsening", opt.testRegridCoarsening);

  opt.vectorHOinterp = false;
  opt.scalarHOinterp = true;
  ppMain.query("vectorHOinterp", opt.vectorHOinterp);
  ppMain.query("scalarHOinterp", opt.scalarHOinterp);

  opt.makeRegridPlots = false;
  ppMain.query("regrid_plots", opt.makeRegridPlots);

  opt.regrid_dt_scale = 0.1;
  ppMain.query("regrid_dt_scale", opt.regrid_dt_scale);

  opt.initLambdaPostRegrid = true;
  ppMain.query("initLambda", opt.initLambdaPostRegrid);

  opt.regrid_freestream_subcycle = true;
  ppMain.query("regrid_freestream_subcycle", opt.regrid_freestream_subcycle);

  opt.regrid_advect_before_freestream=false;
  ppMain.query("regrid_advect_before_freestream", opt.regrid_advect_before_freestream);

  opt.regrid_eta_scale=1.0;
  ppMain.query("regrid_eta_scale", opt.regrid_eta_scale);

  opt.regrid_smoothing_coeff = 0.05;
  ppMain.query("regrid_smoothing_coeff", opt.regrid_smoothing_coeff);

  //Initialization
  if (params.isDarcyBrinkman())
  {
    opt.initialize_pressures = true;
    opt.project_initial_vel = true;
  }
  else
  {
    opt.initialize_pressures = true;
    opt.project_initial_vel = false;
  }

  ppMain.query("initialize_pressure", opt.initialize_pressures);
  ppMain.query("project_initial_vel", opt.project_initial_vel);

  opt.usePrevPressureForUStar = true;
  ppMain.query("addSubtractGradP", opt.usePrevPressureForUStar);

  opt.customInitData = -1;
  ppMain.query("initData", opt.customInitData);

  opt.writePressureInitFields = false;
  ppMain.query("writePressureInitFields", opt.writePressureInitFields);

  opt.initResetStates = true;
  ppMain.query("init_resetStates", opt.initResetStates);

  /**
   * Solver options
   */
  opt.useLimiting = true;
  ppMain.query("use_limiting", opt.useLimiting);

  opt.viscous_solver_tol = 1e-10;
  ppMain.query("viscous_solver_tol", opt.viscous_solver_tol);

  opt.viscous_num_smooth_down = 8;
  ppMain.query("viscous_num_smooth_down", opt.viscous_num_smooth_down);

  opt.viscous_num_smooth_up = 8;
  ppMain.query("viscous_num_smooth_up", opt.viscous_num_smooth_up);

  opt.timeIntegrationOrder = 1;
  ppMain.query("time_integration_order", opt.timeIntegrationOrder);

  int mgtype = MGmethod::MGTypeFAS;

//  opt.AMRMultigrid_verbosity = 0;
//  opt.AMRMultigridRelaxMode = 1; // 1=GSRB, 4=jacobi
  opt.AMRMultigridVerb=0;
  opt.AMRMultigridTolerance=1e-10;
  opt.AMRMultigridHang=1e-10;
  opt.AMRMultigridNormThresh=1e-10;
//  ppAMRMultigrid.query("multigrid", opt.AMRMultigrid_verbosity);
//  ppAMRMultigrid.query("relaxMode", opt.AMRMultigridRelaxMode);
  ppAMRMultigrid.query("hang_eps", opt.AMRMultigridHang);
  ppAMRMultigrid.query("tolerance", opt.AMRMultigridTolerance);
  ppAMRMultigrid.query("norm_thresh", opt.AMRMultigridNormThresh);
  ppAMRMultigrid.query("verbosity", opt.AMRMultigridVerb);
  ppAMRMultigrid.query("MGtype", mgtype);
  opt.MGtype = MGmethod(mgtype);


//  opt.implicitAdvectionSolve = false;
//  ppMain.query("implicitAdvectionSolve", opt.implicitAdvectionSolve);

//  opt.usePhiForImplicitAdvectionSolve = true;
//  ppMain.query("usePhiForImplicitAdvectionSolve", opt.usePhiForImplicitAdvectionSolve);

  opt.maxDivUFace = 1e-8;
  ppMain.query("maxDivUFace", opt.maxDivUFace);

  opt.multiCompUStarSolve = false;
  ppMain.query("multiCompUStarSolve", opt.multiCompUStarSolve);

  opt.nonlinearHCOpSuperOptimised = false;
  ppMain.query("nonlinearHCOpSuperOptimised", opt.nonlinearHCOpSuperOptimised);

//  opt.smoothingCoeff = 0.0; //0.01;
//  ppProjection.query("pre_smoothing", opt.smoothingCoeff);

  opt.velMGNumSmooth=2;
  opt.velMGTolerance=1e-10;
  opt.velMGHang=1e-10;
  opt.velMGNormThresh=1e-10;
  opt.velMGNumMG=1;
  opt.VelMGMaxIter=10;

  ppVelMultigrid.query("num_smooth", opt.velMGNumSmooth);
  ppVelMultigrid.query("tolerance",  opt.velMGTolerance);
  ppVelMultigrid.query("hang_eps",   opt.velMGHang);
  ppVelMultigrid.query("num_mg",     opt.velMGNumMG);
  ppVelMultigrid.query("norm_thresh", opt.velMGNormThresh);
  ppVelMultigrid.query("max_iter",  opt.VelMGMaxIter);

  opt.HCMultigridNumSmoothUp=4;
  opt.HCMultigridNumSmoothDown=1;
  opt.HCMultigridNumMG=1;
  opt.HCMultigridMaxIter=10;
  opt.HCMultigridVerbosity=0;
  opt.HCMultigridBottomSolveIterations=40;
  opt.HCMultigridTolerance=1e-10;
  opt.HCMultigridHang=1e-10;
  opt.HCMultigridNormThresh=1e-10;
  opt.HCMultigridRelaxMode = 1; // 1=GSRB, 4=jacobi
  opt.HCMultigridUseRelaxBottomSolverForHC = true;

  HCMultigrid.query("num_smooth_up", opt.HCMultigridNumSmoothUp);
  HCMultigrid.query("num_mg", opt.HCMultigridNumMG);
  HCMultigrid.query("hang_eps", opt.HCMultigridHang);
  HCMultigrid.query("norm_thresh", opt.HCMultigridNormThresh);
  HCMultigrid.query("tolerance", opt.HCMultigridTolerance);
  HCMultigrid.query("max_iter", opt.HCMultigridMaxIter);
  HCMultigrid.query("verbosity", opt.HCMultigridVerbosity);
  HCMultigrid.query("numSmoothDown", opt.HCMultigridNumSmoothDown);
  HCMultigrid.query("relaxMode", opt.HCMultigridRelaxMode);
  HCMultigrid.query("bottomSolveIterations", opt.HCMultigridBottomSolveIterations);
  HCMultigrid.query("useRelaxBottomSolver", opt.HCMultigridUseRelaxBottomSolverForHC);

//  opt.noMultigrid = false;
//  opt.noMultigridIter = 100;
//  ppMain.query("noMultigrid", opt.noMultigrid);
//  ppMain.query("noMultigridIter", opt.noMultigridIter);

  // 1 -> PLM, 2 -> PPM
  opt.velAdvNormalPredOrder = 1;

  // Use 4th order slope computations
  opt.velAdvUseFourthOrderSlopes = true;


  // No artificial viscosity
  opt.velAdvUseArtVisc = false;
  opt.velAdvArtVisc = -0.0;

  ppPatchGodunov.query("velOrder", opt.velAdvNormalPredOrder);
  ppPatchGodunov.query("velFourthOrderSlopes", opt.velAdvUseFourthOrderSlopes);
  ppPatchGodunov.query("velUseArtVisc", opt.velAdvUseArtVisc);
  ppPatchGodunov.query("velArtVisc", opt.velAdvArtVisc);

  // There is a bug in Chombo at the moment which means we have to use the limiter if we're
  // doing 2nd order slopes
  if (opt.velAdvNormalPredOrder == 2)
  {
    opt.velAdvHigherOrderLimiter = true;
  }
  else
  {
    opt.velAdvHigherOrderLimiter = false;
  }
  ppPatchGodunov.query("higherOrderLimiter", opt.velAdvHigherOrderLimiter);

  opt.HCUseArtVisc = opt.velAdvUseArtVisc;
  opt.HCArtVisc = opt.velAdvArtVisc;

  // This speeds up convergence a little.
  // The eutectic boundary causes issues with higher order methods.
  opt.HCNormalPredOrder = 1;
  opt.HCUseFourthOrderSlopes = false;

  opt.HCHigherOrderLimiter = opt.velAdvHigherOrderLimiter;

  ppPatchGodunov.query("HCOrder", opt.HCNormalPredOrder);
  ppPatchGodunov.query("HCFourthOrderSlopes", opt.HCUseFourthOrderSlopes);
  ppPatchGodunov.query("HCUseArtVisc", opt.HCUseArtVisc);
  ppPatchGodunov.query("HCArtVisc", opt.HCArtVisc);

  /***
   * Initialisation
   */

  opt.num_init_passes = 1;
  ppMain.query("num_init_passes", opt.num_init_passes);

  opt.restart_new_time = -1.0;
  ppMain.query("restart_newTime", opt.restart_new_time);

//  opt.increaseDt = true;
//  ppMain.query("init_increase_dt", opt.increaseDt);

  opt.init_use_prev_pressure_for_Ustar = false;
  ppMain.query("init_add_subtract_grad_p", opt.init_use_prev_pressure_for_Ustar); // legacy option
  ppMain.query("init_use_prev_pressure_for_Ustar", opt.init_use_prev_pressure_for_Ustar);

  opt.init_compute_uDelu = true;
  ppMain.query("init_compute_uDelu", opt.init_compute_uDelu);

  // Restarting
  // Take care of changing dimensionless units if necessary
  // Default reference temperature/salinity
  opt.refTemp = -23;
  opt.refSalinity = 230;

  // Check to see if we've specified other reference values
  ppParams.query("referenceTemperature", opt.refTemp);
  ppParams.query("referenceSalinity", opt.refSalinity);

  // Default reference values
  opt.prevRefTemp = opt.refTemp;
  opt.prevRefSalinity = opt.refSalinity;

  // In case we've specified other reference values
  ppParams.query("prevReferenceTemperature", opt.prevRefTemp);
  ppParams.query("prevReferenceSalinity", opt.prevRefSalinity);

  opt.dtReductionFactor = 1.0;
  ppMain.query("restart_dtReduction", opt.dtReductionFactor);

  /**
   * Initial conditions options
   */

  opt.enforceAnalyticSoln = false;
  ppMain.query("enforceAnalyticSoln", opt.enforceAnalyticSoln);
  int analyticSoln = -1;
  ppMain.query("analyticSoln", analyticSoln);
  if (analyticSoln > -1)
  {
    opt.enforceAnalyticSoln = true;
  }

  // Actual value of the analytic soln is by default the same as this physical problem
  opt.analyticSolution = params.physicalProblem;
  ppMain.query("analyticSoln", opt.analyticSolution);

//  opt.useAnalyticSource = false;
//  ppMain.query("analyticSourceTerm", opt.useAnalyticSource);

  opt.initialPerturbation = 0.0;
  ppMain.query("initialPerturbation", opt.initialPerturbation);

  opt.perturbationPhaseShift = 0.0;
  ppMain.query("perturbationPhaseShift", opt.perturbationPhaseShift);

  opt.delayedPerturbation = 0.0;
  ppMain.query("delayedPerturbaation", opt.delayedPerturbation);

  opt.perturbationTime = 0.0;
  ppMain.query("perturbationTime", opt.perturbationTime);

  opt.perturbationWavenumber = 0.0;
  ppMain.query("perturbationWavenumber", opt.perturbationWavenumber);

  opt.perturbationSin = false;
  ppMain.query("perturbationSin", opt.perturbationSin);

  opt.maxRestartWavenumbers = 50;
  ppMain.query("maxRestartWavenumber", opt.maxRestartWavenumbers);

  opt.fixedPorosity = -1.0;
  ppMain.query("fixed_porosity", opt.fixedPorosity);

  opt.porosityTimescale = 1/params.darcy;
  ppMain.query("porosityTimescale", opt.porosityTimescale);

  opt.initVel = 0.0;
  ppInit.query("initVel", opt.initVel);

  int porosityFunc = PorosityFunctions::constantLiquid;
  ppMain.query("porosity_function", porosityFunc);
  opt.porosityFunction = PorosityFunctions(porosityFunc);

  opt.restartPerturbation = 0.0;
  ppMain.query("restart_perturbation", opt.restartPerturbation);

  opt.porousHoleRadius = 0.1*opt.domainWidth;
  ppMain.query("radius", opt.porousHoleRadius);

  opt.initVelScale = 1e-3;
  ppMain.query("initVelScale", opt.initVelScale);

  opt.summerProfile = -1;
  ppInit.query("summerProfile", opt.summerProfile);

  opt.mushHeight = opt.domainHeight/2;
  ppInit.query("summerProfileMushHeight", opt.mushHeight);

  opt.linearGradient = false;
  ppInit.query("linearGradient", opt.linearGradient);

  opt.meltPondDepth = 0;
  ppMain.query("meltPondDepth", opt.meltPondDepth);

  // 0 is enthalpy
  opt.restart_perturbation_var = ScalarVars::m_enthalpy;
  ppMain.query("restart_perturbation_var", opt.restart_perturbation_var);

  opt.horizAverageRestart = false;
  ppMain.query("horizontallyAverageRestart", opt.horizAverageRestart);

//  ppMain.query("doAutomaticRestart",opt.initiallyDoAutomaticRestart);


  /**
   * Other options
   */

  opt.FixedPorositySTD = 0.005;
  opt.fixedPorosityMaxChi = 1.05; // want to make this a bit more than 1, so we get porosity=1 region of finite size
  ppMain.query("maxChi", opt.fixedPorosityMaxChi);
  ppMain.query("stdev", opt.FixedPorositySTD);

  opt.fixedPorosityFractionalInnerRadius = 0.2;
  ppMain.query("innerRadius", opt.fixedPorosityFractionalInnerRadius);

  opt.fixedPorosityEndTime = -1;
  ppMain.query("porosityEndTime", opt.fixedPorosityEndTime);

  // Try a corner copier?
  opt.scalarExchangeCorners = true;
  ppMain.query("scalarExchangeCorners", opt.scalarExchangeCorners);

  opt.buoyancy_zero_time = -1;
  ppMain.query("turn_off_buoyancy_time", opt.buoyancy_zero_time);



//  opt.maxEta = -1;
//  ppMain.query("max_eta", opt.maxEta); // let user specify different max eta if they want

  opt.minDt = 1e-7;
  ppMain.query("min_dt", opt.minDt);

  opt.max_dt = -1;
  ppMain.query("max_dt", opt.max_dt);

  opt.printAccelDt = false;
  ppMain.query("printAccelDt", opt.printAccelDt);

  opt.useAccelDt = false;
  ppMain.query("useAccelDt", opt.useAccelDt);

  opt.useInitAccelDt = true;
  ppMain.query("useInitAccelDt", opt.useInitAccelDt);

  opt.accelCFL = -1;
  ppMain.query("accelCFL", opt.accelCFL);

  opt.max_init_dt = -1;
  ppMain.query("max_init_dt", opt.max_init_dt);

  opt.max_possible_level = 0;
  ppMain.query("max_level", opt.max_possible_level);

  opt.computeVorticityStreamFunction = true;
  ppMain.query("computeVorticity", opt.computeVorticityStreamFunction);


  opt.useFortranRegularisationFace = false;
  ppMain.query("fortranRegularisationFace", opt.useFortranRegularisationFace);
  opt.useFortranRegularisation = true;
  ppMain.query("fortranRegularisation", opt.useFortranRegularisation);


//  opt.stokesDarcyForcingTimescale = 0.5;
//  ppParams.query("forcing_timescale", opt.stokesDarcyForcingTimescale);

  opt.viscousBCs = params.isViscous();
    ppMain.query("viscousBCs", opt.viscousBCs);



  a_fact = RefCountedPtr<AMRLevelMushyLayerFactory> (new AMRLevelMushyLayerFactory(opt, params));

}
void
defineAMR(AMR&                                          a_amr,
          RefCountedPtr<AMRLevelMushyLayerFactory>&     a_fact,
          const ProblemDomain&                          a_prob_domain,
          const Vector<int>&                            a_refRat)
{
//  pout() << "MushyLayerSubcycleUtils - defineAMR(..) - creating AMR object" << endl;

  ParmParse ppMain("main");
  int max_level = 0;
  ppMain.get("max_level",max_level);

  int num_read_levels = Max(max_level,1);
  std::vector<int> regrid_intervals;
  if (max_level > 0)
  {
    ppMain.getarr("regrid_interval",regrid_intervals,0,num_read_levels);
  }
  else
  {
    regrid_intervals.push_back(-1);
  }

  if (max_level > 0)
  {
    std::vector<int> ref_ratios = std::vector<int>(); // (num_read_levels,1);
    ppMain.getarr("ref_ratio", ref_ratios, 0, num_read_levels);
    for (int i=0; i<=max_level; i++)
    {
      if (ref_ratios[i] == 1)
      {
        MayDay::Error("Refinement ratio can't be 1!");
      }
    }
  }

  Vector<Vector<Box> > fixedGrids;
  bool predefinedGrids = false;
  predefinedGrids = getFixedGrids(fixedGrids, a_prob_domain);

  if (predefinedGrids &&
      fixedGrids.size() != max_level + 1)
  {
    MayDay::Error("Specified grids do not contain the correct number of levels");
  }

  if (predefinedGrids)
  {
    for (int i = 0; i <=max_level; i++)
    {
      regrid_intervals[i] = -1;
    }
  }

  int block_factor = 1;
  ppMain.query("block_factor",block_factor);

  int max_grid_size = 128;
  ppMain.query("max_grid_size",max_grid_size);

  Real fill_ratio = 0.75;
  ppMain.query("fill_ratio",fill_ratio);

  int checkpoint_interval = 0;
  ppMain.query("checkpoint_interval",checkpoint_interval);

  int plot_interval = 0;
  ppMain.query("plot_interval",plot_interval);

  Real plot_period = 0;
  ppMain.query("plot_period", plot_period);

  Real max_dt_growth = 1.1;
  ppMain.query("max_dt_growth",max_dt_growth);

  Real fixed_dt = -1.0;
  ppMain.query("fixed_dt", fixed_dt);

  Real dt_tolerance_factor = 1.1;
  ppMain.query("dt_tolerance_factor",dt_tolerance_factor);
  AMR amr;
  a_amr.define(max_level, a_refRat,
               a_prob_domain,&(*a_fact));


  // set grid generation parameters
  a_amr.maxGridSize(max_grid_size);
  a_amr.blockFactor(block_factor);
  a_amr.fillRatio(fill_ratio);

  // the hyperbolic codes use a grid buffer of 1
  int gridBufferSize = 1;
  ppMain.query("grid_buffer_size",gridBufferSize);
  a_amr.gridBufferSize(gridBufferSize);

  // set output parameters
  a_amr.checkpointInterval(checkpoint_interval);
  a_amr.plotInterval(plot_interval);
  a_amr.plotPeriod(plot_period);
  a_amr.regridIntervals(regrid_intervals);
  a_amr.maxDtGrow(max_dt_growth);
  a_amr.dtToleranceFactor(dt_tolerance_factor);
  // We print diagnostics in the routine which checks for steady state.
  // Therefore must ensure this is true, or won't print diagnostics.
  a_amr.checkForSteadyState(true);

  std::string output_folder = "";
  ppMain.query("output_folder", output_folder);


  if (fixed_dt > 0)
  {
    a_amr.fixedDt(fixed_dt);
  }
  if (ppMain.contains("use_subcycling"))
  {
    bool useSubcycling;
    ppMain.get("use_subcycling", useSubcycling);
    if (!useSubcycling)
    {
      pout() << "SUBCYCLING IN TIME TURNED OFF!!!"  << endl;
    }
    a_amr.useSubcyclingInTime(useSubcycling);
  }
  if (ppMain.contains("plot_prefix"))
  {
    std::string prefix;
    ppMain.get("plot_prefix",prefix);

    // I think this length being too small is the cause of stack smashing
    char(fullPrefix[2000]);
    sprintf(fullPrefix, "%s/%s", output_folder.c_str(),prefix.c_str());

    a_amr.plotPrefix(fullPrefix);
  }

  if (ppMain.contains("chk_prefix"))
  {
    std::string prefix;
    ppMain.get("chk_prefix",prefix);

    char(fullPrefix[2000]);
    sprintf(fullPrefix, "%s/%s", output_folder.c_str(),prefix.c_str());

    a_amr.checkpointPrefix(fullPrefix);
  }

#ifdef CH_FORK
  if (ppMain.contains("subcycled_plots"))
  {
    bool subcycledPlots = false;
    ppMain.get("subcycled_plots", subcycledPlots);
    a_amr.subcycledPlots(subcycledPlots);
  }

  if (ppMain.contains("chombo_chk_not_plot"))
  {
    bool makeCheck = false;
    ppMain.get("chombo_chk_not_plot", makeCheck);
    a_amr.makeChkNotPlot(makeCheck);
  }
#endif

  int verbosity = 1;
  ppMain.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  a_amr.verbosity(verbosity);


}

bool
getFixedGrids(Vector<Vector<Box> >& amrGrids,  ProblemDomain prob_domain, string gridfileParam)
{

  //  pout() << "getFixedGrids" << endl;

  ParmParse ppMain("main");
  int verbosity = 3;

  int max_level;
  int max_grid_size;

  ppMain.get("max_level", max_level);
  ppMain.get("max_grid_size", max_grid_size);


  bool predefinedGrids = ppMain.contains(gridfileParam);
  string gridfile;

  // Don't get grids if we only have one level
  if (predefinedGrids && max_level > 0)
  {
    ppMain.get(gridfileParam.c_str(), gridfile);
  }
  else
  {
    return false;
  }

#ifdef CH_MPI
  // if (procID() ==  uniqueProc(SerialTask::compute))
  MPI_Barrier(Chombo_MPI::comm);
  if (procID() ==  0)
  {
#endif
    // read in predefined grids
    ifstream is(gridfile.c_str(), ios::in);
    if (is.fail())
    {
      pout() << "cannot open grids file " << gridfile << endl;
      MayDay::Error("Cannot open grids file");
    }

#ifdef CH_MPI

    pout() << "procID: " << procID() << "  opening gridfile \n" << endl;

#endif

    // format of file -- number of levels, then for
    // each level starting with level 1, number of
    // grids on level, list of boxes
    int in_numLevels;
    is >> in_numLevels;
    CH_assert (in_numLevels <= max_level+1);
    if (verbosity >= 3)
    {
      pout() << "numLevels = " << in_numLevels << endl;
    }
    while (is.get() != '\n');
    amrGrids.resize(in_numLevels);
    // check to see if coarsest level needs to be broken up
    domainSplit(prob_domain,amrGrids[0], max_grid_size,4);

    if (verbosity >= 3)
    {
      pout() << "level 0: ";
      for (int n=0; n<amrGrids[0].size(); n++)
      {
        pout () << amrGrids[0][n] << endl;
      }
    }
    // now loop over levels, starting with level 1
    int ngrid;
    for (int lev=1; lev<in_numLevels; lev++)
    {
      is >> ngrid;
      if (verbosity >= 3)
      {
        pout() << "level " << lev << " numGrids = " << ngrid << endl;
        pout() << "Grids: ";
      }
      while (is.get() != '\n');
      amrGrids[lev].resize(ngrid);
      for (int i=0; i<ngrid; i++)
      {
        Box bx;
        is >> bx;

        // advance to next box
        while (char ch = is.get())
        {
          if (ch == '#') break;
          if (ch == '\n') break;
        }

        // quick check on box size
        Box bxRef(bx);
        // not really sure why i was doing this (holdover from
        // legacy code)
        //bxRef.refine(ref_ratios[lev-1]);
        if (bxRef.longside() > max_grid_size)
        {
          pout() << "Grid " << bx << " too large" << endl;
          MayDay::Error();
        }
        if (verbosity >= 3)
        {
          pout() << bx << endl;
        }
        amrGrids[lev][i] = bx;
      } // end loop over boxes on this level
    } // end loop over levels

#ifdef CH_MPI
  } // end if procID = 0
  //broadcast (amrGrids, uniqueProc(SerialTask::compute));
  broadcast (amrGrids, 0);
#endif



  return true;
}

void
setupAMRForAMRRun(AMR& a_amr, ProblemDomain prob_domain)
{

  ParmParse ppMain("main");

  // make new blank diagnostics file
  std::ofstream diagnosticsFile ("diagnostics.out");
  diagnosticsFile.close();

  // Check
  Vector<Vector<Box> > fixedGrids;
  bool predefinedGrids = false;
  predefinedGrids = getFixedGrids(fixedGrids, prob_domain);

  if (ppMain.contains("restart_file"))
  {
    std::string restart_file;
    ppMain.get("restart_file",restart_file);
    pout() << " restarting from file " << restart_file << endl;

#ifdef CH_USE_HDF5
    HDF5Handle handle(restart_file,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    pout() << "Opened checkpoint file" << endl;
    a_amr.setupForRestart(handle);
    handle.close();
#else

    MayDay::Error("restart only defined with hdf5");
#endif

    Real resetTime = -1;
    ppMain.query("restart_newTime", resetTime);
    if (resetTime >= 0 )
    {
#ifdef CH_FORK
      a_amr.cur_time(resetTime);
      Vector<AMRLevel*> amrLev = a_amr.getAMRLevels();
      for (int i=0; i < amrLev.size(); i++)
      {
        amrLev[i]->time(resetTime);

      }
      pout() << "Set time = " << resetTime << endl;
#else
      MayDay::Warning("Unable to reset time as you're not using the forked version of Chombo, and therefore your AMR class doesn't have a cur_time() method");
#endif
    }


  }
  else if (predefinedGrids)
  {
    a_amr.setupForFixedHierarchyRun(fixedGrids, 1);
  }
  else
  {
    // initialize from scratch for AMR run
    // initialize hierarchy of levels
    a_amr.setupForNewAMRRun();
  }




}




#include "NamespaceFooter.H"
