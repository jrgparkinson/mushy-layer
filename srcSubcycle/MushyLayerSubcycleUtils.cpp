
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

  ParmParse ppMain("main");

  ParmParse ppParams("parameters");
  ParmParse ppMG("amrmultigrid");
  ParmParse ppProjection("projection");
  ParmParse ppInit("init");
  ParmParse ppAdvsrc("advSrc");

  MushyLayerOptions opt;

  opt.cfl = 0.8;
  ppMain.get("cfl",opt.cfl);

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

  opt.advVelCentering = 0.1;
  ppMain.query("init_advVel_centering", opt.advVelCentering);

  opt.initial_dt_multiplier = 0.1;
  ppMain.query("initial_cfl", opt.initial_dt_multiplier);


  opt.computeDiagnostics = true;
  ppMain.query("computeDiagnostics", opt.computeDiagnostics);

  opt.doEulerPart = true;
  ppMain.query("doEuler", opt.doEulerPart);



  opt.doScalarAdvectionDiffusion = true;
  ppMain.query("doScalarAdvectionDiffusion", opt.doScalarAdvectionDiffusion);

  int advectionMethod =  velocityAdvectionTypes::m_porosityOutsideAdvection;
  ppMain.query("advectionMethod", advectionMethod);
  opt.advectionMethod = velocityAdvectionTypes(advectionMethod);

  opt.skipTrickySourceTerm = -1;
  ppMain.query("skipTrickySourceTime", opt.skipTrickySourceTerm);

  opt.allowMulticompAdvection = true;
  ppMain.query("allowMulticompAdvection", opt.allowMulticompAdvection);


  opt.scaleP_MAC = true;
  ppMain.query("scalePwithChi", opt.scaleP_MAC);
  // Makes more sense to store this parameter under projection
  ppProjection.query("scalePressureWithPorosity", opt.scaleP_MAC);

  opt.scaleP_CC = opt.scaleP_MAC;
  ppProjection.query("scaleCCPressure", opt.scaleP_CC);

  opt.explicitDarcyTerm = false;
  ppMain.query("explicitDarcyTerm", opt.explicitDarcyTerm);

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

  // by default make this tiny (so essentially turned off)
  opt.uDelU_porosityLimit = 10*opt.lowerPorosityLimit; //1e-15;
  ppMain.query("uDelu_porosity", opt.uDelU_porosityLimit);

  opt.lapVelBCOrder = 0;
  ppMain.query("lap_vel_bc_order", opt.lapVelBCOrder);

  // Only do multiple projections on the base level
  opt.maxProjBaseLevel = 1;
  ppMain.query("maxProjections", opt.maxProjBaseLevel);

  opt.enforceAnalyticVel = false;
  ppMain.query("analyticVel", opt.enforceAnalyticVel);

  opt.useOldAdvVel = false;
  ppMain.query("useOldAdvVelForTracing", opt.useOldAdvVel);

  opt.lapVelNumSmooth = 0;
  ppMain.query("lapVelNumSmooth", opt.lapVelNumSmooth);

  opt.lapVelSmoothScale = 0.0;
  ppMain.query("lapVelSmoothScale", opt.lapVelSmoothScale);

  opt.CCVelSrcTermCentering = 0.5;
  ppMain.query("vel_src_centering", opt.CCVelSrcTermCentering);

  opt.advSrcAllowLaggedLapVel = false;
  ppAdvsrc.query("allow_lagged_lap_vel", opt.advSrcAllowLaggedLapVel);


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

  // Cell-centred velocity source term
  ParmParse ppCCSrc("ccSrc");

  opt.CCAdvSrc = true;
  ppCCSrc.query("advection", opt.CCAdvSrc);

  // This line must come after quering explicit Darcy Term
  opt.CCDarcySrc = opt.explicitDarcyTerm;
  ppCCSrc.query("darcy", opt.CCDarcySrc);

  opt.CCBuoyancySrc = true;
  ppCCSrc.query("buoyancy", opt.CCBuoyancySrc);



  /**
   * Physics related options
   */

  opt.iter_plot_interval = -1;
  ppMain.query("iter_plot_interval", opt.iter_plot_interval);

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
  ppMain.get ("refine_thresh", opt.refineThresh);

  opt.tagBufferSize = 4;
  ppMain.get("tag_buffer_size", opt.tagBufferSize);

  opt.doProjection=true;
  ppMain.query("doProjection", opt.doProjection);

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

  opt.variable_eta_factor = 1.0;
  ppMain.query("variable_eta_factor", opt.variable_eta_factor);
  CH_assert(opt.variable_eta_factor >= 1); // This must be >= 1, else eta will increase when it should be decreasing (and vice-versa)

  opt.compute_initial_VD_corr = true;
  ppMain.query("initialize_VD_corr", opt.compute_initial_VD_corr);

  /**
   * Solver options
   */
  opt.useLimiting = true;
  ppMain.get("use_limiting", opt.useLimiting);

  opt.viscous_solver_tol = 1e-10;
  ppMain.query("viscous_solver_tol", opt.viscous_solver_tol);

  opt.viscous_num_smooth_down = 8;
  ppMain.query("viscous_num_smooth_down", opt.viscous_num_smooth_down);

  opt.viscous_num_smooth_up = 8;
  ppMain.query("viscous_num_smooth_up", opt.viscous_num_smooth_up);


  opt.timeIntegrationOrder = 1;
  ppMain.query("time_integration_order", opt.timeIntegrationOrder);

  opt.verbosity_multigrid = 0;
  ppMG.query("multigrid", opt.verbosity_multigrid);

  opt.implicitAdvectionSolve = false;
  ppMain.query("implicitAdvectionSolve", opt.implicitAdvectionSolve);

  opt.maxDivUFace = 1e-8;
  ppMain.query("maxDivUFace", opt.maxDivUFace);

  int mgtype = MGmethod::MGTypeFAS;
  ppMG.query("MGtype", mgtype);
  opt.MGtype = MGmethod(mgtype);

  opt.multiCompUStarSolve = false;
  ppMain.query("multiCompUStarSolve", opt.multiCompUStarSolve);

  opt.nonlinearHCOpSuperOptimised = false;
  ppMain.query("nonlinearHCOpSuperOptimised", opt.nonlinearHCOpSuperOptimised);

  opt.smoothingCoeff = 0.0; //0.01;
  ppProjection.query("pre_smoothing", opt.smoothingCoeff);

  /***
   * Initialisation
   */

  opt.num_init_passes = 1;
  ppMain.query("num_init_passes", opt.num_init_passes);

  opt.restart_new_time = -1.0;
  ppMain.query("restart_newTime", opt.restart_new_time);

  opt.increaseDt = true;
  ppMain.query("init_increase_dt", opt.increaseDt);

  opt.init_add_subtract_grad_p = false;
  ppMain.query("init_add_subtract_grad_p", opt.init_add_subtract_grad_p);

  opt.init_compute_uDelu = true;
  ppMain.query("init_compute_uDelu", opt.init_compute_uDelu);

  /**
   * Initial conditions options
   */

  opt.enforceAnalyticSoln = false;
  ppMain.query("enforceAnalyticSoln", opt.enforceAnalyticSoln);
  Real analyticSoln=-1;
  ppMain.query("analyticSoln", analyticSoln);
  if (analyticSoln > -1)
  {
    opt.enforceAnalyticSoln = true;
  }

  opt.initialPerturbation = 0.0;
  ppMain.query("initialPerturbation", opt.initialPerturbation);


  opt.delayedPerturbation = 0.0;
  ppMain.query("delayedPerturbaation", opt.delayedPerturbation);

  opt.perturbationTime = 0.0;
  ppMain.query("perturbationTime", opt.perturbationTime);

  opt.perturbationWavenumber = 1.0;
  ppMain.query("perturbationWavenumber", opt.perturbationWavenumber);

  opt.perturbationSin = false;
  ppMain.query("perturbationSin", opt.perturbationSin);

  opt.fixedPorosity = -1.0;
  ppMain.query("fixed_porosity", opt.fixedPorosity);

  opt.initVel = 0.0;
  ppInit.query("initVel", opt.initVel);

  int porosityFunc = PorosityFunctions::constantLiquid;
  ppMain.query("porosity_function", porosityFunc);
  opt.porosityFunction = PorosityFunctions(porosityFunc);

  opt.restartPerturbation = 0.0;
  ppMain.query("restart_perturbation", opt.restartPerturbation);

  opt.radius = 0.1*opt.domainWidth;
  ppMain.query("radius", opt.radius);

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

  a_fact = RefCountedPtr<AMRLevelMushyLayerFactory> (new AMRLevelMushyLayerFactory(opt));

}
void
defineAMR(AMR&                                          a_amr,
          RefCountedPtr<AMRLevelMushyLayerFactory>&  a_fact,
          const ProblemDomain&                          a_prob_domain,
          const Vector<int>&                            a_refRat)
{
  pout() << "MushyLayerSubscyleUtils - defineAMR(..) - creating AMR object" << endl;

  ParmParse ppMain("main");
  int max_level = 0;
  ppMain.get("max_level",max_level);

  int num_read_levels = Max(max_level,1);
  std::vector<int> regrid_intervals; // (num_read_levels,1);
  ppMain.getarr("regrid_interval",regrid_intervals,0,num_read_levels);

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
  ppMain.get("block_factor",block_factor);

  int max_grid_size = 32;
  ppMain.get("max_grid_size",max_grid_size);

  Real fill_ratio = 0.75;
  ppMain.get("fill_ratio",fill_ratio);

  int checkpoint_interval = 0;
  ppMain.get("checkpoint_interval",checkpoint_interval);

  int plot_interval = 0;
  ppMain.query("plot_interval",plot_interval);

  Real plot_period = 0;
  ppMain.query("plot_period", plot_period);

  Real max_dt_growth = 1.1;
  ppMain.get("max_dt_growth",max_dt_growth);

  Real fixed_dt = -1.0;
  ppMain.query("fixed_dt", fixed_dt);

  Real dt_tolerance_factor = 1.1;
  ppMain.get("dt_tolerance_factor",dt_tolerance_factor);
  AMR amr;
  a_amr.define(max_level, a_refRat,
               a_prob_domain,&(*a_fact));



  // set grid generation parameters
  a_amr.maxGridSize(max_grid_size);
  a_amr.blockFactor(block_factor);
  a_amr.fillRatio(fill_ratio);

  // the hyperbolic codes use a grid buffer of 1
  int gridBufferSize;
  ppMain.get("grid_buffer_size",gridBufferSize);
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

  int verbosity;
  ppMain.get("verbosity",verbosity);
  CH_assert(verbosity >= 0);


  a_amr.verbosity(verbosity);

  pout() << "MushyLayerSubscyleUtils - defineAMR(..) - finished setting up AMR" << endl;
  //	 pout() << "AMR - set verbosity" << endl;

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
  //  pout() << "setupAMRForAMRRun" << endl;

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
