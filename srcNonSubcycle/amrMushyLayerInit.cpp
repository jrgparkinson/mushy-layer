#include "amrMushyLayer.H"

#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "CoarseAverage.H"
#include "ParmParse.H"
#include "CellToEdge.H"
#include "AdvectIBC.H"
#include "LevelAdvect.H"
#include "ExtrapFillPatch.H"
#include "FineInterp.H"
//#include "LevelFluxRegisterEdge.H"


void
amrMushyLayer::setDefaults()
{
  // set some bogus values as defaults
  m_is_defined = false;
  m_max_level = -1;
  m_finest_level = -1;
  m_block_factor = -1;
  m_fill_ratio = -1;
  m_do_restart = false;
  m_restart_step = -1;

  // set the rest of these to reasonable defaults
  m_nesting_radius = 1;
  m_tagging_val = 0.1;
  m_tags_grow = 1;
  m_cfl = 0.25;
  m_max_dt_grow = 1.5;
  m_max_dt = 0.001;
  m_dt = 1e-2;
  m_time = 0.0;
  m_cur_step = 0;
  m_max_box_size = 64;
  m_refine_initial_domain = false;
  // domain size default is unit square/cube
  m_domainSize = RealVect(D_DECL(1.0,1.0,1.0));
  m_steady_state = false;
  //todo make this 4 unless periodic
  m_num_ghost = 4; // Warning - think this wants to be 4, or else stuff breaks
  m_ghostVect = m_num_ghost*IntVect::Unit;
  m_parameters.getParameters();
 // m_parameters.printParameters();
  m_timeIntegrationOrder = 1;
  m_enforceAnalyticSoln = true;
  m_printAnalyticSoln = false;
  m_frameAdvectionMethod = m_godunov;
  m_minTimestep = 0;
  m_fluidVelDiffOrder = 2;
  m_fluidVelInterpOrder = 1;
  m_ownVelCalc = true;
  m_tagging_scalar_var = 5;
  m_tagging_vector_var = -1;

  m_debugging_resetTemperatureField = false;
  m_debugging_fixVelocity = false;

  m_output_folder = "output";
  m_plot_prefix = "plt";
  m_plot_interval = 10000000;
  m_check_prefix = "chk";
  m_check_interval = -1;
  m_iteration_plot_interval = -1;

  m_total_num_iterations = 1;


  m_physBCPtr = NULL; //new PhysBCUtil(m_parameters, m_amrDx[0]);

  // MAKE SURE YOU UPDATE THIS IF YOU ADD ANOTHER (SCALAR) VARIABLE!
  //  m_numVars = 48;

  m_varNames = Vector<string>(m_numVars, string("Placeholder name"));


  m_varNames[m_enthalpy] = string("Enthalpy");
  m_varNames[m_enthalpySolidus] = string("Enthalpy solidus");
  m_varNames[m_enthalpyLiquidus] = string("Enthalpy liquidus");
  m_varNames[m_enthalpyEutectic] = string("Enthalpy eutectic");
  m_varNames[m_bulkConcentration] = string("Composition");
  m_varNames[m_theta] = string("temperature");
  m_varNames[m_thetaLiquidus] = string("temperature liquidus");
  m_varNames[m_thetaSolidus] = string("temperature solidus");
  m_varNames[m_porosity] = string("Porosity");
  m_varNames[m_compositionLiquid] = string("Liquid Composition");
  m_varNames[m_compositionSolid] = string("Solid Composition");
  m_varNames[m_porosityEutectic] = string("Eutectic Porosity");
  m_varNames[m_thetaForcing] = string("theta forcing");
  m_varNames[m_liquidCompositionGrad] = string("Liquid Composition Gradient");
  m_varNames[m_solidFraction] = string("Solid Fraction");
  m_varNames[m_steadyStateImbalance] = string("Steady state imbalance");
  m_varNames[m_Tanalytic] = string("Analytic Temperature");
  m_varNames[m_solidFractionTrue] = string("Analytic Solid Fraction");
  m_varNames[m_thetaTrue] = string("Analytic theta");
  m_varNames[m_enthalpyAdvection] = string("Enthalpy advection");
  m_varNames[m_thetaLaplacian] = string("theta laplacian");
  m_varNames[m_streamFunction] = string("Stream function");
  m_varNames[m_ThetaLAnalytic] = string("Analytic liquid concentration");
  m_varNames[m_resid] = string("Residual");
  m_varNames[m_concSource] = string("Concentration Eqn Source");
  m_varNames[m_ThetaDiffusion] = string("Theta Diffusion");
  m_varNames[m_ThetaDiffusionN] = string("Theta Diffusion at n");
  m_varNames[m_thetaBcoef] = string("theta B coeff");
  m_varNames[m_permeability] = string("Permeability");
  m_varNames[m_pressure] = string("Pressure");
  m_varNames[m_ThetaBCoef] = string("Theta B coeff");
  m_varNames[m_divU] = string("Divergence U");
  m_varNames[m_enthalpyAnalytic] = string("Enthalpy analytic");
  m_varNames[m_ThetaAnalytic] = string("Theta analytic");
  m_varNames[m_ThetaFrameAdvection] = string("Theta frame adv");
  m_varNames[m_ThetaSSource] = string("Theta solid source");
  m_varNames[m_ThetaPorositySource] = string("Theta porosity source");
  m_varNames[m_pressureError] = string("Pressure Error");
  m_varNames[m_pressureAnalytic] = string("Pressure analytic");
  m_varNames[m_divUerror] = string("Div U rror");
  m_varNames[m_thetaForcingError] = string("theta forcing error");
  m_varNames[m_thetaForcingAnalytic] = string("theta forcing analytic");
  m_varNames[m_CFInterpTest] = string("CF interp test");
  m_varNames[m_divUstar] = string("div U star");
  m_varNames[m_divUstarErr] = string("div U star err");
  m_varNames[m_lambda] = string("lambda passive tracer");
  m_varNames[m_lambdaPostCorr] = string("lambda after correction");
  m_varNames[m_buoyancyFilter] = string("buoyancy Filter");
  //m_varNames[m_] = string("");


  //If adding to this list, make sure you update m_numVars above


  m_numVectorVars = 6;

  m_vectorVarNames = Vector<string>(m_numVectorVars, string("Placeholder name"));
  m_vectorVarNames[m_fluidVel] = string("Fluid velocity");
  m_vectorVarNames[m_fluidVelPred] = string("Fluid velocity prediction");
  m_vectorVarNames[m_gradPressure] = string("Grad pressure");
  m_vectorVarNames[m_fluidVelAnalytic] = string("Fluid vel analytic");
  m_vectorVarNames[m_fluidVelErr] = string("Fluid vel error");
  m_vectorVarNames[m_gradPressureErr] = string("Grad Pressure Err");
}


void
amrMushyLayer::initialize()
{
  CH_TIME("amrMushyLayer::initialize()");


  pout() << "amrMushyLayer::initialize" << endl;


  // first, read in info from parmParse file
  ParmParse pp("main");

  Vector<Real> domsize(SpaceDim);

  // assumption is that domains are not periodic

  for (int dir=0; dir<SpaceDim; dir++)
    m_is_periodic[dir] = false;
  Vector<int> is_periodic_int(SpaceDim, 0);

  int tempBool = m_refine_initial_domain;
  pp.query("refineInitialDomain", tempBool);
  m_refine_initial_domain = (tempBool == 1);

  //  pp.query("solverType", m_solver_type);

  pp.get("max_level", m_max_level);

  m_ncells.resize(SpaceDim);
  pp.getarr("num_cells", m_ncells, 0, SpaceDim);

  if (pp.contains("domain_size"))
  {
    pp.getarr("domain_size", domsize, 0, SpaceDim);
    m_domainSize = RealVect(D_DECL(domsize[0], domsize[1], domsize[2]));
  }


  pp.getarr("periodic_bc", is_periodic_int, 0, SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++)
  {
    m_is_periodic[dir] = (is_periodic_int[dir] == 1);
  }

  pp.query("cfl", m_cfl);

  m_initial_cfl = m_cfl;
  pp.query("initial_cfl", m_initial_cfl);

  pp.query("max_dt_grow_factor", m_max_dt_grow);

  pp.query("max_dt", m_max_dt);
  m_dt = m_max_dt; // set default value

  pp.query("plot_interval", m_plot_interval);

  pp.query("output_folder", m_output_folder);

  pp.query("plot_prefix", m_plot_prefix);

  pp.query("check_interval", m_check_interval);

  pp.query("check_prefix", m_check_prefix);

  pp.query("iteration_plot_interval", m_iteration_plot_interval);

  pp.query("minStep", m_minTimestep);




  if (m_max_level > 0)
  {
    pp.getarr("ref_ratio", m_refinement_ratios, 0, m_max_level);
  }
  else
  {
    m_refinement_ratios.resize(1);
    m_refinement_ratios[0] = -1;
  }

  pp.query("verbosity", s_verbosity);

  pp.get("regrid_interval", m_regrid_interval);

  pp.get("block_factor", m_block_factor);

  pp.get("fill_ratio", m_fill_ratio);

  pp.query("nestingRadius", m_nesting_radius);

  pp.query("tagging_val", m_tagging_val);

  pp.query("tagging_scalar_var", m_tagging_scalar_var);

  pp.query("tagging_vector_var", m_tagging_vector_var);

  pp.query("tags_grow", m_tags_grow);

  pp.query("max_box_size", m_max_box_size);

  pp.query("printAnalyticSoln", m_printAnalyticSoln);

  pp.query("enforceAnalyticSoln", m_enforceAnalyticSoln);

  pp.query("timeIntegrationOrder", m_timeIntegrationOrder);

  pp.query("frameAdvectionMethod", m_frameAdvectionMethod);

  pp.query("fluidVelDiffOrder", m_fluidVelDiffOrder);

  pp.query("fluidVelInterpOrder", m_fluidVelInterpOrder);

  pp.query("ownVelCalc", m_ownVelCalc);

  ParmParse ppDebug("debugging");
  ppDebug.query("fixVelocity", m_debugging_fixVelocity);
  ppDebug.query("resetTemperatureField", m_debugging_resetTemperatureField);

  // now set up problem domains
  {
    IntVect loVect = IntVect::Zero;
    IntVect hiVect(D_DECL(m_ncells[0]-1, m_ncells[1]-1, m_ncells[2]-1));
    ProblemDomain baseDomain(loVect, hiVect);
    // now set periodicity
    for (int dir=0; dir<SpaceDim; dir++)
      baseDomain.setPeriodic(dir, m_is_periodic[dir]);

    // now set up vector of domains
    m_amrDomains.resize(m_max_level+1);
    m_amrDx.resize(m_max_level+1);

    m_amrDomains[0] = baseDomain;
    m_amrDx[0] = m_domainSize[0]/baseDomain.domainBox().longside();

    for (int lev=1; lev<= m_max_level; lev++)
    {
      m_amrDomains[lev] = refine(m_amrDomains[lev-1],
                                 m_refinement_ratios[lev-1]);
      m_amrDx[lev] = m_amrDx[lev-1]/m_refinement_ratios[lev-1];
    }

  } // leaving problem domain setup scope

  // Now we have amr dx we can do this
  m_physBCPtr = new PhysBCUtil(m_parameters, m_amrDx[0]);

  // check to see if we're using predefined grids
  bool usePredefinedGrids = false;
  std::string gridFile;
  if (pp.contains("gridsFile"))
  {
    usePredefinedGrids = true;
    pp.get("gridsFile",gridFile);
  }

  // check to see if we're restarting from a checkpoint file
  if (!pp.contains("restart_file"))
  {
    // if we're not restarting

    // now set up data holders
    //      m_old_phi.resize(m_max_level+1, NULL);
    //      m_new_phi.resize(m_max_level+1, NULL);

    m_scalarOld.resize(m_numVars);
    m_scalarNew.resize(m_numVars);
    m_dScalar.resize(m_numVars);

    m_vectorOld.resize(m_numVectorVars);
    m_vectorNew.resize(m_numVectorVars);
    m_dVector.resize(m_numVectorVars);

    m_frameAdv.resize(m_max_level+1, NULL);
    m_fluidAdv.resize(m_max_level+1, NULL);

    m_fluxRegister.resize(m_numVars);
    m_vectorFluxRegister.resize(m_numVectorVars);
    //    m_velFluxReg.resize(m_max_level+1);

    m_scalarDiffusionCoeffs.resize(m_numVars, 0);

    for (int a_var=0; a_var<m_numVars; a_var++)
    {
      m_scalarOld[a_var].resize(m_max_level+1);
      m_scalarNew[a_var].resize(m_max_level+1);
      m_dScalar[a_var].resize(m_max_level+1);
      m_fluxRegister[a_var].resize(m_max_level+1);

      m_scalarDiffusionCoeffs[m_enthalpy] = 1.0;
      m_scalarDiffusionCoeffs[m_bulkConcentration] = 1/m_parameters.lewis;
    }

    m_HC.resize(m_max_level+1, NULL);

    for (int a_var=0; a_var<m_numVectorVars; a_var++)
    {
      m_vectorOld[a_var].resize(m_max_level+1, NULL);
      m_vectorNew[a_var].resize(m_max_level+1, NULL);
      m_dVector[a_var].resize(m_max_level+1, NULL);
      m_vectorFluxRegister[a_var].resize(m_max_level+1);
    }





    int finest_level = -1;
    if (usePredefinedGrids)
    {
      setupFixedGrids(gridFile);
    }
    else
    {
      // now create initial grids
      initGrids(finest_level);

      // initGrids creates the wrong grid structure for AMRPoissonOp (it was designed for OldAMRPoissonOp).
      // Specifically, m_amrGrids should only be defined for the levels we are using. Any levels finer
      // than those that we're using should just be ignored.
      // I have fixed this in regrid (which is called after each time step)
      // so just going to call it again here rather than rewriting initGrids.

      //todo - might have to turn regrid back on
      //			regrid();
    }

    // that should be it
  }
  else
  {
    // we're restarting from a checkpoint file
    string restart_file;
    pp.get("restart_file", restart_file);
    m_do_restart = true;

    restart(restart_file);

    for (int a_var=0; a_var <m_numVars; a_var++)
    {
      applyBCs(a_var);
      for (int lev=0; lev<=m_finest_level; lev++)
      {
        m_scalarNew[a_var][lev]->exchange();
      }
      applyBCs(a_var);
    }

    setupFluxRegisters();



    Vector<DisjointBoxLayout> activeGrids;
    activeGrids.resize(m_finest_level+1);

    for (int lev=0; lev<=m_finest_level; lev++)
    {
      activeGrids[lev] = m_amrGrids[lev];
    }

    // Need to initialise velocity

    //		updateVelocity(activeGrids);

  } //end if restarting from a checkpoint file




  // set up counter of number of cells
  m_num_cells.resize(m_max_level+1, 0);
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
    LayoutIterator lit = levelGrids.layoutIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& thisBox = levelGrids.get(lit());
      m_num_cells[lev] += thisBox.numPts();
    }
  }


  // finally, set up covered_level flags
  m_covered_level.resize(m_max_level+1, 0);

  // note that finest level can't be covered.
  for (int lev=m_finest_level-1; lev>=0; lev--)
  {

    // if the next finer level is covered, then this one is too.
    if (m_covered_level[lev+1] == 1)
    {
      m_covered_level[lev] = 1;
    }
    else
    {
      // see if the grids finer than this level completely cover it
      IntVectSet fineUncovered(m_amrDomains[lev+1].domainBox());
      const DisjointBoxLayout& fineGrids = m_amrGrids[lev+1];

      LayoutIterator lit = fineGrids.layoutIterator();
      for (lit.begin(); lit.ok(); ++lit)
      {
        const Box& thisBox = fineGrids.get(lit());
        fineUncovered.minus_box(thisBox);
      }

      if (fineUncovered.isEmpty())
      {
        m_covered_level[lev] = 1;
      }
    }
  } // end loop over levels to determine covered levels

  setupIBCS();

  calculateAnalyticSolution();

  if (!m_do_restart)
  {
    if (m_enforceAnalyticSoln)
    {
      enforceAnalyticSolution();

    }

    //May need to regrid now (unless we have a fixed grid)
    if (!usePredefinedGrids)
    {
      //TODO - turn on regridding during initialisation (if not using predefined grids)?
      //regrid();
    }
  }







}


void amrMushyLayer::initVars(const int lev)
{
  for (int a_var = 0; a_var<m_numVars; a_var++)
  {
    m_scalarNew[a_var][lev] = RefCountedPtr<LevelData<FArrayBox> >
    (new LevelData<FArrayBox>(m_amrGrids[lev], 1,m_ghostVect));
    m_scalarOld[a_var][lev] = RefCountedPtr<LevelData<FArrayBox> >
    (new LevelData<FArrayBox>(m_amrGrids[lev], 1,m_ghostVect));
    m_dScalar[a_var][lev] = RefCountedPtr<LevelData<FArrayBox> >
    (new LevelData<FArrayBox>(m_amrGrids[lev], 1,m_ghostVect));
  }

  // HC has two components
  m_HC[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 2 ,m_ghostVect);

  for (int a_var=0; a_var<m_numVectorVars; a_var++)
  {
    m_vectorNew[a_var][lev] = new LevelData<FArrayBox>(m_amrGrids[lev], SpaceDim,m_ghostVect);
    m_vectorOld[a_var][lev] = new LevelData<FArrayBox>(m_amrGrids[lev], SpaceDim,m_ghostVect);
    m_dVector[a_var][lev] = new LevelData<FArrayBox>(m_amrGrids[lev], SpaceDim,m_ghostVect);
  }

  // Changed to 1 component (is this right?)
  m_frameAdv[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, m_ghostVect);
  m_fluidAdv[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, m_ghostVect);
}


void
amrMushyLayer::initGrids(int a_finest_level)
{

  if (s_verbosity > 3)
  {
    pout() << "amrMushyLayer::initGrids" << endl;
  }


  m_finest_level = 0;
  // first create base level
  Vector<Box> baseBoxes;
  domainSplit(m_amrDomains[0], baseBoxes, m_max_box_size,
              m_block_factor);

  Vector<int> procAssign(baseBoxes.size());
  LoadBalance(procAssign,baseBoxes);

  DisjointBoxLayout baseGrids(baseBoxes, procAssign, m_amrDomains[0]);

  m_amrGrids.resize(m_max_level+1);
  m_amrGrids[0] = baseGrids;

  //Initialise variables on level 0
  initVars(0);

  // initialize base level data
  initData();

  //	int numLevels = 1;
  bool moreLevels = (m_max_level > 0);

  int baseLevel = 0;
  //int topLevel = m_finest_level;


  BRMeshRefine meshrefine;
  if (moreLevels)
  {
    meshrefine.define(m_amrDomains[0], m_refinement_ratios,
                      m_fill_ratio, m_block_factor,
                      m_nesting_radius, m_max_box_size);
  }

  Vector<IntVectSet> tagVect(m_max_level);

  Vector<Vector<Box> > oldBoxes(1);
  Vector<Vector<Box> > newBoxes;
  oldBoxes[0] = baseBoxes;
  newBoxes = oldBoxes;
  int new_finest_level = 0;

  while (moreLevels)
  {
    // default is moreLevels = false
    // (only repeat loop in the case where a new level is generated
    // which is still coarser than maxLevel)
    moreLevels = false;
    tagCellsInit(tagVect);

    // two possibilities -- need to generate grids
    // level-by-level, or we are refining all the
    // way up for the initial time.  check to
    // see which it is by seeing if the finest-level
    // tags are empty
    if (tagVect[m_max_level-1].isEmpty())
    {
      int top_level = m_finest_level;
      int old_top_level = top_level;
      new_finest_level = meshrefine.regrid(newBoxes,
                                           tagVect, baseLevel,
                                           top_level,
                                           oldBoxes);

      if (new_finest_level > top_level) top_level++;
      oldBoxes = newBoxes;

      // now see if we need another pass through grid generation
      if ((top_level < m_max_level) && (top_level > old_top_level))
      {
        moreLevels = true;
      }

    }
    else
    {

      // for now, define old_grids as just domains
      oldBoxes.resize(m_max_level+1);
      for (int lev=1; lev<=m_max_level; lev++)
      {
        oldBoxes[lev].push_back(m_amrDomains[lev].domainBox());
      }

      int top_level = m_max_level -1;
      new_finest_level = meshrefine.regrid(newBoxes,
                                           tagVect, baseLevel,
                                           top_level,
                                           oldBoxes);
    }


    broadcast(new_finest_level, uniqueProc(SerialTask::compute));

    //		int numLevels = Min(new_finest_level, m_max_level)+1;

    broadcast(newBoxes, uniqueProc(SerialTask::compute));


    // now loop through levels and define
    for (int lev=baseLevel+1; lev<= new_finest_level; ++lev)
    {
      int numGridsNew = newBoxes[lev].size();
      Vector<int> procIDs(numGridsNew);
      LoadBalance(procIDs, newBoxes[lev]);
      const DisjointBoxLayout newDBL(newBoxes[lev], procIDs,
                                     m_amrDomains[lev]);
      m_amrGrids[lev] = newDBL;

      initVars(lev);
    }

    m_finest_level = new_finest_level;

    // finally, initialize data on final hierarchy
    // only do this if we've created new levels
    if (m_finest_level > 0)
    {
      initData();
    }
  } // end while more levels to do



}



void
amrMushyLayer::initData()
{

  setupFluxRegisters();

  if (s_verbosity > 3)
  {
    pout() << "amrMushyLayer::initData" << endl;
  }

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    initScalarVars(lev);
    initVectorVars(lev);
  }

  // Copy initialised values into old data holders
  for (int a_var=0; a_var<m_numVars; a_var++)
  {
    assert (m_scalarOld[a_var].size() == m_scalarNew[a_var].size());
    activeLevelCopy(m_scalarNew[a_var], m_scalarOld[a_var]);
  }


  Vector<DisjointBoxLayout> activeGrids;
  activeGrids.resize(m_finest_level+1);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    activeGrids[lev] = m_amrGrids[lev];
  }


  // may be necessary to average down here

  //	averageCoarseToFineSolutions();

  // Calculate initial velocity field from the just initialised scalarVars

  setupIBCS();

  updateVelocityComp(activeGrids, 1);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    getFrameAdvection(lev);
  }
}

void amrMushyLayer::
initVectorVars(const int lev)
{
  if (s_verbosity > 3)
    {
      pout() << "amrMushyLayer::initVectorVars" << endl;
    }

  //Initialise vector variables
  for (int a_var=0; a_var<m_numVectorVars; a_var++)
  {
    LevelData<FArrayBox>& levelPhiNew = *(m_vectorNew[a_var][lev]);
    DataIterator levelDit = levelPhiNew.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      FArrayBox& thisPhiNew = levelPhiNew[levelDit()];
      Box thisBox = thisPhiNew.box();

      BoxIterator bit(thisBox);
      for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        //				RealVect loc(iv);
        //				loc *= m_amrDx[lev];
        //				loc += ccOffset;

        //Just set everything to zero for now
        for (int dir=0; dir<SpaceDim; dir++)
        {
          thisPhiNew(iv,dir) = 0;
        }

      }
    }

  }
}

void amrMushyLayer::getLocation(const IntVect iv, const int lev, RealVect &loc, Real ccOffsetX, Real ccOffsetY)
{
  RealVect offset(ccOffsetX, ccOffsetY);
  ::getLocation(iv, loc, m_amrDx[lev], offset);
}

void amrMushyLayer::
initScalarVars(const int lev)
{
  if (s_verbosity > 3)
    {
      pout() << "amrMushyLayer::initScalarVars" << endl;
    }

  vector<int> ignore;

  //Initialise scalar variables

    LevelData<FArrayBox>& levelHC = *(m_HC[lev]);
    DataIterator levelDit = levelHC.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      levelHC[levelDit()].setVal(m_parameters.Hinitial, 0);
      levelHC[levelDit()].setVal(m_parameters.bcValBulkConcentrationLo[SpaceDim-1], 0);

//      Box thisBox = thisPhiNew.box();
//      IntVect box_iv = thisBox.smallEnd();
//
//      BoxIterator bit(thisBox);
//      for (bit.begin(); bit.ok(); ++bit)
//      {
//        IntVect iv = bit();
//        RealVect loc;
//        getLocation(iv, lev, loc);
//        Real x = loc[0];
//        Real z = loc[1];
//
//        if (a_var == m_HC)
//        {
//
//          thisPhiNew(iv,0) = ;
//          thisPhiNew(iv,1) = ;
//
//        }
//        else
//        {
//          thisPhiNew(iv,0) = 0;
//        }
//      } //end loop over intvects


    } //end loop over boxes

  updateEnthalpyVariables();


  //Copy new to old vars
  //	for (int a_var=0; a_var<m_numVars; a_var++)
  //	{
  //		activeLevelCopy(m_scalarNew[a_var], m_scalarOld[a_var]);
  //	}
}


void amrMushyLayer::setupFluxRegisters()
{
  //Setup flux registers
  // Flux register lives on the coarser level, so don't have one for the max level
  // as there isn't a max_level+1
  for (int lev=0; lev<m_finest_level; lev++)
  {
    for (int a_var=0; a_var<m_numVars; a_var++)
    {
      m_fluxRegister[a_var][lev] = RefCountedPtr<LevelFluxRegister>
      (new LevelFluxRegister(m_amrGrids[lev+1], m_amrGrids[lev], m_amrDomains[lev+1], m_refinement_ratios[lev], 1));
      m_fluxRegister[a_var][lev]->setToZero();
    }

    for (int a_var=0; a_var<m_numVectorVars; a_var++)
    {
      m_vectorFluxRegister[a_var][lev] = RefCountedPtr<LevelFluxRegister>
      (new LevelFluxRegister(m_amrGrids[lev+1], m_amrGrids[lev], m_amrDomains[lev+1], m_refinement_ratios[lev], SpaceDim));
      m_vectorFluxRegister[a_var][lev]->setToZero();
    }

//    m_velFluxReg[lev] = RefCountedPtr<LevelFluxRegisterEdge>
//    (new LevelFluxRegisterEdge(m_amrGrids[lev+1], m_amrGrids[lev], m_amrDomains[lev+1], m_refinement_ratios[lev], 1));
//    m_velFluxReg[lev]->setToZero();
  }
}

void
amrMushyLayer::setupFixedGrids(const std::string& a_gridFile)
{
  Vector<Vector<Box> > gridvect;

  if (procID() == uniqueProc(SerialTask::compute))
  {
    gridvect.push_back(Vector<Box>(1,m_amrDomains[0].domainBox()));

    // read in predefined grids
    ifstream is(a_gridFile.c_str(), ios::in);

    if (is.fail())
    {
      MayDay::Error("Cannot open grids file");
    }

    // format of file:
    //   number of levels, then for each level (starting with level 1):
    //   number of grids on level, list of boxes
    int inNumLevels;
    is >> inNumLevels;

    if (s_verbosity > 3)
    {
      pout() << "numLevels = " << inNumLevels << endl;
      pout() << "max_level = " << m_max_level << endl;
    }

    assert (inNumLevels <= m_max_level+1);


    while (is.get() != '\n');

    gridvect.resize(inNumLevels);

    // check to see if coarsest level needs to be broken up
    domainSplit(m_amrDomains[0],gridvect[0], m_max_box_size,
                m_block_factor);

    if (s_verbosity >= 3)
    {
      pout() << "level 0: ";
      for (int n=0; n < gridvect[0].size(); n++)
      {
        pout() << gridvect[0][n] << endl;
      }
    }



    // now loop over levels, starting with level 1
    int numGrids = 0;
    for (int lev=1; lev<inNumLevels; lev++)
    {
      is >> numGrids;

      if (s_verbosity >= 3)
      {
        pout() << "level " << lev << " numGrids = "
            << numGrids <<  endl;
        pout() << "Grids: ";
      }

      while (is.get() != '\n');

      gridvect[lev].resize(numGrids);

      for (int i=0; i<numGrids; i++)
      {
        Box bx;
        is >> bx;

        while (is.get() != '\n');

        // quick check on box size
        Box bxRef(bx);

        if (bxRef.longside() > m_max_box_size)
        {
          pout() << "Grid " << bx << " too large" << endl;
          MayDay::Error();
        }

        if (s_verbosity >= 3)
        {
          pout() << bx << endl;
        }

        gridvect[lev][i] = bx;
      } // end loop over boxes on this level



    } // end loop over levels
  } // end if serial proc

  // broadcast results
  broadcast(gridvect, uniqueProc(SerialTask::compute));

  // now create disjointBoxLayouts and allocate grids

  m_amrGrids.resize(m_max_level+1);
  for (int lev=0; lev<gridvect.size(); lev++)
  {
    int numGridsLev = gridvect[lev].size();
    Vector<int> procIDs(numGridsLev);
    LoadBalance(procIDs, gridvect[lev]);
    const DisjointBoxLayout newDBL(gridvect[lev],
                                   procIDs,
                                   m_amrDomains[lev]);

    m_amrGrids[lev] = newDBL;

    // Initialise global variables on base level
    initVars(lev);

  }

  // finally set finest level and initialize data on hierarchy
  m_finest_level = gridvect.size() -1;

  // NB this must come before initData()
  setupFluxRegisters();
  setupIBCS();

  initData();



}






void
amrMushyLayer::setupAdvectionSolvers()
{
  CH_TIME("amrMushyLayer::setupAdvectionSolvers()");

  m_godunovtheta.resize(m_finest_level+1);
  m_godunovConcLiquid.resize(m_finest_level+1);
  m_godunovConc.resize(m_finest_level+1);
  m_godunovEnthalpy.resize(m_finest_level+1);

  for (int lev = 0; lev <= m_finest_level; lev++)
  {
    bool hasCoarser = (lev > 0);
    bool hasFiner = (lev < m_finest_level);

    bool use_limiting = true; // whether to do slope limiting in the primitive variable

    //This isn't very clear code, but can't declare a reference value without initializing so...
    const DisjointBoxLayout& coarseGrid = (hasCoarser) ? m_amrGrids[lev-1] :
        DisjointBoxLayout();
    int refRatio = (hasCoarser) ? m_refinement_ratios[lev-1]:
        m_refinement_ratios[lev];

    //LevelAdvect won't accept refRatio < 0, but this code sets refRatio = -1 if m_max_level = 0, so this is a hack
    //Don't want to change original code for fear of breaking stuff.
    if (m_max_level == 0)
    {
      refRatio = 2;
    }


    m_godunovtheta[lev] = RefCountedPtr<LevelAdvect>(new LevelAdvect());
    m_godunovtheta[lev]->define(*m_advPhystheta[lev],
                                m_amrGrids[lev],
                                coarseGrid,
                                m_amrDomains[lev],
                                refRatio,
                                use_limiting,
                                m_amrDx[lev],
                                hasCoarser,
                                hasFiner,
                                m_num_ghost);


    m_godunovConcLiquid[lev]= RefCountedPtr<LevelAdvect>(new LevelAdvect());
    m_godunovConcLiquid[lev]->define(*m_advPhysConcLiquid[lev],
                                     m_amrGrids[lev],
                                     coarseGrid,
                                     m_amrDomains[lev],
                                     refRatio,
                                     use_limiting,
                                     m_amrDx[lev],
                                     hasCoarser,
                                     hasFiner,
                                     m_num_ghost);

    m_godunovConc[lev]= RefCountedPtr<LevelAdvect>(new LevelAdvect());
    m_godunovConc[lev]->define(*m_advPhysConc[lev],
                               m_amrGrids[lev],
                               coarseGrid,
                               m_amrDomains[lev],
                               refRatio,
                               use_limiting,
                               m_amrDx[lev],
                               hasCoarser,
                               hasFiner,
                               m_num_ghost);

    m_godunovEnthalpy[lev]= RefCountedPtr<LevelAdvect>(new LevelAdvect());
    m_godunovEnthalpy[lev]->define(*m_advPhysEnthalpy[lev],
                                   m_amrGrids[lev],
                                   coarseGrid,
                                   m_amrDomains[lev],
                                   refRatio,
                                   use_limiting,
                                   m_amrDx[lev],
                                   hasCoarser,
                                   hasFiner,
                                   m_num_ghost);





  }
}

void
amrMushyLayer::setupIBCS()
{
  //todo - fix this to do advection
  if (s_verbosity > 3)
  {
    pout() << "amrMushyLayerInit::setupIBCs" << endl;
  }


  m_advPhystheta.resize(m_max_level+1);
  m_advPhysConcLiquid.resize(m_max_level+1);
  m_advPhysConc.resize(m_max_level+1);
  m_advPhysEnthalpy.resize(m_max_level+1);
  m_advPhysLambda.resize(m_max_level+1);

  for (int lev=0; lev<= m_max_level; lev++)
  {
    /*
    //Create IBC's on this level
    RefCountedPtr<AdvectConcIBC> ibcConcLiquid = RefCountedPtr<AdvectConcIBC>
    (new AdvectConcIBC(m_parameters.ThetaLInitial, m_parameters.ThetaLBottom, m_parameters.ThetaLTop));
    RefCountedPtr<AdvectHeatIBC> ibctheta = RefCountedPtr<AdvectHeatIBC>
    (new AdvectHeatIBC(m_parameters.thetaInitial, m_parameters.thetaBottom, m_parameters.thetaTop));
    RefCountedPtr<AdvectHeatIBC> ibcEnthalpy = RefCountedPtr<AdvectHeatIBC>
    (new AdvectHeatIBC(m_parameters.Hinitial, m_parameters.HBottom, m_parameters.HTop));
    RefCountedPtr<AdvectHeatIBC> ibcConc = RefCountedPtr<AdvectHeatIBC>
    (new AdvectHeatIBC(m_parameters.ThetaInitial, m_parameters.ThetaTop, m_parameters.ThetaBottom));
    RefCountedPtr<AdvectHeatIBC> ibcLambda = RefCountedPtr<AdvectHeatIBC>
    (new AdvectHeatIBC(1,1,1));

    ibcConcLiquid->define(m_amrDomains[lev], m_amrDx[lev]);
    ibctheta->define(m_amrDomains[lev], m_amrDx[lev]);
    ibcEnthalpy->define(m_amrDomains[lev], m_amrDx[lev]);
    ibcConc->define(m_amrDomains[lev], m_amrDx[lev]);
    ibcLambda->define(m_amrDomains[lev], m_amrDx[lev]);

    //Create advection physics and apply IBC's
    AdvectPhysics advPhystheta, advPhysConcLiquid, advPhysEnthalpy, advPhysConc, advPhysLambda;
    advPhystheta.define(m_amrDomains[lev], m_amrDx[lev]);
    advPhysConcLiquid.define(m_amrDomains[lev], m_amrDx[lev]);
    advPhysConc.define(m_amrDomains[lev], m_amrDx[lev]);
    advPhysEnthalpy.define(m_amrDomains[lev], m_amrDx[lev]);
    advPhysLambda.define(m_amrDomains[lev], m_amrDx[lev]);

    advPhystheta.setPhysIBC(&(*ibctheta));
    advPhysConcLiquid.setPhysIBC(&(*ibcConcLiquid));
    advPhysEnthalpy.setPhysIBC(&(*ibcEnthalpy));
    advPhysConc.setPhysIBC(&(*ibcConc));
    advPhysLambda.setPhysIBC(&(*ibcLambda));

    // Derive godunov physics
    GodunovPhysics* gphysPtrtheta = advPhystheta.new_godunovPhysics();
    GodunovPhysics* gphysPtrConcLiquid = advPhysConcLiquid.new_godunovPhysics();
    GodunovPhysics* gphysPtrEnthalpy = advPhysEnthalpy.new_godunovPhysics();
    GodunovPhysics* gphysPtrConc = advPhysConc.new_godunovPhysics();
    GodunovPhysics* gphysPtrLambda = advPhysLambda.new_godunovPhysics();

    // Save objects to vector
    m_advPhystheta[lev] = RefCountedPtr<AdvectPhysics>((AdvectPhysics*)gphysPtrtheta);
    m_advPhysConcLiquid[lev] = RefCountedPtr<AdvectPhysics>((AdvectPhysics*)gphysPtrConcLiquid);
    m_advPhysConc[lev] = RefCountedPtr<AdvectPhysics>((AdvectPhysics*)gphysPtrConc);
    m_advPhysEnthalpy[lev] = RefCountedPtr<AdvectPhysics>((AdvectPhysics*)gphysPtrEnthalpy);
    m_advPhysLambda[lev] = RefCountedPtr<AdvectPhysics>((AdvectPhysics*)gphysPtrLambda);
*/
  }

}

void
amrMushyLayer::setupPoissonSolvers(
    Vector<LevelData<FArrayBox>* >& a_thetaNew,
    Vector<LevelData<FArrayBox>* >& thetaSource,
    Vector<LevelData<FArrayBox>* >& a_ThetaLNew,
    Vector<LevelData<FArrayBox>* >& ThetaLSource,
    Vector<DisjointBoxLayout>& activeGrids)
{
  CH_TIME("amrMushyLayer::setupPoissonSolvers()");

  //  s_bottomSolver.m_verbosity = 0;
  //
  //  //Poisson solver for theta
  //  Real alpha = 1.0;
  //  Real beta  = 1.0;
  //
  //  //This object is used for the BE solver
  //  AMRPoissonOpFactory* opFact = new AMRPoissonOpFactory();
  //  opFact->define(m_amrDomains[0],
  //                 activeGrids,
  //                 m_refinement_ratios,
  //                 m_amrDx[0],
  //                 m_physBCPtr->thetaFuncBC(), alpha, beta);
  //
  //  //This object is used for calculating residuals
  //  thetaPoissonOpFact = RefCountedPtr<AMRPoissonOpFactory>(opFact);
  //  thetaPoissonOpFact->define(m_amrDomains[0],
  //                             activeGrids,
  //                             m_refinement_ratios,
  //                             m_amrDx[0],
  //                             m_physBCPtr->thetaFuncBC(), alpha, beta);
  //
  //
  //  // Variable coefficient solver for theta
  //  Real thetaVCAlpha, thetaVCBeta;
  //  Vector<RefCountedPtr<LevelData<FArrayBox> > > thetaVCa;
  //  Vector<RefCountedPtr<LevelData<FluxBox> > >  thetaVCb;
  //  thetaVCOpCoeffs(thetaVCAlpha, thetaVCBeta, thetaVCa, thetaVCb);
  //  thetaVCOpFact = RefCountedPtr<VCAMRPoissonOp2Factory>(new VCAMRPoissonOp2Factory());
  //  thetaVCOpFact->define(m_amrDomains[0],
  //                        activeGrids,
  //                        m_refinement_ratios,
  //                        m_amrDx[0],
  //                        m_physBCPtr->thetaFuncBC(),
  //                        thetaVCAlpha, thetaVCa,
  //                        thetaVCBeta, thetaVCb);
  //
  //
  //  // Variable Coefficient solver for Theta_l
  //  ThetaLVCOpFact = RefCountedPtr<VCAMRPoissonOp2Factory>(new VCAMRPoissonOp2Factory());
  //
  //  Real ThetaLAlpha, ThetaLBeta;
  //  Vector<RefCountedPtr<LevelData<FArrayBox> > > ThetaLa;
  //  Vector<RefCountedPtr<LevelData<FluxBox> > > ThetaLb;
  //  getThetaLOpCoeffs(ThetaLAlpha, ThetaLBeta, ThetaLa, ThetaLb);
  //
  //  ThetaLVCOpFact->define(m_amrDomains[0],
  //                         activeGrids,
  //                         m_refinement_ratios,
  //                         m_amrDx[0],
  //                         //ADParseBCConcLiquid,
  //                         m_physBCPtr->ThetaLFuncBC(),
  //                         ThetaLAlpha, ThetaLa,
  //                         ThetaLBeta, ThetaLb);
  //
  //
  //
  //
  //





}

void
amrMushyLayer::setupMultigrid(Vector<DisjointBoxLayout>& activeGrids)
{

  CH_TIME("amrMushyLayer::setupMultigrid()");

  //Get general parameters for multigrid solve
  int lbase = 0;
//  int lmax = m_finest_level;

  int numSmooth, numMG, maxIter, mgverb, numLevels;
  Real tolerance, hang, normThresh;

  numLevels = m_finest_level+1;

  ParmParse pp("amrmultigrid");
  pp.get("num_smooth", numSmooth);
  pp.get("num_mg",     numMG);
  pp.get("hang_eps",   hang);
  pp.get("norm_thresh",normThresh);
  pp.get("tolerance",  tolerance);
  pp.get("max_iter",   maxIter);
  pp.get("verbosity",  mgverb);

//  LinearSolver<LevelData<FArrayBox> >* bottomSolverPtrVC = &s_bottomSolverVC;
  LinearSolver<LevelData<FArrayBox> >* bottomSolver = &s_bottomSolver;

  //Multigrid solver for calculating reflux correction
  m_refluxOpFact = RefCountedPtr<AMRPoissonOpFactory>(new AMRPoissonOpFactory());
  m_refluxOpFact->define(m_amrDomains[0],
                         activeGrids,
                         m_refinement_ratios,
                         m_amrDx[0],
                         m_physBCPtr->noFluxBC());

  m_refluxAMRMG = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > > (new AMRMultiGrid<LevelData<FArrayBox> >);
  m_refluxAMRMG->define(m_amrDomains[0],
                        *m_refluxOpFact,
                        bottomSolver,
                        numLevels);
  m_refluxAMRMG->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG,
                                     maxIter, tolerance, hang, normThresh);
  m_refluxAMRMG->m_verbosity=mgverb;



  //	Multigrid solver for theta variable coefficient operator
  //  m_amrSolverVCtheta = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > > (new AMRMultiGrid<LevelData<FArrayBox> >);
  //  //	m_amrSolverVCtheta->m_maxDepth = 1; // Trying to deal with periodic BCs and > 2 ghost cells
  //  m_amrSolverVCtheta->define(m_amrDomains[0],
  //                             *thetaVCOpFact,
  //                             bottomSolverPtrVC,
  //                             numLevels);
  //  m_amrSolverVCtheta->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG,
  //                                          maxIter, tolerance, hang, normThresh);
  //  m_amrSolverVCtheta->m_verbosity=mgverb;

  // TODO - fill these
  Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(numLevels);
  Vector<RefCountedPtr<LevelData<FluxBox> > > porosityFace(numLevels);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(numLevels);


  EdgeVelBCHolder porosityEdgeBC(m_physBCPtr->porosityFaceBC());

//  updateEnthalpyVariables();

  // two components: enthalpy and salinity
  int numComps = 2;
  int Hcomp = 0;
  int Ccomp = 1;

  //todo - check if we need this ghost vector
  IntVect ivGhost = IntVect::Unit;


  for (int lev=0; lev<numLevels; lev++)
  {

    bCoef[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(activeGrids[lev], numComps, ivGhost));
    aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(activeGrids[lev], numComps, ivGhost));
    porosityFace[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(activeGrids[lev], 1, ivGhost));

    CellToEdge(*m_scalarNew[m_porosity][lev], *porosityFace[lev]);

    for (DataIterator dit = bCoef[lev]->dataIterator(); dit.ok(); ++dit)
    {
      (*aCoef[lev])[dit].setVal(1.0);
      (*bCoef[lev])[dit].setVal(1.0);

      // bCoef for salt solve
//      (*bCoef[lev])[dit].mult((*porosityFace[lev])[dit], (*porosityFace[lev])[dit].box(), 0, Ccomp);

      // for heat solve
//      (*bCoef[lev])[dit].minus((*porosityFace[lev])[dit], (*porosityFace[lev])[dit].box(), 0, Hcomp);
//      for (int dir=0; dir<SpaceDim; dir++)
//      {
//        (*bCoef[lev])[dit][dir].mult(m_parameters.heatConductivityRatio, Hcomp);
//      }
//      (*bCoef[lev])[dit].plus((*porosityFace[lev])[dit], (*porosityFace[lev])[dit].box(), 0, Hcomp);
//

      for (int dir=0; dir<SpaceDim; dir++)
      {
        (*bCoef[lev])[dit][dir].mult(-m_scalarDiffusionCoeffs[m_enthalpy], Hcomp);
        (*bCoef[lev])[dit][dir].mult(-m_scalarDiffusionCoeffs[m_bulkConcentration], Ccomp);
      }
    }

  } // end loop over levels

  MushyLayerParams* mlParamsPtr = &m_parameters;

  //        EnthalpyVariable calcTemperature = computeTemperatureFunc;
  BCHolder temperature_Sl_BC; // = m_physBCPtr->BasicthetaFuncBC();
  temperature_Sl_BC = m_physBCPtr->temperatureLiquidSalinityBC();
  BCHolder HC_BC  = m_physBCPtr->enthalpySalinityBC();

  // todo - check these values are correct
  Real alpha = 1.0; Real beta = 1.0;
  int relaxMode = 1;

  AMRNonLinearMultiCompOpFactory* HCop = new AMRNonLinearMultiCompOpFactory();
  HCop->define(m_amrDomains[0], activeGrids, m_refinement_ratios, m_amrDx[0], HC_BC,
               alpha, aCoef, beta, bCoef,
               m_scalarNew[m_enthalpySolidus], m_scalarNew[m_enthalpyLiquidus],
               m_scalarNew[m_enthalpyEutectic],
               mlParamsPtr, temperature_Sl_BC,
               relaxMode, porosityEdgeBC);

  RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >HCOpFact;
  HCOpFact = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(
      HCop);

//  int maxAMRlevels = numLevels;

  s_multiCompFASMG = RefCountedPtr<AMRFASMultiGrid<LevelData<FArrayBox> > >(
      new AMRFASMultiGrid<LevelData<FArrayBox> >());

//  s_multiCompFASMG->define(m_amrDomains[0], *HCOpFact, &bottomSolve,
//                           maxAMRlevels);
//  s_multiCompFASMG->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG,
//                                        maxIter, tolerance, hang, normThresh);
//  s_multiCompFASMG->m_verbosity = mgverb;



// Doesn't currently work. Potential issue is that LevelTGA dosn't like my custom operator.
//  s_enthalpySalinityTGA = RefCountedPtr<LevelTGA>(
//      new LevelTGA(activeGrids, m_refinement_ratios, m_amrDomains[0], HCOpFact,
//                   s_multiCompFASMG));


  s_enthalpySalinityBE = RefCountedPtr<LevelBackwardEuler>(
      new LevelBackwardEuler(activeGrids, m_refinement_ratios, m_amrDomains[0], HCOpFact,
                             s_multiCompFASMG));


  //Need to do this for TGA solver
  if (m_timeIntegrationOrder == 2)
  {
    //todo- init enthalpy salinity solve
    //		amrSolver->init(a_thetaNew,
    //				thetaSource,
    //				lmax,
    //				lbase);

    //    m_amrSolverVCtheta->init(a_thetaNew,
    //                             thetaSource,
    //                             lmax,
    //                             lbase);
    //
    //    m_amrSolverVCThetaL->init(a_ThetaLNew,
    //                              ThetaLSource,
    //                              lmax,
    //                              lbase);
  }

  //Backward Euler solver for timestepping theta
  m_BEEnthalpySalinity = RefCountedPtr<BackwardEuler>
  (new BackwardEuler(
      s_enthalpySalinityTGA,
      *HCop,
      m_amrDomains[0],
      m_refinement_ratios,
      numLevels,
      s_verbosity-1));


  //TGA solver for timestepping theta
  m_TGAEnthalpySalinity =  RefCountedPtr<AMRTGA<LevelData<FArrayBox> > >
  (new AMRTGA<LevelData<FArrayBox> >(
      s_enthalpySalinityTGA,
      *HCop,
      m_amrDomains[0],
      m_refinement_ratios,
      numLevels,
      s_verbosity-1));




}


void amrMushyLayer::
getThetaLOpCoeffs(Real& alpha, Real& beta,
                  Vector<RefCountedPtr<LevelData<FArrayBox> > >& aCoef,
                  Vector<RefCountedPtr<LevelData<FluxBox> > >& bCoef)
{
  alpha = 1;
  beta = 1;

  Vector<LevelData<FArrayBox>* > b;

  aCoef.resize(m_finest_level+1);
  bCoef.resize(m_finest_level+1);
  b.resize(m_finest_level+1, NULL);

  Vector<LevelData<FArrayBox>* > averagePorosity = timeCenteredScalar(m_porosity);

  for (int lev=0; lev<=m_finest_level; lev++)
  {

    aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >
    (new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_ghostVect));

    bCoef[lev] = RefCountedPtr<LevelData<FluxBox> >
    (new LevelData<FluxBox>(m_amrGrids[lev], 1, m_ghostVect));

    //Copy porosity to the B coefficient. Do this explicitly to ensure ghost cells are copied.
    for (DataIterator dit=m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    {

      FArrayBox& fab = (*m_scalarNew[m_ThetaBCoef][lev])[dit];
      FArrayBox& fab2 = (*m_scalarNew[m_porosity][lev])[dit];
      fab.copy(fab2);
    }

    for (DataIterator dit=m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    {
      (*aCoef[lev])[dit].setVal(0);
      (*aCoef[lev])[dit] += (*averagePorosity[lev])[dit];

      (*m_scalarNew[m_ThetaBCoef][lev])[dit].mult(-1/m_parameters.lewis);

      FluxBox& fbox = (*bCoef[lev])[dit];
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
      {
        FArrayBox& fluxDir = fbox[faceDir];
        fluxDir.setVal(0);
      }

      //			Box bBox = (*bCoef[lev])[dit].box();
      //			int tempVal=0;

    }

    // This deals with boundaries between boxes on a level
    aCoef[lev]->exchange();
    bCoef[lev]->exchange();
    m_scalarNew[m_ThetaBCoef][lev]->exchange();


    //Turn FArrayBox into FluxBox for B coefficient
    //LevelData<FluxBox> *fluxBox = &(*bCoef[lev]);
    //		for (DataIterator dit=m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    //		{
    //			//					FluxBox& fbox = (*fluxBox)[dit];
    //			FluxBox& fbox = (*bCoef[lev])[dit];
    //
    //			FArrayBox& fab = (*m_scalarNew[m_ThetaBCoef][lev])[dit];
    //			FArrayBox& fab2 = (*m_scalarNew[m_porosity][lev])[dit];
    //
    //			int temp=1;
    //		}


    CellToEdge(*m_scalarNew[m_ThetaBCoef][lev], *bCoef[lev]);

    //		for (DataIterator dit=m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    //		{
    //			//			FluxBox& fbox = (*fluxBox)[dit];
    //			FluxBox& fbox = (*bCoef[lev])[dit];
    //			int temp=1;
    //		}


  }

  //Clean up memory
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    if (averagePorosity[lev] != NULL)
    {
      delete averagePorosity[lev];
      averagePorosity[lev] = NULL;
    }
  }
}

/*
 * initialise ops for solving (alpha * aCoef(x) * I - beta * Div(bCoef(x).Grad))phi = rho
 * note that later, alpha and beta are reset by the BE/TGA solvers.
 * We set
 * 	aCoef = (c_p - 1) d(porosity)/dt,
 * 	bCoef = -(porosity + k*(1-porosity))
 * 	alpha = beta = 1
 *
 * Such that we are solving
 * 	phi + Div( (porosity + k(1-porosity)) Grad) phi = rhs
 */
void amrMushyLayer::
thetaVCOpCoeffs(Real& alpha, Real& beta,
                Vector<RefCountedPtr<LevelData<FArrayBox> > >& aCoef,
                Vector<RefCountedPtr<LevelData<FluxBox> > >& bCoef)
{

  //I think these are reset later anyway, so doesn't really matter what we put here
  alpha=1;
  beta=1;

  aCoef.resize(m_finest_level+1);
  bCoef.resize(m_finest_level+1);

  Vector<LevelData<FArrayBox>* > chi = timeCenteredScalar(m_porosity);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >
    (new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_num_ghost*IntVect::Unit ));

    bCoef[lev] = RefCountedPtr<LevelData<FluxBox> >
    (new LevelData<FluxBox>(m_amrGrids[lev], 1, m_num_ghost*IntVect::Unit ));

    //		bCoef[lev].neverDelete();

    LevelData<FluxBox> *fluxBox = &(*bCoef[lev]);

    //Create FArray box of a and b coefficients
    for (DataIterator dit=m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    {
      Real k = m_parameters.heatConductivityRatio;
      (*m_scalarNew[m_thetaBcoef][lev])[dit].setVal(1);
      (*m_scalarNew[m_thetaBcoef][lev])[dit] -= (*m_scalarNew[m_porosity][lev])[dit];
      (*m_scalarNew[m_thetaBcoef][lev])[dit].mult(k);
      (*m_scalarNew[m_thetaBcoef][lev])[dit] += (*m_scalarNew[m_porosity][lev])[dit];

      //Ensure we have the right minus sign!
      (*m_scalarNew[m_thetaBcoef][lev])[dit].mult(-1);

      //aCoef = chi + c_p(1-chi)
      (*aCoef[lev])[dit].setVal(1);
      (*aCoef[lev])[dit] -= (*chi[lev])[dit];
      (*aCoef[lev])[dit].mult(m_parameters.specificHeatRatio);
      (*aCoef[lev])[dit] += (*chi[lev])[dit];

      FluxBox& fbox = (*fluxBox)[dit];
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
      {
        FArrayBox& fluxDir = fbox[faceDir];
        fluxDir.setVal(0);
      }
    }

    aCoef[lev]->exchange();
    bCoef[lev]->exchange();

    //Turn FArrayBox into FluxBox for B coefficient
    //LevelData<FluxBox> *fluxBox = &(*bCoef[lev]);
    CellToEdge(*m_scalarNew[m_thetaBcoef][lev], *fluxBox);


  }

  //Clean up memory
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    if (chi[lev] != NULL)
    {
      delete chi[lev];
      chi[lev] = NULL;
    }
  }
}



void
amrMushyLayer::getLevelAdvectionsVars(const int a_var,
                                      const int lev,
                                      int& refRatio,
                                      const bool hasCoarser,
                                      const bool hasFiner,
                                      const bool use_limiting,
                                      RefCountedPtr<LevelData<FArrayBox> >& coarserDataOldPtr,
                                      RefCountedPtr<LevelData<FArrayBox> >& coarserDataNewPtr,
                                      RefCountedPtr<LevelFluxRegister>& coarserFRPtr,
                                      RefCountedPtr<LevelFluxRegister>& finerFRPtr,
                                      DisjointBoxLayout& coarseGrid,
                                      LevelAdvect& levAdvect,
                                      Real& tCoarserNew,
                                      Real& tCoarserOld)
{

  //	refRatio = (lev == 0) ? m_refinement_ratios[0] : m_refinement_ratios[lev-1];
  //
  //	if (refRatio < 0)
  //	{
  //		refRatio = 2;
  //	}
  refRatio = getRefinementRatio(lev-1);

  tCoarserNew = m_time+m_dt;
  tCoarserOld = m_time;

  //Ensure these pointers aren't null
  coarserFRPtr = RefCountedPtr<LevelFluxRegister>(new LevelFluxRegister());
  finerFRPtr = RefCountedPtr<LevelFluxRegister>(new LevelFluxRegister());

  coarserDataOldPtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>());
  coarserDataNewPtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>());

  if (hasCoarser)
  {
    coarserDataOldPtr = m_scalarOld[a_var][lev-1];
    coarserDataNewPtr = m_scalarNew[a_var][lev-1];
    coarseGrid = m_amrGrids[lev-1];
    coarserFRPtr = m_fluxRegister[a_var][lev-1];
  }

  if (hasFiner)
  {
    finerFRPtr = m_fluxRegister[a_var][lev];
  }

  RefCountedPtr<AdvectPhysics> advPhys;
  if (a_var == m_theta)
  {
    advPhys = m_advPhystheta[lev];
  }
  else if(a_var == m_compositionLiquid)
  {
    advPhys = m_advPhysConcLiquid[lev];
  }
  else if(a_var == m_bulkConcentration)
  {
    advPhys = m_advPhysConc[lev];
  }
  else if(a_var == m_enthalpy)
  {
    advPhys = m_advPhysEnthalpy[lev];
  }
  else if (a_var == m_lambda or a_var==m_lambdaPostCorr)
  {
    advPhys = m_advPhysLambda[lev];
  }
  else
  {
    MayDay::Error("amrMushyLayerInit::getLevelAdvectionVars() - unknown variable");
  }

  levAdvect.define(*advPhys,
                   m_amrGrids[lev],
                   coarseGrid,
                   m_amrDomains[lev],
                   refRatio,
                   use_limiting,
                   m_amrDx[lev],
                   hasCoarser,
                   hasFiner,
                   m_num_ghost);
}

// Load data for a field from a file
void amrMushyLayer::loadFromFile(const int a_var, const string filename, const int num_ghost, const int comp)
{

//  vector <vector <double> > data;
//  char *fname = const_cast<char*>(filename.c_str());
//  ifstream infile( fname );
//
//  while (infile)
//  {
//    string s;
//    if (!getline( infile, s )) break;
//
//    istringstream ss( s );
//    vector <double> record;
//
//
//    while (ss)
//    {
//      string s;
//      if (!getline( ss, s, ',' )) break;
//
//      record.push_back( atof(s.c_str()) );
//
//    }
//
//    data.push_back( record );
//
//  }
//
//  //If our data doesn't have the same number of ghost cells as our FAB, we need to introduce an offset
//  int ghostx = (*m_scalarNew[a_var][0]).ghostVect()[0];
//  int ghosty = (*m_scalarNew[a_var][0]).ghostVect()[1];
//
//  int offsetx = ghostx - num_ghost;
//  int offsety = ghosty - num_ghost;
//
//  int maxRow = data[0].size() - 1;
//  int maxCol = data.size() - 1;
//
//  DataIterator dit = (*m_scalarNew[a_var][0]).dataIterator();
//  for (dit.reset(); dit.ok(); ++dit)
//  {
//    FArrayBox& fab = (*m_scalarNew[a_var][0])[dit()];
//
//    // Only fill interior cells
//    BoxIterator bit = BoxIterator(m_amrGrids[0][dit]);
//
//    for(bit.reset(); bit.ok(); ++bit)
//    {
//      IntVect iv = bit();
//
//      int row = iv[0] + num_ghost;
//      int col = iv[1] + num_ghost;
//
//      if (row < 0 || col < 0 || row > maxRow || col > maxCol)
//      {
//        continue;
//      }
//
//      fab(iv, comp) = data[col][row];
//
//    }
//  }
//
//  if (!infile.eof())
//  {
//    cerr << "Error reading in file!\n";
//  }
//
//  // Fill finer levels by refinement if necessary
//  for (int lev=1; lev<=m_finest_level; lev++)
//  {
//    FineInterp interpolatorScalar(m_amrGrids[lev], 1,
//                                  m_refinement_ratios[lev-1],
//                                  m_amrDomains[lev]);
//
//    interpolatorScalar.interpToFine(*m_scalarNew[a_var][lev], *m_scalarNew[a_var][lev-1]);
//  }
//
//  // Fill ghost cells
//  for (int lev=0; lev <= m_finest_level; lev++)
//  {
//    for (DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
//    {
//      FArrayBox& phi = (*m_scalarNew[a_var][lev])[dit];
//      BCHolder BC;
//
//      if (a_var == m_porosity)
//      {
//        BC = m_physBCPtr->BasicPorosityFuncBC();
//        BC(phi, phi.box(), m_amrDomains[lev], m_amrDx[lev], false);
//      }
//      if (a_var == m_theta)
//      {
//        BC = m_physBCPtr->BasicthetaFuncBC();
//        BC(phi, phi.box(), m_amrDomains[lev], m_amrDx[lev], false);
//      }
//      if (a_var == m_HC)
//      {
//        BC = m_physBCPtr->BasicEnthalpyFuncBC();
//        BC(phi, phi.box(), m_amrDomains[lev], m_amrDx[lev], false);
//      }
//      if (a_var == m_bulkConcentration)
//      {
//        BC = m_physBCPtr->BasicThetaFuncBC();
//        BC(phi, phi.box(), m_amrDomains[lev], m_amrDx[lev], false);
//      }
//
//
//
//
//
//    }
//
//    Box b = m_amrDomains[lev].domainBox();
//
//    // Grow box to include cells we wish to fill
//    b.grow(m_num_ghost);
//
//    ExtrapFillPatch filler(m_scalarNew[a_var][lev]->disjointBoxLayout(), b, Interval(1,m_num_ghost));
//
//
//    for (int dir=0; dir<SpaceDim; dir++)
//    {
//      filler.fillExtrap(*m_scalarNew[a_var][lev], dir, 0, 1);
//    }
//  }
//
//  if (a_var == m_porosity)
//  {
//    for (int lev=0; lev<=m_finest_level; lev++)
//      for (DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
//      {
//
//        {
//          (*m_scalarNew[m_solidFraction][lev])[dit].setVal(1.0);
//          (*m_scalarNew[m_solidFraction][lev])[dit] -= (*m_scalarNew[m_porosity][lev])[dit];
//        }
//      }
//  }


}

//Add a transient initial period before we reach steady state porosity that gives channels
Real amrMushyLayer::timeDepPorosity(Real x, Real z, Real t)
{
  Real porosityTimescale = 0.0;
  ParmParse pp("parameters");
  pp.query("porosityTimescale", porosityTimescale);

  //t=t*1e8;
  t=(t+1e-9)/porosityTimescale;

  Real max_x = m_domainSize[0];
  Real chimney_width = min(0.05, pow(log(1+t/15),1.1)); //0.05
  Real chimneyOne = 0.25*max_x;
  Real chimneyTwo = 0.75*max_x;

  Real h0 = 0.25;
  h0 = max(1-sqrt(0.75*t), 0.25); // Grow mushy layer at diffusive timescale

  Real chim_height = 2;

  Real h = h0 - (pow(chimney_width,2)*(1/pow(0.05,2)))*(1/pow(chim_height, 5))*exp(-chim_height*cos(M_PI*x));

  //					Real  vertical = 1-tanh(5* pow((z-h0)/h,3) );
  //					Real vertical = 1-( 1-exp(-pow(z, 2)) )*tanh(5* pow((z-h0)/h,3) );
  Real vertical = 1-( log(3*z+4)/3 )*tanh(5* pow((z-h0)/h,3) );

  Real porosity = vertical;

  Real channelVertical = max(0.0, min(1.0, sqrt(1+2*t-z)-1));
  if (porosityTimescale==0)
  {
    channelVertical=1;
  }

  if (z > 0.8)
  {
    channelVertical = 0;
  }

  porosity = porosity + 1.2*exp(-pow( (x-chimneyOne)/chimney_width,4) )*channelVertical;
  porosity = porosity + 1.2*exp(-pow( (x-chimneyTwo)/chimney_width,4) )*channelVertical;


  porosity = min(porosity, min(0.99, 0.5+t/2));
  porosity = max(porosity, 1e-4);

  //porosity = pow(porosity, 0.5);
  porosity = pow(porosity, 0.5);

  //porosity = 0.2;


  //						porosity = 1.0;

  return porosity;

}
