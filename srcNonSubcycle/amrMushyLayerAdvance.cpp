#include "amrMushyLayer.H"

#include "computeNorm.H"
#include "Gradient.H"
#include "Divergence.H"
#include "AMRTGA.H"
#include "CH_assert.H"
#include "MultilevelLinearOp.H"
#include "RelaxSolver.H"
#include "CellToEdge.H"
#include "ExtrapFillPatch.H"
#include "VelBCHolder.H"
#include "PiecewiseLinearFillPatchFace.H"
#include "SetValLevel.H"
#include "EdgeToCell.H"
#include "CoarseAverageFace.H"


void amrMushyLayer::enforcePorosity()
{
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    LevelData<FArrayBox>& levelPhiNew = *(m_scalarNew[m_porosity][lev]);
    DataIterator levelDit = levelPhiNew.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      FArrayBox& porosity = levelPhiNew[levelDit()];
      FArrayBox& solidFrac = (*m_scalarNew[m_solidFraction][lev])[levelDit()];
      Box thisBox = porosity.box();
      IntVect box_iv = thisBox.smallEnd();

      BoxIterator bit(thisBox);
      for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, lev, loc);
        Real x = loc[0];
        Real z = loc[1];

        porosity(iv) = timeDepPorosity(x,z,m_time);
        solidFrac(iv) = 1-porosity(iv);
      }
    }
  }
}

void
amrMushyLayer::timeStep()
{
  CH_TIME("amrMushyLayer::timestep()");

  if (m_parameters.physicalProblem == m_enforcedPorosity)
  {
    enforcePorosity();
  }

  if (s_verbosity >=2)
  {
    pout() << "Timestep " << m_cur_step
        << " Advancing solution from time "
        << m_time << " with dt = " << m_dt << endl;
    pout() << "Finest level = " << m_finest_level << endl;
  }
  // copy new data to old data
  for (int a_var=0; a_var<m_numVars; a_var++)
  {
    assert (m_scalarOld[a_var].size() == m_scalarNew[a_var].size());
    activeLevelCopy(m_scalarNew[a_var], m_scalarOld[a_var]);
  }

  //Generate vectors of grids and variables which only span the active levels
  //This is the data structure that AMRElliptic expects
  Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > a_scalarNew, a_scalarOld;

  Vector<DisjointBoxLayout> amrGrids;
  getActiveGridsVars(a_scalarNew, a_scalarOld, amrGrids);

  //We need a different data structure for the BE/TGA solver (this is a bit silly)
  Vector<LevelData<FArrayBox>* >  a_thetaNew, a_thetaOld, a_thetaDiff, a_thetaOldTemp, a_ThetaLNew, a_ThetaLOld;

  //Our guess of theta^(n+1) at the previous iteration, used to check for convergence of the iteration
  Vector<LevelData<FArrayBox>* > a_thetaNewPrev;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > Theta_prev, H_prev;

  /// values at n-2 timestep, for predictor in non-linear update
  Vector<RefCountedPtr<LevelData<FArrayBox> > > H_twoPrev, Theta_twoPrev;


  a_thetaNew.resize(m_finest_level+1, NULL);
  a_thetaOld.resize(m_finest_level+1, NULL);
  a_thetaNewPrev.resize(m_finest_level+1, NULL);
  Theta_prev.resize(m_finest_level+1);
  H_prev.resize(m_finest_level+1);
  Theta_twoPrev.resize(m_finest_level+1);
  H_twoPrev.resize(m_finest_level+1);

  a_ThetaLNew.resize(m_finest_level+1, NULL);
  a_ThetaLOld.resize(m_finest_level+1, NULL);

  //Keep these in case they're useful for more debugging
  a_thetaDiff.resize(m_finest_level+1, NULL);
  a_thetaOldTemp.resize(m_finest_level+1, NULL);

  IntVect zeroGhost = IntVect::Zero;

//  for (int lev=0; lev<=m_finest_level; lev++)
//  {
//    a_thetaNew[lev] = &(*a_scalarNew[m_theta][lev]);
//    a_thetaOld[lev] = &(*a_scalarOld[m_theta][lev]);
//
//    a_thetaNewPrev[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_ghostVect);
//    a_thetaNew[lev]->copyTo(*a_thetaNewPrev[lev]);
//
//    a_ThetaLNew[lev] = &(*a_scalarNew[m_compositionLiquid][lev]);
//    a_ThetaLOld[lev] = &(*a_scalarOld[m_compositionLiquid][lev]);
//
//    Theta_prev[lev] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(m_amrGrids[lev], 1, zeroGhost));
//    m_scalarNew[m_bulkConcentration][lev]->copyTo(*Theta_prev[lev]);
//
//    H_prev[lev] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(m_amrGrids[lev], 1, zeroGhost));
//    m_scalarNew[m_HC][lev]->copyTo(*H_prev[lev]);
//
//    H_twoPrev[lev] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(m_amrGrids[lev], 1, zeroGhost));
//    m_scalarNew[m_HC][lev]->copyTo(*H_twoPrev[lev]);
//
//    Theta_twoPrev[lev] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(m_amrGrids[lev], 1, zeroGhost));
//    m_scalarNew[m_bulkConcentration][lev]->copyTo(*Theta_twoPrev[lev]);
//
//  }


  logMessage(8, "    amrMushyLayer::timestep() - setup solvers");


  Vector<LevelData<FArrayBox>* > HC_src(m_finest_level+1,NULL);
  Vector<LevelData<FArrayBox>* > ThetaLSource(m_finest_level+1,NULL);
  Vector<LevelData<FArrayBox>* > zeroSource(m_finest_level+1,NULL);
  Vector<LevelData<FArrayBox>* > ThetaDiffusion(m_finest_level+1,NULL);
  Vector<LevelData<FArrayBox>* > V_gradH_n(m_finest_level+1,NULL);

  Vector<LevelData<FArrayBox>* > V_dThetadz_n(m_finest_level+1,NULL);
  Vector<LevelData<FArrayBox>* > ThetaDiffusion_n(m_finest_level+1,NULL);


  for (int lev=0; lev<=m_finest_level; lev++)
  {
    const DisjointBoxLayout& grids = m_amrGrids[lev];
    HC_src[lev] = new LevelData<FArrayBox>(grids, 1,m_ghostVect);
    ThetaLSource[lev] = new LevelData<FArrayBox>(grids, 1,m_ghostVect);
    zeroSource[lev] = new LevelData<FArrayBox>(grids, 1,m_ghostVect);
    V_gradH_n[lev] = new LevelData<FArrayBox>(grids, 1,m_ghostVect);
    V_dThetadz_n[lev] = new LevelData<FArrayBox>(grids, 1,m_ghostVect);
    ThetaDiffusion_n[lev] = new LevelData<FArrayBox>(grids, 1,m_ghostVect);

    DataIterator dit = HC_src[lev]->dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      (*HC_src[lev])[dit()].setVal(0);
//      (*ThetaLSource[lev])[dit()].setVal(0);
//      (*zeroSource[lev])[dit()].setVal(0);
//      (*V_gradH_n[lev])[dit()].setVal(0);
//      (*V_dThetadz_n[lev])[dit()].setVal(0);
//      (*ThetaDiffusion_n[lev])[dit()].setVal(0);
//
//      (*a_thetaNewPrev[lev])[dit].setVal(0);
//      (*a_thetaNewPrev[lev])[dit] += (*a_thetaNew[lev])[dit];

    }
  }


  //Setup solvers for this timestep
  // Not doing advection for now
//  setupAdvectionSolvers();

  // This is consumed into multigrid setup
  //  setupPoissonSolvers(a_thetaNew, thetaSource, a_ThetaLNew, zeroSource, a_amrGrids);

  //Must be done after setupPoissonSolvers()
  updateEnthalpyVariables();
  setupMultigrid(amrGrids);

  //Do fluid advection
  //Get fluid velocities based on our current best estimate of the
  //fields at time n+1/2
  /*  if (m_parameters.physicalProblem != m_bmDiffusiveSolidification && m_parameters.physicalProblem != m_bmDiffusion)
  {
    Real timeInterpolant = 0.5;
    updateVelocityComp(a_amrGrids, timeInterpolant);
    makeAdvectionTerm(m_theta, zeroSource, m_fluidAdv, m_dScalar[m_theta]);
  }
  else
  {
    for (int lev=0; lev<=m_finest_level; lev++)
    {
      for (DataIterator dit=m_dScalar[m_theta][lev]->dataIterator(); dit.ok(); ++dit)
      {
        (*m_dScalar[m_theta][lev])[dit].setVal(0.0);
      }
    }
  }

  ///Calculate V*dH/dz at time n
  if (m_frameAdvectionMethod == m_finiteDifference)
  {
    calculateFrameAdvection(m_HC, V_gradH_n, m_parameters.nonDimVel);
    calculateFrameAdvection(m_composition, V_dThetadz_n, m_parameters.nonDimVel);
  }
  else if(m_frameAdvectionMethod == m_godunov)
  {
    calculateFrameAdvectionGodunov(m_HC);
    calculateFrameAdvectionGodunov(m_composition);
  }
  else
  {
    MayDay::Error("Unrecognized m_frameAdvectionMethod specified");
  }*/

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    for (DataIterator dit=m_dScalar[m_theta][lev]->dataIterator(); dit.ok(); ++dit)
    {
      (*HC_src[lev])[dit].setVal(0.0);
    }
  }


  logMessage(8, "amrMushyLayer::timestep() - begin heat solve");

  int lbase = 0;
  int lmax = m_finest_level;

  if (m_timeIntegrationOrder == 1)
  {
    m_BEEnthalpySalinity->oneStep(m_HC,
                                  m_HC,
                              HC_src,
                              m_dt,
                              lbase,
                              lmax,
                              false); //don't set phi=zero
  }
  else if(m_timeIntegrationOrder == 2)
  {
    m_TGAEnthalpySalinity->setTime(m_time);
    m_TGAEnthalpySalinity->resetAlphaAndBeta(1,1); //think this might be unnecessary

    m_TGAEnthalpySalinity->oneStep(a_thetaNew,
                               a_thetaOld,
                               HC_src,
                               m_dt,
                               lbase,
                               lmax,
                               m_time);
  }
  else
  {
    MayDay::Error("amrMushyLayer::timestep() - Invalid m_timeIntegrationOrder");
  }


  updateEnthalpyVariables();



  //Post timestep operations:


  //Do refluxing if we've been doing advection
  //				if (m_parameters.problemType != m_bmDiffusion and m_parameters.problemType != m_bmDiffusiveSolidification)
  //				{
  //					for (int lev=0; lev<m_finest_level; lev++)
  //					{
  //						doReflux(m_compositionLiquid, lev);
  //					}
  //
  //				}

  //For debugging
  //Turn off to speed up run time
  //	calculateThetaDiffusion(ThetaDiffusion_n);

  // average solutions down to coarser levels
  averageCoarseToFineSolutions();


  // finally, update to new time and increment current step
  m_time += m_dt;
  m_cur_step += 1;

  // clean up temp storage

  for (int lev=0; lev<=m_finest_level; lev++)
  {

    if (a_thetaNewPrev[lev] != NULL)
    {
      delete a_thetaNewPrev[lev];
      a_thetaNewPrev[lev] = NULL;
    }
    if (HC_src[lev] != NULL)
    {
      delete HC_src[lev];
      HC_src[lev] = NULL;
    }
    if (zeroSource[lev] != NULL)
    {
      delete zeroSource[lev];
      zeroSource[lev] = NULL;
    }
    if (ThetaLSource[lev] != NULL)
    {
      delete ThetaLSource[lev];
      ThetaLSource[lev] = NULL;
    }
    if (ThetaDiffusion[lev] != NULL)
    {
      delete ThetaDiffusion[lev];
      ThetaDiffusion[lev] = NULL;
    }
    if (V_gradH_n[lev] != NULL)
    {
      delete V_gradH_n[lev];
      V_gradH_n[lev] = NULL;
    }
    if (V_dThetadz_n[lev] != NULL)
    {
      delete V_dThetadz_n[lev];
      V_dThetadz_n[lev] = NULL;
    }
    if (ThetaDiffusion_n[lev] != NULL)
    {
      delete ThetaDiffusion_n[lev];
      ThetaDiffusion_n[lev] = NULL;
    }
  }

  pout () << "amrMushyLayer::timestep " << m_cur_step
      << ",      end time = "
      << setiosflags(ios::fixed) << setprecision(6) << setw(12) << m_time << "( " << m_time*m_parameters.timescale << " secs) "
      << ", dt = "
      << setiosflags(ios::fixed) << setprecision(6) << setw(12) << m_dt << "( " << m_dt * m_parameters.timescale << " secs) "
      << endl;


  int totalCellsAdvanced = 0;
  for (int lev=0; lev<m_num_cells.size(); lev++)
  {
    totalCellsAdvanced += m_num_cells[lev];
  }

  //	pout() << "Time = " << m_time << " cells advanced = "
  //			<< totalCellsAdvanced << endl;
  for (int lev=0; lev<m_num_cells.size(); lev++)
  {
    //		pout () << "Time = " << m_time
    //				<< "  level " << lev << " cells advanced = "
    //				<< m_num_cells[lev] << endl;
  }


}

void amrMushyLayer::
doReflux(int a_var, int refluxLevel)
{
  if (m_parameters.physicalProblem == m_fullProblem or m_parameters.physicalProblem == m_bmHRL)
  {
    doImplicitReflux(a_var, refluxLevel);
  }
  else if (m_parameters.physicalProblem == m_bmAdvection)
  {
    // explicit Reflux
    Real scale = -1.0/m_amrDx[refluxLevel];
    m_fluxRegister[a_var][refluxLevel]->reflux(*m_scalarNew[a_var][refluxLevel],scale);
  }
}

void
amrMushyLayer::
doImplicitReflux(int a_var, int refluxLevel)
{

  CH_assert(refluxLevel < m_finest_level);

  // now do implicit refluxing
  // Vector of pointers to LevelData of FABS
  Vector<LevelData<FArrayBox>* > refluxCor(m_finest_level+1, NULL);
  Vector<LevelData<FArrayBox>* > refluxRHS(m_finest_level+1, NULL);
  // collectUN: AMR vector containing soln at new time
  Vector<LevelData<FArrayBox>* > collectUN(m_finest_level+1, NULL);

  // loop over levels, allocate storage, set up for AMRSolve
  // if coarser level exists, define it as well for BCs.
  int startlev = Max(refluxLevel-1, 0);

  for (int ilev = startlev; ilev<= m_finest_level; ilev++)
  {
    // rhs has no ghost cells, correction does
    refluxRHS[ilev]  = new LevelData<FArrayBox>(m_amrGrids[ilev], 1, IntVect::Zero);
    refluxCor[ilev]  = new LevelData<FArrayBox>(m_amrGrids[ilev], 1, IntVect::Unit);
    collectUN[ilev]  = &(*m_scalarNew[a_var][ilev]);
    for (DataIterator dit = m_amrGrids[ilev].dataIterator(); dit.ok(); ++dit)
    {
      (*refluxRHS[ilev])[dit()].setVal(0.0);
      (*refluxCor[ilev])[dit()].setVal(0.0);
    }
  } // end loop over levels for setup.

  // now loop over levels and set up RHS
  // note that we don't need to look at finest level here,
  // since it will have no finer level to get a reflux correction
  // from.   Also this starts at refluxLevel instead of startLev since,
  // in the case refluxLevel > 0, refluxLevel-1 is only used for boundary conditions
  for (int lev= refluxLevel; lev < m_finest_level; lev++)
  {
    //TODO: should this be negative or positive?
    Real dxLev = m_amrDx[lev];
    Real refluxScale = 1.0/dxLev;

    m_fluxRegister[a_var][lev]->reflux(*refluxRHS[lev], refluxScale);
  }


  int lbase = refluxLevel;
  int lmax  = m_finest_level;

  // this resets the coeffients including eta, alpha, beta
  Real alpha = 1.0;
  Real beta  = -m_dt;

  Vector<MGLevelOp<LevelData<FArrayBox> >* > ops = m_refluxAMRMG->getAllOperators();
  for (int iop = 0; iop < ops.size(); iop++)
  {
    LevelTGAHelmOp<FArrayBox,FluxBox>* helmop = (LevelTGAHelmOp<FArrayBox,FluxBox>*) ops[iop];
    helmop->setAlphaAndBeta(alpha, beta);
  }

  m_refluxAMRMG->solve(refluxCor, refluxRHS, lmax, lbase);

  for (int lev= refluxLevel; lev <= m_finest_level; lev++)
  {
    for (DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    {
      (*m_scalarNew[a_var][lev])[dit()] += (*refluxCor[lev])[dit()];
    }
  }

  //remember that startlev can be different from refluxLevel
  for (int lev = startlev; lev<= m_finest_level; lev++)
  {
    delete refluxRHS[lev];
    delete refluxCor[lev];
    refluxRHS[lev] = NULL;
    refluxCor[lev] = NULL;
  }
}

void amrMushyLayer::
calculatePermeability(bool oldTime)
{
  CH_TIME("amrMushyLayer::calculatePermeability()");

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    for (DataIterator dit=m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    {


      if (oldTime)
      {
        ::calculatePermeability((*m_scalarOld[m_permeability][lev])[dit],
                                (*m_scalarOld[m_solidFraction][lev])[dit],
                                m_parameters, m_amrDx[lev]);
      }
      else
      {
        ::calculatePermeability((*m_scalarNew[m_permeability][lev])[dit],
                                (*m_scalarNew[m_solidFraction][lev])[dit],
                                m_parameters, m_amrDx[lev]);
      }


    } // end data iterator
  }
}
//


void amrMushyLayer::
doImplicitVelocityReflux(int m_level)
{

  // loop over levels and compute RHS
  Vector<LevelData<FArrayBox>* > refluxRHS(m_finest_level+1,NULL);
  Vector<LevelData<FArrayBox>* > refluxCorr(m_finest_level+1,NULL);
  // this is necessary because while solver only can do
  // one component, levelfluxRegister::reflux can only
  // do all of them at once.
  Vector<LevelData<FArrayBox>* > tempRefluxData(m_finest_level+1,NULL);
  // loop over levels, allocate storage, set up for AMRMultiGrid
  // solve

  Real refluxScale;

  int startLev=m_level;
  // if crser level exists, define it as well for BC's
  if (startLev > 0)
  {
    startLev = startLev-1;
  }

  for (int lev=startLev; lev<=m_finest_level; lev++)
  {
    //		const DisjointBoxLayout& levelGrids = thisNSPtr->newVel().getBoxes();
    const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

    // recall that AMRMultiGrid can only do one component
    // rhs has no ghost cells
    refluxRHS[lev] = new LevelData<FArrayBox>(levelGrids, 1);
    tempRefluxData[lev] = new LevelData<FArrayBox>(levelGrids,
                                                   SpaceDim);
    //soln has one layer of ghost cells
    IntVect ghostVect(D_DECL(1,1,1));
    refluxCorr[lev] = new LevelData<FArrayBox>(levelGrids,1,
                                               ghostVect);

    // initialize rhs to 0
    DataIterator levelDit = tempRefluxData[lev]->dataIterator();
    LevelData<FArrayBox>& levelRefluxData = *(tempRefluxData[lev]);
    for (levelDit.reset(); levelDit.ok(); ++levelDit)
    {
      levelRefluxData[levelDit()].setVal(0.0);
    }

    // while we're here, do refluxing.
    // (recall that startLev may be coarser than m_level
    // for BC's,
    // however, don't want to do anything to that level)
    if  ((lev >= m_level) && (lev < m_finest_level))
    {
      refluxScale = -1.0/m_amrDx[lev]; // petermc made negative, 7 Dec 07
      m_vectorFluxRegister[m_fluidVel][lev]->reflux(levelRefluxData,
                                                    refluxScale);
    }

    //		int temp=0;

  } // end loop over levels

  // coarse BC is 0
  if (m_level > 0)
  {
    LevelData<FArrayBox>& crseBCData = *refluxCorr[m_level-1];
    DataIterator crseDit = crseBCData.dataIterator();
    for (crseDit.reset(); crseDit.ok(); ++crseDit)
    {
      crseBCData[crseDit()].setVal(0.0);
    }
  }

  // set convergence metric to be norm of velocity
  // for now just do component-wise maxnorm
  // compute norm over all directions and then use for
  // all velocity components (in order to be consistent)
  //	Vector<LevelData<FArrayBox>* > vectVel(m_finest_level+1, NULL);
  //	thisNSPtr = this;
  //	for (int lev=m_level; lev<=m_finest_level; lev++)
  //	{
  //		vectVel[lev] = thisNSPtr->m_vel_new_ptr;
  //		if (lev < m_finest_level)
  //		{
  //			thisNSPtr = thisNSPtr->finerNSPtr();
  //		}
  //	}

  int normType = 0;
  Interval allVelComps(0, SpaceDim-1);

  Real velNorm = computeNorm(m_vectorNew[m_fluidVel], m_refinement_ratios,
                             m_amrDx[m_level], allVelComps, normType, m_level);

  // now loop over directions
  Interval solverComps(0,0);
  for (int dir=0; dir<SpaceDim; dir++)
  {
    Interval velComps(dir,dir);
    for (int lev=m_level; lev<=m_finest_level; ++lev)
    {
      // copy rhs to single-component holder
      tempRefluxData[lev]->copyTo(velComps,*(refluxRHS[lev]),
                                  solverComps);
      // initial guess for correction is RHS.
      LevelData<FArrayBox>& levelCorr = *(refluxCorr[lev]);
      DataIterator levelDit = levelCorr.dataIterator();
      for (levelDit.reset(); levelDit.ok(); ++levelDit)
      {
        levelCorr[levelDit()].setVal(0.0);
      }
      refluxRHS[lev]->copyTo(solverComps,*(refluxCorr[lev]),
                             solverComps);
    } // end set initial guess

    // now set up solver
    int numLevels = m_finest_level+1;

    // This is a Helmholtz operator
    Real alpha = 1.0;
    /* viscosity - set to 0 as we don't (currently) have a viscous term in our momentum equation
     * \mathbf{U} = - \nabla P + (Ra_T \theta) \mathbf{k}
     */
    Real nu = 0.0;

    Real beta = -nu*m_dt;

    AMRPoissonOpFactory viscousOpFactory;
    viscousOpFactory.define(m_amrDomains[0],
                            m_amrGrids,
                            m_refinement_ratios,
                            m_amrDx[0],
                            m_physBCPtr->viscousRefluxBC(dir),
                            alpha,
                            beta);

    RelaxSolver<LevelData<FArrayBox> > bottomSolver;
    bottomSolver.m_verbosity = s_verbosity;

    AMRMultiGrid<LevelData<FArrayBox> > viscousSolver;
    AMRLevelOpFactory<LevelData<FArrayBox> >& viscCastFact = (AMRLevelOpFactory<LevelData<FArrayBox> >&) viscousOpFactory;
    viscousSolver.define(m_amrDomains[0],
                         viscCastFact,
                         &bottomSolver,
                         numLevels);

    viscousSolver.m_verbosity = s_verbosity;

    //		int  s_viscous_solver_type = 2;
    Real s_viscous_solver_tol = 1e-7;
    int s_viscous_num_smooth_up = 4;
    int s_viscous_num_smooth_down = 4;

    viscousSolver.m_eps = s_viscous_solver_tol;

    viscousSolver.m_pre  = s_viscous_num_smooth_down;
    viscousSolver.m_post = s_viscous_num_smooth_up;

    // This needs to be fixed when AMRMultiGrid has this
    // capability.
    //
    // viscousSolver.setConvergenceMetric(velNorm, 0);
    viscousSolver.m_convergenceMetric = velNorm;

    // now solve
    viscousSolver.solve(refluxCorr, refluxRHS,
                        m_finest_level, m_level,
                        false); // don't initialize to zero

    // now increment velocity with reflux correction
    for (int lev=m_level; lev<=m_finest_level; ++lev)
    {
      //			LevelData<FArrayBox>& levelVel = *(m_vectorNew[m_fluidVel][lev]);
      LevelData<FluxBox>& levelVel = *(m_fluidAdv[lev]);
      LevelData<FArrayBox>& levelCorr = *(refluxCorr[lev]);
      DataIterator levelDit = levelCorr.dataIterator();
      for (levelDit.reset(); levelDit.ok(); ++levelDit)
      {
        FArrayBox& levelVelDir = levelVel[levelDit][dir];
        levelVelDir.plus(levelCorr[levelDit()],0,dir,1);
      }
    }
  } // end loop over directions

  // clean up storage
  for (int lev=startLev; lev<=m_finest_level; lev++)
  {
    if (refluxRHS[lev] != NULL)
    {
      delete refluxRHS[lev];
      refluxRHS[lev] = NULL;
    }

    if (refluxCorr[lev] != NULL)
    {
      delete refluxCorr[lev];
      refluxCorr[lev] = NULL;
    }

    if (tempRefluxData[lev] != NULL)
    {
      delete tempRefluxData[lev];
      tempRefluxData[lev] = NULL;
    }
  }
}



void amrMushyLayer::
updateVelocityComp(Vector<DisjointBoxLayout> activeGrids, Real alpha)
{

  logMessage(8, "    amrMushyLayerAdvance::updateVelocity");
  // Need the latest permeability in order to calculate U^*
  calculatePermeability(); //Calculate permeability at new timestep

  if (alpha < 1)
  {
    calculatePermeability(true);//Calculate permeability at old timestep
  }

  //For debugging
  Vector<RefCountedPtr<LevelData<FArrayBox> > > zFluidVel;
  zFluidVel.resize(m_finest_level+1);

  PiecewiseLinearFillPatch fillPatch;

  logMessage(10, "    amrMushyLayerAdvance::updateVelocity - apply BCs and fill ghost cells");
  bool isViscous = false;
  //	VelBCHolder velBC(m_physBCPtr->velFuncBC(isViscous));
  VelBCHolder velBCExtrap(m_physBCPtr->velExtrapBC(isViscous));
  //	EdgeVelBCHolder edgeVelExtrapBC(m_physBCPtr->velExtrapBC(isViscous));
  EdgeVelBCHolder edgeVelBC(m_physBCPtr->edgeVelFuncBC(isViscous));

  //	VelBCHolder velInteriorExtrap(m_physBCPtr->velExtrapInterior(isViscous, 4));

  //Calculate the predicted velocity, U^*
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    zFluidVel[lev] = RefCountedPtr<LevelData<FArrayBox> >
    (new LevelData<FArrayBox>(m_amrGrids[lev], 1,m_ghostVect));
    m_vectorNew[m_fluidVel][lev]->copyTo(Interval(1,1), *zFluidVel[lev], Interval(0,0));

    if (lev > 0)
    {
      fillPatch.define(m_amrGrids[lev], m_amrGrids[lev-1], SpaceDim, m_amrDomains[lev-1], getRefinementRatio(lev-1), m_num_ghost);
    }

    LevelData<FArrayBox>& uStar = *m_vectorNew[m_fluidVelPred][lev];
    LevelData<FArrayBox>& uAnalytic = *m_vectorNew[m_fluidVelAnalytic][lev];
    LevelData<FArrayBox>& gradPErr = *m_vectorNew[m_gradPressureErr][lev];
    LevelData<FArrayBox>& divUstarAnalytic = *m_scalarNew[m_divUstarErr][lev];
    LevelData<FArrayBox>& pressureErr = *m_scalarNew[m_pressureError][lev];
    LevelData<FArrayBox>& thetaAdvectionAnalytic = *m_scalarNew[m_thetaForcingAnalytic][lev];

    logMessage(10, "    amrMushyLayerAdvance::updateVelocity - make U^*");
    for (DataIterator dit = m_vectorNew[m_fluidVelPred][lev]->dataIterator(); dit.ok(); ++dit)
    {
      logMessage(10, "    amrMushyLayerAdvance::updateVelocity - make U^* on grid");

      Box b = uStar[dit].box();

      FArrayBox thetaInit(b, 1);

      for (BoxIterator bit(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, lev, loc);

        Real permeability, theta, Theta;

        //Time average permeability
        Real permeabilityNew = (*m_scalarNew[m_permeability][lev])[dit](iv,0);
        Real permeabilityOld = (*m_scalarOld[m_permeability][lev])[dit](iv,0);


        Real thetaNew = (*m_scalarNew[m_theta][lev])[dit](iv,0);
        Real thetaOld = (*m_scalarOld[m_theta][lev])[dit](iv,0);


        Real ThetaNew = (*m_scalarNew[m_bulkConcentration][lev])[dit](iv,0);
        Real ThetaOld = (*m_scalarOld[m_bulkConcentration][lev])[dit](iv,0);

        permeability = (1-alpha) * permeabilityOld + alpha * permeabilityNew;
        theta = (1-alpha) * thetaOld + alpha * thetaNew;
        Theta = (1-alpha) * ThetaOld + alpha * ThetaNew;


        Real x = (iv[0]+0.5)*m_amrDx[lev];
        Real y = (iv[1]+0.5)*m_amrDx[lev];

        //				uStar[dit](iv, 0) = 0;
        //				uStar[dit](iv, 1) = 0;
        //
        //				if (iv[0] == 10 && iv[0] == 10)
        //				{
        //					uStar[dit](iv, 0) = 1;
        //				}
        uStar[dit](iv, 0) = 0;
        uStar[dit](iv, 1) = 0;

        //				if (iv[0] < 30  && iv[0] > 5 &&
        //					iv[1] < 30  && iv[1] > 5 )
        //				{
        //					uStar[dit](iv, 1) = theta;
        //				}

        /*
         * Test case for BCs:
         * Actually now dp/dx = 0 (x=0,1), dp/dz = 1 (z=0,1)
         * P = cos(pi*x) * cos(pi*z) analytic solution
         */
        uStar[dit](iv, 0) = - M_PI * cos(M_PI*y)*sin(M_PI*x);
        uStar[dit](iv, 1) =  - M_PI * sin(M_PI*y)*cos(M_PI*x);
        divUstarAnalytic[dit](iv, 0) = -2*M_PI*M_PI*cos(M_PI*y)*cos(M_PI*x);

        pressureErr[dit](iv,0) = cos(M_PI*x)*cos(M_PI*y) + 0.0232083;

        gradPErr[dit](iv,0) = -M_PI*sin(M_PI*x)*cos(M_PI*y);
        gradPErr[dit](iv,1) = -M_PI*cos(M_PI*x)*sin(M_PI*y);

        uAnalytic[dit](iv, 0) = 0;
        uAnalytic[dit](iv, 1) =  0;

        /*
         * Test case for dp/dx = 0 (x=0,1) and dp/dz=-pi*cos(pi*x) (z=0,1) BCs
         *cos(M_PI*x)				 * P = cos(pi*x) * sin(pi*z) analytic solution
         */
        //				uStar[dit](iv, 0) = - M_PI * sin(M_PI*x)*sin(M_PI*y);
        //				uStar[dit](iv, 1) =  M_PI * cos(M_PI*y)*cos(M_PI*x);

        /*
         * Test case for zero divergence
         */
        //				uStar[dit](iv, 0) = cos(M_PI*4*x)*cos(M_PI*y);
        //				uStar[dit](iv, 1) = 4*sin(M_PI*y)*sin(M_PI*4*x);

        /*
         * Actual problem
         * U^*_x = 0
         * U^*_z = permeability * (Ra_t * theta + Ra_C * Theta)
         */

        //Perturbation comes from init
        Real perturbation = 0.1;

        uStar[dit](iv, 0)  = 0;
        uStar[dit](iv, 1) = permeability * ( m_parameters.rayleighTemp * theta +
            m_parameters.rayleighComposition * Theta);

        /*
         * Create fields of analytic solutions
         */

        divUstarAnalytic[dit](iv, 0) = m_parameters.rayleighTemp *(-1+perturbation*M_PI*cos(M_PI*(x))*cos(M_PI*(y)));

        pressureErr[dit](iv,0) = m_parameters.rayleighTemp *(y-0.5*y*y - (perturbation/(2*M_PI))*cos(M_PI*(x))*cos(M_PI*(y)));

        // Have to be careful here as we're comparing against edge centred data.
        x = x - m_amrDx[lev]/2;
        gradPErr[dit](iv,0) = m_parameters.rayleighTemp *(0.5*perturbation*sin(M_PI*(x))*cos(M_PI*(y)));

        x = x + m_amrDx[lev]/2; y = y - m_amrDx[lev]/2;
        gradPErr[dit](iv,1) = m_parameters.rayleighTemp *(1-y+0.5*perturbation*cos(M_PI*(x))*sin(M_PI*(y)));

        y = y + m_amrDx[lev]/2;



        //				thetaAdvectionError[dit](iv, 0) = m_parameters.rayleighTemp * (-1 + cos(M_PI*x) *
        //						(-M_PI * alpha * cos(M_PI*y) + 0.5  * alpha*sin(M_PI*y) )
        //						+ (M_PI/4) * alpha*alpha*sin(2*M_PI*y) );
        thetaAdvectionAnalytic[dit](iv, 0) = - m_parameters.rayleighTemp * perturbation * (
            - 0.5 * cos(M_PI*x) + 0.5 * M_PI * perturbation * cos(M_PI*y)) * sin(M_PI*y) ;



      } //finish loop over intvects in box

      logMessage(10, "    amrMushyLayerAdvance::updateVelocity - finished make U^* on grid");



    } // finish loop over boxes

    //Apply extrap BCs to U^* to ensure that div(U^*) is calculated correctly at the boundaries
    velBCExtrap.applyBCs(uStar, activeGrids[lev],
                         m_amrDomains[lev], m_amrDx[lev],
                         false); // inhomogeneous

  } // end loop over levels to setup U^*

  logMessage(10, "    amrMushyLayerAdvance::updateVelocity - do projection");
  // Project the predicted velocity onto the space of divergence free vectors

  CCProjectorComp proj;
  proj.define(activeGrids, m_amrDomains, m_amrDx, m_refinement_ratios, m_finest_level,
              *m_physBCPtr, m_scalarNew[m_theta], m_scalarNew[m_permeability], m_num_ghost);

  proj.projectVelocity(m_vectorNew[m_fluidVel], m_scalarNew[m_permeability],
                       m_fluidAdv, m_vectorNew[m_fluidVelPred], m_fluidVelDiffOrder);

  logMessage(10, "    amrMushyLayerAdvance::updateVelocity - apply corrections");

  // Synchronization step
  //	logMessage(2, " updateVelocity - do synchronisation");

  // 1. Do velocity refluxing
  for (int lev=0; lev <= m_finest_level; lev++)
  {
    //		doImplicitVelocityReflux(lev);
    //		doExplicitVelocityReflux();
  }

  // 2. Ensure velocity is divergence free

  // we should not need a lambda correction here?
  // Set lambda = 1, advect it, subtract 1. We should get zero - anything else is an error.
 /* Vector<LevelData<FArrayBox>* > zeroSource = newEmptyPtr();

  for (int lev=0; lev <= m_finest_level; lev++)
  {
    //		setValLevel(*m_scalarNew[m_lambda][lev], 1.0);
    for (DataIterator dit = m_scalarNew[m_lambda][lev]->dataIterator(); dit.ok(); ++dit)
    {
      (*m_scalarNew[m_lambda][lev])[dit].setVal(1.0);
      (*m_scalarOld[m_lambda][lev])[dit].setVal(1.0);
      (*m_scalarNew[m_lambdaPostCorr][lev])[dit].setVal(1.0);
      (*m_scalarOld[m_lambdaPostCorr][lev])[dit].setVal(1.0);
      (*zeroSource[lev])[dit].setVal(0.0);
    }
  }

  advectScalar(m_lambda, zeroSource, m_fluidAdv);

  Vector<LevelData<FArrayBox>* > a_lambda(m_finest_level+1, NULL);
  refCountedPtrToPtr(m_scalarNew[m_lambda], a_lambda);

  proj.doLambdaCorrection(m_fluidAdv, a_lambda, m_time, m_dt);*/

  //Re apply BCs after correction
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    edgeVelBC.applyBCs(*m_fluidAdv[lev], activeGrids[lev],
                       m_amrDomains[lev], m_amrDx[lev],
                       false); // inhomogeneous
  }

//  advectScalar(m_lambdaPostCorr, zeroSource, m_fluidAdv);

//  for (int lev=0; lev <= m_finest_level; lev++)
//  {
//    for (DataIterator dit = m_scalarNew[m_lambdaPostCorr][lev]->dataIterator(); dit.ok(); ++dit)
//    {
//      (*m_scalarNew[m_lambda][lev])[dit].plus(-1);
//      (*m_scalarNew[m_lambdaPostCorr][lev])[dit].plus(-1);
//    }
//  }


  Vector<LevelData<FluxBox>* > gradPressureEdge;
  gradPressureEdge.resize(m_finest_level+1, NULL);
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    gradPressureEdge[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, m_ghostVect - IntVect::Unit);
  }

  proj.getDivUStar(m_scalarNew[m_divUstar]);
  proj.getPressure(m_scalarNew[m_pressure]);
  proj.getGradPressure(m_vectorNew[m_gradPressure]);
  proj.getGradPressureEdge(gradPressureEdge);



  //Apply BCs
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    edgeVelBC.applyBCs(*m_fluidAdv[lev], activeGrids[lev],
                       m_amrDomains[lev], m_amrDx[lev],
                       false); // inhomogeneous
  }

  logMessage(10, "    amrMushyLayerAdvance::updateVelocity - calculate errors");



  // Need to clean up memory
  proj.undefine();


  //CF Interp test
  doCFInterpTest();

  //Calculate div(U)
  LevelData<FluxBox>* uEdgeFinePtr = NULL;
  Real* dxFine = NULL;
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    if (lev < m_finest_level)
    {
      uEdgeFinePtr = m_fluidAdv[lev+1];
      dxFine = &(m_amrDx[lev+1]);
    }
    else
    {
      uEdgeFinePtr = NULL;
      dxFine = NULL;
    }
    Divergence::compDivergenceMAC(*m_scalarNew[m_divU][lev], *m_fluidAdv[lev], uEdgeFinePtr,
                                  m_amrDx[lev],
                                  dxFine,
                                  getRefinementRatio(lev),
                                  m_amrDomains[lev]);
  }


  logMessage(10, "    amrMushyLayerAdvance::updateVelocity - Finished all levels");


  // Clean up memory
  for (int lev=0; lev<=m_finest_level; lev++){
    if(gradPressureEdge[lev]!=NULL)
    {
      delete gradPressureEdge[lev];
      gradPressureEdge[lev]  = NULL;
    }


  }

}

void amrMushyLayer::
doCFInterpTest()
{

  //Set up initial field
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    for (DataIterator dit=m_scalarNew[m_CFInterpTest][lev]->dataIterator(); dit.ok(); ++dit)
    {
      Box b = (*m_scalarNew[m_CFInterpTest][lev])[dit].box();
      FArrayBox& CFInterpTest = (*m_scalarNew[m_CFInterpTest][lev])[dit];
      for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, lev, loc);
        Real x = loc[0]; Real y = loc[1];
        CFInterpTest(iv,0) = sin(M_PI*x)*sin(M_PI*y);
      }
    }
  }

  // Do CF interp
  for (int lev=1; lev<=m_finest_level; lev++)
  {
    QuadCFInterp CFInterp(m_scalarNew[m_CFInterpTest][lev]->disjointBoxLayout(), &(m_scalarNew[m_CFInterpTest][lev-1]->disjointBoxLayout()),
                          m_amrDx[lev], getRefinementRatio(lev-1),
                          1, m_amrDomains[lev]);

    CFInterp.coarseFineInterp(*m_scalarNew[m_CFInterpTest][lev], *m_scalarNew[m_CFInterpTest][lev-1]);
  }

  // Calculate error
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    for (DataIterator dit=m_scalarNew[m_CFInterpTest][lev]->dataIterator(); dit.ok(); ++dit)
    {
      Box b = (*m_scalarNew[m_CFInterpTest][lev])[dit].box();
      FArrayBox& CFInterpTest = (*m_scalarNew[m_CFInterpTest][lev])[dit];
      for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, lev, loc);
        Real x = loc[0]; Real y = loc[1];
        CFInterpTest(iv,0) -= sin(M_PI*x)*sin(M_PI*y);
      }

      if (lev > 0)
      {
        //				int temp=0;
      }
    }
  }
}





Real amrMushyLayer::
calculateMaxDiff(Vector<LevelData<FArrayBox> * > a_phiNew, Vector<LevelData<FArrayBox> * > a_phiOld)
{
  Vector<RefCountedPtr<LevelData<FArrayBox> > > refPtrNew, refPtrOld;

  ptrToRefCountedPtr(a_phiNew, refPtrNew);
  ptrToRefCountedPtr(a_phiOld, refPtrOld);

  return calculateMaxDiff(refPtrNew, refPtrOld);
}

Real amrMushyLayer::
calculateDiffAtPoint(Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phiNew,
                     Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phiOld,
                     int z_i)
{
  Real diff = 0.0;
  for(int lev=0; lev<=m_finest_level; lev++)
  {
    for(DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    {
      for (BoxIterator bit = BoxIterator(m_amrGrids[lev][dit]); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        if (iv[1] == z_i)
        {
          Real phiNewVal =  (*a_phiNew[lev])[dit](iv,0);
          Real phiOldVal =  (*a_phiOld[lev])[dit](iv,0);
          diff = phiNewVal - phiOldVal;
        }
      }

    }
  }
  return diff;

}
Real amrMushyLayer::
calculateMaxDiff(Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phiNew, Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phiOld)
{
  Vector<LevelData<FArrayBox>* > diff(m_finest_level+1, NULL);

  //Calculate difference between two consecutive guesses of theta
  for(int lev=0; lev<=m_finest_level; lev++)
  {
    diff[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_ghostVect);
    a_phiNew[lev]->copyTo(*diff[lev]);

    for(DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    {
      (*diff[lev])[dit] -= (*a_phiOld[lev])[dit];
    }
  }

  Interval comps(0,0);
  Real norm = computeNorm(diff, m_refinement_ratios, m_amrDx[0], comps, 0);

  return norm;
}




void amrMushyLayer::
applyBCs(const int a_var, const int lev)
{
  CH_TIME("amrMushyLayer::applyBCs");
  const DisjointBoxLayout& dbl = (*m_scalarNew[a_var][lev]).disjointBoxLayout();
  bool a_homogeneous = false;

  BCHolder bcHolder;

  //Get the BC function for this variable.
//  if (a_var==m_HC)
//  {
//    bcHolder = m_physBCPtr->enthalpySalinityBC(a_homogeneous);
//  }
//  else if(a_var == m_composition)
//  {
//    bcHolder = m_physBCPtr->ThetaFuncBC();
//  }
//  else if(a_var == m_porosity)
//  {
//    bcHolder = m_physBCPtr->BasicPorosityFuncBC(a_homogeneous);
//    //				ABParseBCPorosity( (*m_scalarNew[a_var][lev])[dit], dbl[dit], m_amrDomains[lev], m_amrDx[lev], a_homogeneous);
//  }
//  else if(a_var == m_theta)
//  {
//    bcHolder = m_physBCPtr->BasicthetaFuncBC(a_homogeneous, NULL);
//  }
//  else if(a_var == m_compositionLiquid)
//  {
//    bcHolder = m_physBCPtr->ThetaLFuncBC();
//  }
//  else if(a_var == m_compositionSolid)
//  {
//    bcHolder = m_physBCPtr->BasicScalarFuncBC();
//  }
//  else
//  {
    //let's not throw an error anymore
    MayDay::Error("applyBCs() - Can't calculate apply bcs for a_var specified");
    return;
//  }



  DataIterator dit = m_amrGrids[lev].dataIterator();
  {
    for (dit.begin(); dit.ok(); ++dit)
    {
      Box b = dbl[dit];

      bcHolder.operator()((*m_scalarNew[a_var][lev])[dit],
                          b,
                          m_amrDomains[lev],
                          m_amrDx[lev],
                          a_homogeneous);

      bcHolder.operator()((*m_scalarOld[a_var][lev])[dit],
                          b,
                          m_amrDomains[lev],
                          m_amrDx[lev],
                          a_homogeneous);



      //Grow box b to include the ghost cells we want to fill
      b.grow(m_num_ghost);
      //Interval is the grid cells we want to interp into
      //0 is the ghost cell adjacent to boundary, which we already have
      ExtrapFillPatch patch(dbl, b, Interval(1,m_num_ghost));
      for (int dir=0; dir<SpaceDim; dir++)
      {
        patch.fillExtrap((*m_scalarNew[a_var][lev]), dir, 0, 1);
        patch.fillExtrap((*m_scalarOld[a_var][lev]), dir, 0, 1);
      }

    }
  }

}

void amrMushyLayer::applyBCs(int a_var)
{
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    applyBCs(a_var, lev);
  }
}

void amrMushyLayer::
applyVectorBCs(const int a_var, const int lev)
{
  CH_TIME("amrMushyLayer::applyVectorBCs");
  const DisjointBoxLayout& dbl = (*m_vectorNew[a_var][lev]).disjointBoxLayout();
  //	bool a_homogeneous = false;

  VelBCHolder bcHolder;

  //Get the BC function for this variable.
  if (a_var==m_fluidVel)
  {
    bcHolder = m_physBCPtr->uStarFuncBC(false);
  }
  else
  {
    //let's not throw an error anymore
    //				MayDay::Error("applyBCs() - Can't calculate apply bcs for a_var specified");
    return;
  }

  bcHolder.applyBCs(*m_vectorNew[a_var][lev], m_amrGrids[lev],
                    m_amrDomains[lev], m_amrDx[lev],
                    false);


  DataIterator dit = m_amrGrids[lev].dataIterator();
  {

    Box b = m_amrGrids[lev][dit];
    //Grow box b to include the ghost cells we want to fill
    b.grow(m_num_ghost);
    //Interval is the grid cells we want to interp into
    //0 is the ghost cell adjacent to boundary, which we already have
    ExtrapFillPatch patch(dbl, b, Interval(1,m_num_ghost));
    for (int dir=0; dir<SpaceDim; dir++)
    {
      patch.fillExtrap((*m_vectorNew[a_var][lev]), dir, 0, SpaceDim);
      patch.fillExtrap((*m_vectorOld[a_var][lev]), dir, 0, SpaceDim);
    }


  }

}

void amrMushyLayer::
calculateFrameAdvectionGodunov(const int a_var)
{
  Vector<LevelData<FArrayBox>* > zeroSource = newEmptyPtr();

  advectScalar(a_var, zeroSource, m_frameAdv);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    for (DataIterator dit = m_dScalar[a_var][lev]->dataIterator(); dit.ok(); ++dit)
    {
      (*m_dScalar[a_var][lev])[dit].setVal(0);
      (*m_dScalar[a_var][lev])[dit] += (*m_scalarNew[a_var][lev])[dit];
      (*m_dScalar[a_var][lev])[dit] -= (*m_scalarOld[a_var][lev])[dit];
      (*m_dScalar[a_var][lev])[dit] /= m_dt;
    }

    //Reset scalar new to undo the effects of advectScalar()
    m_scalarOld[a_var][lev]->copyTo(*m_scalarNew[a_var][lev]);
  }
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    if(zeroSource[lev]!=NULL)
    {
      delete zeroSource[lev];
      zeroSource[lev] = NULL;
    }
  }
}


//This is a new method to explicitly calculate V*d(phi)/dz, as the
//method using a godunov update appears to be first order
void amrMushyLayer::
calculateFrameAdvection(const int a_var, Vector<LevelData<FArrayBox>* >& a_d_dz, Real frameAdvVel)
{
  CH_TIME("amrMushyLayer::calculateFrameAdvection()");

  Vector<LevelData<FArrayBox>* > gradient(m_finest_level+1,NULL);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    const DisjointBoxLayout& grids = m_amrGrids[lev];
    gradient[lev] = new LevelData<FArrayBox>(grids, 2,m_ghostVect);
    DataIterator dit = gradient[lev]->dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      (*gradient[lev])[dit()].setVal(0);
    }
  }

  for (int lev=0; lev<=m_finest_level; lev++)
  {

    //		for (DataIterator dit = m_dScalar[a_var][lev]->dataIterator(); dit.ok(); ++dit)
    //		{
    //			FArrayBox& box = (*m_scalarNew[a_var][lev])[dit];
    //			int bal = 1;
    //		}

    //Make sure ghost cells are filled correctly before calculating gradient
    m_scalarNew[a_var][lev]->exchange();

    //Asumes BCS are set
    //Gradient::levelGradientCC(*gradient[lev], *m_scalarNew[a_var][lev], m_amrDx[lev]);

    //The upwinded version ensures we don't consider downwind information, which shouldn't be affecting
    //the solution
    Gradient::levelGradientCCUpwind(*gradient[lev], *m_scalarNew[a_var][lev], m_amrDx[lev]);

    for (DataIterator dit = m_dScalar[a_var][lev]->dataIterator(); dit.ok(); ++dit)
    {

      //multiply by non dimensional velocity
      (*gradient[lev])[dit].mult(frameAdvVel);

      //Take z-component of gradient and put in our dvar/dz
      BoxIterator bit(m_amrGrids[lev][dit]);
      for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        (*a_d_dz[lev])[dit](iv, 0) += (*gradient[lev])[dit](iv, 1);
      }
    }
  }

  //Clean up memory
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    if (gradient[lev] != NULL)
    {
      delete gradient[lev];
      gradient[lev] = NULL;
    }
  }
}

Vector<LevelData<FArrayBox>* > amrMushyLayer::
newEmptyPtr()
{

  Vector<LevelData<FArrayBox>* > phi(m_finest_level+1,NULL);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    const DisjointBoxLayout& grids = m_amrGrids[lev];
    phi[lev] = new LevelData<FArrayBox>(grids, 1,m_ghostVect);

    DataIterator dit = phi[lev]->dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      (*phi[lev])[dit()].setVal(0);
    }
  }

  return phi;
}
Vector<RefCountedPtr<LevelData<FArrayBox> > > amrMushyLayer::
newEmptyRefPtr()
{

  Vector<RefCountedPtr<LevelData<FArrayBox> > > phi(m_finest_level+1);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    const DisjointBoxLayout& grids = m_amrGrids[lev];
    phi[lev] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(grids, 1,m_ghostVect));

    DataIterator dit = phi[lev]->dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      (*phi[lev])[dit()].setVal(0);
    }
  }

  return phi;
}

void amrMushyLayer::
calculateThetaLSourceTerm(Vector<LevelData<FArrayBox>* > V_dThetadz_n, Vector<LevelData<FArrayBox>* > ThetaLSource)
{
  CH_TIME("amrMushLayer:calculateThetaLSourceTerm");
  //Best guess at frame advection at n+1
  Vector<LevelData<FArrayBox>* > V_dThetadz_nPlusOne = newEmptyPtr();

  //Theta_l^(n+1/2)
  Vector<LevelData<FArrayBox>* > averageThetaL = timeCenteredScalar(m_compositionLiquid);

  //Solid concentration source term
  Vector<LevelData<FArrayBox>* > ThetaSSource = newEmptyPtr();


  if (m_frameAdvectionMethod == m_finiteDifference)
  {
    calculateFrameAdvection(m_bulkConcentration, V_dThetadz_nPlusOne, m_parameters.nonDimVel);
  }

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    for (DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    {
      (*m_dScalar[m_bulkConcentration][lev])[dit].setVal(0);

      // Frame advection, averaged over times n+1 and n
      if (m_frameAdvectionMethod == m_finiteDifference)
      {
        (*m_dScalar[m_bulkConcentration][lev])[dit].setVal(0);
        (*m_dScalar[m_bulkConcentration][lev])[dit] += (*V_dThetadz_n[lev])[dit];
        (*m_dScalar[m_bulkConcentration][lev])[dit] += (*V_dThetadz_nPlusOne[lev])[dit];
        (*m_dScalar[m_bulkConcentration][lev])[dit] /= 2;
      }

      //Theta_s source term
      (*m_dScalar[m_compositionSolid][lev])[dit].setVal(1);
      (*m_dScalar[m_compositionSolid][lev])[dit].minus((*m_scalarNew[m_porosity][lev])[dit]);
      (*m_dScalar[m_compositionSolid][lev])[dit].mult((*m_scalarNew[m_compositionSolid][lev])[dit]);

      (*ThetaSSource[lev])[dit].setVal(1);
      (*ThetaSSource[lev])[dit].minus((*m_scalarOld[m_porosity][lev])[dit]);
      (*ThetaSSource[lev])[dit].mult((*m_scalarOld[m_compositionSolid][lev])[dit]);

      (*m_scalarNew[m_ThetaSSource][lev])[dit].setVal(0);
      (*m_scalarNew[m_ThetaSSource][lev])[dit] += (*m_dScalar[m_compositionSolid][lev])[dit];
      (*m_scalarNew[m_ThetaSSource][lev])[dit] -= (*ThetaSSource[lev])[dit];
      (*m_scalarNew[m_ThetaSSource][lev])[dit].divide(m_dt);
      (*m_scalarNew[m_ThetaSSource][lev])[dit].mult(-1);

      //Porosity source term
      (*m_scalarNew[m_ThetaPorositySource][lev])[dit].setVal(0);
      (*m_scalarNew[m_ThetaPorositySource][lev])[dit] += (*m_scalarNew[m_porosity][lev])[dit];
      (*m_scalarNew[m_ThetaPorositySource][lev])[dit] -= (*m_scalarOld[m_porosity][lev])[dit];
      (*m_scalarNew[m_ThetaPorositySource][lev])[dit] /= m_dt;
      (*m_scalarNew[m_ThetaPorositySource][lev])[dit].mult((*averageThetaL[lev])[dit]);
      (*m_scalarNew[m_ThetaPorositySource][lev])[dit].mult(-1);

      //Put it all together, praying that all the signs are correct
      (*ThetaLSource[lev])[dit].setVal(0);
      (*ThetaLSource[lev])[dit] += (*m_dScalar[m_bulkConcentration][lev])[dit];
      (*ThetaLSource[lev])[dit] += (*m_scalarNew[m_ThetaPorositySource][lev])[dit];
      (*ThetaLSource[lev])[dit] += (*m_scalarNew[m_ThetaSSource][lev])[dit];

      //For debugging
      (*m_scalarNew[m_concSource][lev])[dit].setVal(0);
      (*m_scalarNew[m_concSource][lev])[dit] += (*ThetaLSource[lev])[dit];

    }
    ThetaLSource[lev]->exchange();
  }

  //Clean up memory
  for (int lev=0; lev <= m_finest_level; lev++)
  {
    if (V_dThetadz_nPlusOne[lev] != NULL)
    {
      delete V_dThetadz_nPlusOne[lev];
      V_dThetadz_nPlusOne[lev] = NULL;
    }
    if (averageThetaL[lev] != NULL)
    {
      delete averageThetaL[lev];
      averageThetaL[lev] = NULL;
    }
    if (ThetaSSource[lev] != NULL)
    {
      delete ThetaSSource[lev];
      ThetaSSource[lev] = NULL;
    }

  }


}


//Calculate (1/Le) * grad (porosity* grad Theta_l)



void amrMushyLayer::makeAdvectionTerm(const int a_var,
                                      const Vector<LevelData<FArrayBox>* > a_source,
                                      const Vector<LevelData<FluxBox>* > a_advVel,
                                      Vector<RefCountedPtr<LevelData<FArrayBox> > > a_advTerm)
{
  Interval interv(0,0);

  //First get a backup of m_scalarNew, and put m_scalarOld into m_scalarNew so it can be advected
  Vector<RefCountedPtr<LevelData<FArrayBox> > > scalarNewBackup;
  scalarNewBackup.resize(m_finest_level+1);
  for (int lev=0; lev<=m_finest_level; lev++)
  {
    scalarNewBackup[lev] = RefCountedPtr<LevelData<FArrayBox> >
    (new LevelData<FArrayBox>(m_amrGrids[lev], 1,m_ghostVect));

    m_scalarNew[a_var][lev]->copyTo(interv, *scalarNewBackup[lev], interv);
    m_scalarOld[a_var][lev]->copyTo(interv, *m_scalarNew[a_var][lev], interv);
  }

  advectScalar(a_var, a_source, a_advVel);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    //TODO: can we remove these exchange calls?
    m_scalarNew[a_var][lev]->exchange();
    m_scalarOld[a_var][lev]->exchange();

    m_scalarNew[a_var][lev]->copyTo(interv, *a_advTerm[lev], interv);

    //for (DataIterator dit = m_dScalar[a_var][lev]->dataIterator(); dit.ok(); ++dit)
    for (DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    {
      (*a_advTerm[lev])[dit] -= (*m_scalarOld[a_var][lev])[dit];
      (*a_advTerm[lev])[dit] /= m_dt;
    }
    a_advTerm[lev]->exchange();
    //Put things back where they should be
    //			m_scalarOld[a_var][lev]->copyTo(interv, *m_scalarNew[a_var][lev], interv);
    scalarNewBackup[lev]->copyTo(interv, *m_scalarNew[a_var][lev], interv);
  }

}


void
amrMushyLayer::advectScalar(int a_var,
                            Vector<LevelData<FArrayBox>* > a_source,
                            Vector<LevelData<FluxBox>* > a_advVel)
{
  for (int lev = 0; lev<=m_finest_level; lev++)
  {

    RefCountedPtr<LevelFluxRegister> coarserFRPtr, finerFRPtr;
    RefCountedPtr<LevelData<FArrayBox> > coarserDataOldPtr, coarserDataNewPtr;
    Real tCoarserNew, tCoarserOld;

    LevelAdvect levAdvect;
    DisjointBoxLayout coarseGrid = DisjointBoxLayout();
    int refRatio; //These live on the coarser level
    bool use_limiting = true;
    bool hasCoarser = (lev > 0);
    bool hasFiner = (lev < m_finest_level);

    ParmParse pp("godunov");
    pp.query("use_limiting", use_limiting);

    getLevelAdvectionsVars(a_var, lev, refRatio, hasCoarser, hasFiner, use_limiting,
                           coarserDataOldPtr, coarserDataNewPtr,
                           coarserFRPtr, finerFRPtr, coarseGrid, levAdvect,
                           tCoarserNew, tCoarserOld);

    m_scalarNew[a_var][lev]->exchange();

    levAdvect.step(*m_scalarNew[a_var][lev],
                   *finerFRPtr,
                   *coarserFRPtr,
                   *a_advVel[lev],
                   *a_source[lev],
                   *coarserDataOldPtr,
                   tCoarserOld,
                   *coarserDataNewPtr,
                   tCoarserNew,
                   m_time,
                   m_dt);

  }


}

void
amrMushyLayer::getActiveGridsVars(Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > >& a_scalarNew,
                                  Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > >& a_scalarOld,
                                  Vector<DisjointBoxLayout>&a_amrGrids)
{
  CH_TIME("amrMushyLayer::getActiveGridsVars");
  a_scalarNew.resize(m_numVars);
  a_scalarOld.resize(m_numVars);
  a_amrGrids.resize(m_finest_level+1);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    a_amrGrids[lev] = m_amrGrids[lev];
  }

  for (int a_var=0; a_var<m_numVars; a_var++)
  {
    a_scalarNew[a_var].resize(m_finest_level+1);
    a_scalarOld[a_var].resize(m_finest_level+1);

    for (int lev=0; lev<=m_finest_level; lev++)
    {
      //These are just pointers
      a_scalarNew[a_var][lev] = m_scalarNew[a_var][lev];
      a_scalarOld[a_var][lev] = m_scalarOld[a_var][lev];
    }
  }
}


void
amrMushyLayer::getFrameAdvection(const int lev)
{
  //Delete the old data to clean up memory
  if (m_frameAdv[lev] != NULL)
  {
    delete m_frameAdv[lev];
    m_frameAdv[lev] = NULL;
  }

  // Now create new data
  m_frameAdv[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, m_num_ghost*IntVect::Unit);

  for (DataIterator dit = m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
  {
    // The current box
    Box curBox = m_amrGrids[lev].get(dit());

    FluxBox& advVelFluxBox = (*m_frameAdv[lev])[dit];
    for (int faceDir=0; faceDir < SpaceDim; faceDir++)
    {
      Real vel;
      if (faceDir == 1)
      {
        // Not the minus sign here
        // Fluid is pulled down towards the cold heat exchanger
        vel = - m_parameters.nonDimVel;
      }
      else
      {
        vel = 0;
      }
      FArrayBox& velDir = advVelFluxBox[faceDir];
      velDir.setVal(vel);
    }
  }
}


//a_coef = 1
//b_coef = permeability
void amrMushyLayer::
getPressureVCcoefs(Vector<RefCountedPtr<LevelData<FArrayBox> > >& aCoef,
                   Vector<RefCountedPtr<LevelData<FluxBox> > >& bCoef)
{
  aCoef.resize(m_finest_level+1);
  bCoef.resize(m_finest_level+1);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >
    (new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_ghostVect ));

    bCoef[lev] = RefCountedPtr<LevelData<FluxBox> >
    (new LevelData<FluxBox>(m_amrGrids[lev], 1, m_ghostVect));

    bCoef[lev].neverDelete();

    LevelData<FluxBox> *fluxBox = &(*bCoef[lev]);

    //Create FArray box of a and b coefficients
    for (DataIterator dit=m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    {
      (*aCoef[lev])[dit].setVal(0);

      FluxBox& fbox = (*fluxBox)[dit];
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
      {
        FArrayBox& fluxDir = fbox[faceDir];
        fluxDir.setVal(-1);
      }
    }

    aCoef[lev]->exchange();
    bCoef[lev]->exchange();

    //		//Turn FArrayBox into FluxBox for B coefficient
    //LevelData<FluxBox> *fluxBox = &(*bCoef[lev]);
    CellToEdge(*m_scalarNew[m_permeability][lev], *fluxBox);


    // Multiply bCoef fluxbox by -1 so signs are correct
    for (DataIterator dit=m_amrGrids[lev].dataIterator(); dit.ok(); ++dit)
    {
      (*fluxBox)[dit] *= -1;
    }

  }
}




void
amrMushyLayer::getFluidAdvection(const int lev)
{
  logMessage(10, "    amrMushyLayerAdvance::getFluidAdvection");
  //Delete the old data to clean up memory
  if (m_fluidAdv[lev] != NULL)
  {
    delete m_fluidAdv[lev];
    m_fluidAdv[lev] = NULL;
  }

  // Now create new data
  m_fluidAdv[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, m_ghostVect);

  //Set everything to zero
  for (DataIterator dit = m_vectorNew[m_fluidVel][lev]->dataIterator(); dit.ok(); ++dit)
  {
    FluxBox& temp_advVel = (*m_fluidAdv[lev])[dit];
    for (int faceDir=0; faceDir < SpaceDim; faceDir++)
    {
      Real vel = 0;

      FArrayBox& velDir = temp_advVel[faceDir];
      velDir.setVal(vel);
    }
    //		int temp=0;
  }

  // Get velocities from the FArrayBox
  m_vectorNew[m_fluidVel][lev]->exchange();
  CellToEdge(*m_vectorNew[m_fluidVel][lev], *m_fluidAdv[lev]);

  m_fluidAdv[lev]->exchange();

  //Do my own quadratic interpolation!!

  for (DataIterator dit = m_vectorNew[m_fluidVel][lev]->dataIterator(); dit.ok(); ++dit)
  {
    //			Box b =  (*m_vectorNew[m_fluidVel][lev])[dit].box();
    Box b = m_amrGrids[lev][dit];

    FluxBox& edgeData  = (*m_fluidAdv[lev])[dit];

    for (int idir=0; idir<SpaceDim; idir++)
    {
      b.growHi(idir, -1);

      FArrayBox& fluidAdvDir =edgeData[idir];

      for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        Real y0 = (*m_vectorNew[m_fluidVel][lev])[dit](iv - BASISV(idir), idir);
        Real y1 = (*m_vectorNew[m_fluidVel][lev])[dit](iv, idir);
        Real y2 = (*m_vectorNew[m_fluidVel][lev])[dit](iv + BASISV(idir), idir);

        if (m_fluidVelInterpOrder == 2)
        {
          fluidAdvDir(iv, 0) =  (3*y0 + 6*y1 - y2)/8;  //2nd order
        }
        else if (m_fluidVelInterpOrder == 1)
        {
          fluidAdvDir(iv, 0) =  (y0+y1)/2; //1st order
        }
        else
        {
          fluidAdvDir(iv, 0) =  y1; //0th order
        }

      }

      Box end = adjCellBox(b, idir, Side::Hi , 2);

      for (BoxIterator bit = BoxIterator(end); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        Real y0 = (*m_vectorNew[m_fluidVel][lev])[dit](iv - 2*BASISV(idir), idir);
        Real y1 = (*m_vectorNew[m_fluidVel][lev])[dit](iv - BASISV(idir), idir);
        Real y2 = (*m_vectorNew[m_fluidVel][lev])[dit](iv, idir);

        if (m_fluidVelInterpOrder == 2)
        {
          fluidAdvDir(iv, 0) =  (-y0 + 6*y1 + 3*y2)/8;  //2nd order
        }
        else if (m_fluidVelInterpOrder == 1)
        {
          fluidAdvDir(iv, 0) =  (y1+y2)/2; //1st order
        }
        else
        {
          fluidAdvDir(iv, 0) =  y1; //0th order
        }

      }

      b.growHi(idir, 1);

    }
  }




  //	for (DataIterator dit = m_vectorNew[m_fluidVel][lev]->dataIterator(); dit.ok(); ++dit)
  //		{
  //
  //			FluxBox& advFlux = (*m_fluidAdv[lev])[dit];
  //			FArrayBox& Ux = (*m_vectorNew[m_fluidVel][lev])[dit];
  //			FArrayBox& Uz = (*m_vectorNew[m_fluidVel][lev])[dit];
  //			Uz.copy((*m_vectorNew[m_fluidVel][lev])[dit], 1,0,1);
  //				int temp=0;
  //		}

}



void amrMushyLayer::removeNanValues(FArrayBox& fab)
{
  //todo - fix
  Box b = fab.box();

  //  for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
  //  {
  //    IntVect iv = bit();
  //    if (isnan(fab(iv,0)))
  //    {
  //      fab(iv,0) = 0;
  //    }
  //  }
}
void amrMushyLayer::enforceAbsCap(Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phi, Real a_cap)
{
  //Ensure this is positive
  a_cap = abs(a_cap);

  for (int lev=0; lev<=m_finest_level; lev++)
  {
    for (DataIterator dit=a_phi[lev]->dataIterator(); dit.ok(); ++dit)
    {
      Box b = m_amrGrids[lev][dit];
      for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        Real phiVal = (*a_phi[lev])[dit](iv,0);
        if (abs(phiVal) > a_cap)
        {
          //Need to enforce the cap with the correct sign
          Real sign = 1;
          if (phiVal < 0)
          {
            sign = -1;
          }

          (*a_phi[lev])[dit](iv,0) = sign*a_cap;
        }
        //        else if (isnan(phiVal))
        //        {
        //          (*a_phi[lev])[dit](iv,0) = 0;
        //        }
      }
    }
  }
}

