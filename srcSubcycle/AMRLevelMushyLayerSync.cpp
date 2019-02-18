#include "AMRLevelMushyLayer.H"


void AMRLevelMushyLayer::doMomentumReflux(Vector<LevelData<FArrayBox>*>& compVel)
{
  CH_TIME("AMRLevelMushyLayer::doMomentumReflux");

  int finest_level = getFinestLevel();

  // amr grid info for solvers
  Vector<DisjointBoxLayout> AmrGrids(finest_level + 1);
  Vector<int> AmrRefRatios(finest_level + 1);
  Vector<Real> AmrDx(finest_level + 1);
  ProblemDomain baseDomain;

  // loop over levels, set up for AMRMultiGrid solve
  AMRLevelMushyLayer* thisMLPtr = this;
  int startLev = m_level;
  // if crser level exists, define it as well for BC's
  if (startLev > 0)
  {
    startLev = startLev - 1;
    thisMLPtr = thisMLPtr->getCoarserLevel();
  }
  //      AMRLevelMushyLayer* startLevelPtr = thisMLPtr;

  // (DFM 4/22/15) -- PoissonOpFactory needs grids for all levels
  // so do coarser levels first
  if (startLev > 0)
  {
    AMRLevelMushyLayer* crseLevPtr = thisMLPtr->getCoarserLevel();
    for (int lev = startLev - 1; lev >= 0; lev--)
    {
      const DisjointBoxLayout& levelGrids = crseLevPtr->m_grids;
      AmrGrids[lev] = levelGrids;
      if (lev == 0)
        baseDomain = crseLevPtr->problemDomain();
      AmrRefRatios[lev] = crseLevPtr->refRatio();
      AmrDx[lev] = crseLevPtr->m_dx;
      if (lev > 0)
        crseLevPtr = crseLevPtr->getCoarserLevel();
    }
  } // end if coarser levels exist

  for (int lev = startLev; lev <= finest_level; lev++)
  {
    const DisjointBoxLayout& levelGrids = thisMLPtr->m_grids;
    AmrGrids[lev] = levelGrids;
    AmrRefRatios[lev] = thisMLPtr->refRatio();
    AmrDx[lev] = thisMLPtr->m_dx;
    if (lev == 0)
      baseDomain = thisMLPtr->problemDomain();

    thisMLPtr = thisMLPtr->getFinerLevel();
  } // end loop over levels involved in the solve




  // loop over levels and compute RHS
  Vector<LevelData<FArrayBox>*> refluxRHS(finest_level + 1, NULL);
  Vector<LevelData<FArrayBox>*> refluxCorr(finest_level + 1,
                                           NULL);
  // this is necessary because while solver only can do
  // one component, levelfluxRegister::reflux can only
  // do all of them at once.
  Vector<LevelData<FArrayBox>*> tempRefluxData(finest_level + 1,
                                               NULL);
  // loop over levels, allocate storage, set up for AMRMultiGrid
  // solve
  thisMLPtr = this;
  startLev = m_level;
  // if crser level exists, define it as well for BC's
  if (startLev > 0)
  {
    startLev = startLev - 1;
    thisMLPtr = thisMLPtr->getCoarserLevel();
  }

  for (int lev = startLev; lev <= finest_level; lev++)
  {
    const DisjointBoxLayout& levelGrids = thisMLPtr->m_grids;
    // recall that AMRMultiGrid can only do one component
    // rhs has no ghost cells
    refluxRHS[lev] = new LevelData<FArrayBox>(levelGrids, 1);
    tempRefluxData[lev] = new LevelData<FArrayBox>(levelGrids,
                                                   SpaceDim);
    //soln has one layer of ghost cells
    IntVect ghostVect(D_DECL(1, 1, 1));
    refluxCorr[lev] = new LevelData<FArrayBox>(levelGrids, 1,
                                               ghostVect);

    // initialize rhs to 0
    DataIterator levelDit = tempRefluxData[lev]->dataIterator();
    LevelData<FArrayBox>& levelRefluxData =
        *(tempRefluxData[lev]);
    for (levelDit.reset(); levelDit.ok(); ++levelDit)
    {
      levelRefluxData[levelDit()].setVal(0.0);
    }

    // while we're here, do refluxing.
    // (recall that startLev may be coarser than m_level
    // for BC's,
    // however, don't want to do anything to that level)
    if ((lev >= m_level) && (lev < finest_level))
    {
      Real refluxScale = -1.0 / AmrDx[lev]; // petermc made negative, 7 Dec 07
      //
      thisMLPtr->m_vectorFluxRegisters[VectorVars::m_fluidVel]->reflux(
          levelRefluxData, refluxScale);
    }

    thisMLPtr = thisMLPtr->getFinerLevel();
  } // end loop over levels

  // coarse BC is 0
  if (m_level > 0)
  {
    LevelData<FArrayBox>& crseBCData = *refluxCorr[m_level - 1];
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
  Vector<LevelData<FArrayBox>*> vectVel(finest_level + 1, NULL);
  thisMLPtr = this;
  for (int lev = m_level; lev <= finest_level; lev++)
  {
    vectVel[lev] = &(*thisMLPtr->m_vectorNew[VectorVars::m_fluidVel]);
    if (lev < finest_level)
    {
      thisMLPtr = thisMLPtr->getFinerLevel();
    }
  }

  int normType = 0;
  Interval allVelComps(0, SpaceDim - 1);

  Real velNorm = computeNorm(vectVel, AmrRefRatios, m_dx,
                             allVelComps, normType, m_level);

  // now loop over directions
  Interval solverComps(0, 0);
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    Interval velComps(dir, dir);
    for (int lev = m_level; lev <= finest_level; ++lev)
    {
      // copy rhs to single-component holder
      tempRefluxData[lev]->copyTo(velComps, *(refluxRHS[lev]),
                                  solverComps);
      // initial guess for correction is RHS.
      LevelData<FArrayBox>& levelCorr = *(refluxCorr[lev]);
      DataIterator levelDit = levelCorr.dataIterator();

      for (levelDit.reset(); levelDit.ok(); ++levelDit)
      {
        levelCorr[levelDit()].setVal(0.0);
      }
      refluxRHS[lev]->copyTo(solverComps, *(refluxCorr[lev]),
                             solverComps);
    } // end set initial guess

    // now set up solver
    int numLevels = finest_level + 1;

    // This is a Helmholtz operator
    Real alpha = 1.0; //1.0
    Real beta = -m_parameters.m_viscosityCoeff * m_dt; // -m_parameters.prandtl * m_dt;

    AMRPoissonOpFactory viscousOpFactory;
    viscousOpFactory.define(baseDomain, AmrGrids, AmrRefRatios,
                            AmrDx[0], m_physBCPtr->viscousRefluxBC(dir, m_viscousBCs), alpha,
                            beta);

    RelaxSolver<LevelData<FArrayBox> > bottomSolver;
    bottomSolver.m_verbosity = m_opt.verbosity_multigrid;

    AMRMultiGrid<LevelData<FArrayBox> > viscousSolver;
    AMRLevelOpFactory<LevelData<FArrayBox> >& viscCastFact =
        (AMRLevelOpFactory<LevelData<FArrayBox> >&) viscousOpFactory;
    viscousSolver.define(baseDomain, viscCastFact,
                         &bottomSolver, numLevels);

    viscousSolver.m_verbosity = m_opt.verbosity_multigrid;

    viscousSolver.m_eps = m_opt.viscous_solver_tol;

    viscousSolver.m_pre = m_opt.viscous_num_smooth_down;
    viscousSolver.m_post = m_opt.viscous_num_smooth_up;

    viscousSolver.m_convergenceMetric = velNorm;

    // now solve
    viscousSolver.solve(refluxCorr, refluxRHS, finest_level,
                        m_level, false); // don't initialize to zero

    // now increment velocity with reflux correction
    AMRLevelMushyLayer* amrMLptr = this;
    while(!amrMLptr->finestLevel())
    {
      amrMLptr = amrMLptr->getFinerLevel();
    }

    for (int lev = finest_level; lev >= m_level; --lev)
    {
      LevelData<FArrayBox>& levelVel = *(compVel[lev]);
      LevelData<FArrayBox>& levelCorr = *(refluxCorr[lev]);
      DataIterator levelDit = levelCorr.dataIterator();
      for (levelDit.reset(); levelDit.ok(); ++levelDit)
      {
        levelVel[levelDit()].plus(levelCorr[levelDit()], 0,  dir, 1);
        //        (*amrMLptr->m_vectorNew[m_fluidRefluxCorr])[levelDit].plus(levelCorr[levelDit()], 0, dir, 1);
      }

      // do average down
      //      bool doAvgDown = true;
      //      pp.query("reflux_average_down", doAvgDown);
      //      if (doAvgDown)
      //      {
      //        if (lev < finest_level)
      //        {
      //          AMRLevelMushyLayer& fineML = *(amrMLptr->getFinerLevel());
      //
      //          // quick sanity check
      //          CH_assert (fineML.m_level == lev+1);
      //
      //          CoarseAverage& avgDown = fineML.m_coarseAverageVector;
      //          LevelData<FArrayBox>& fineVel = *fineML.m_vectorNew[VectorVars::m_fluidVel];
      //          avgDown.averageToCoarse(levelVel, fineVel);
      //        }
      //      }

      amrMLptr = amrMLptr->getCoarserLevel();
    }
  } // end loop over directions

  // clean up storage
  for (int lev = startLev; lev <= finest_level; lev++)
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

Real AMRLevelMushyLayer::maxAllowedEta()
{
  Real maxEta = 2*(1-m_computedCFL);
  if (m_opt.maxEta > 0)
  {
    maxEta = m_opt.maxEta;
  }

  if (m_projection.etaLambda() == 0)
  {
    maxEta = 0;
  }

  return maxEta;
}

void AMRLevelMushyLayer::setAdvVelCentering(Real a_fraction)
{
  CH_assert(m_level==0);
  AMRLevelMushyLayer* mlPtr = this;
  while(mlPtr)
  {
    mlPtr->m_adv_vel_centering = a_fraction;
    mlPtr = mlPtr->getFinerLevel();
  }
}

void AMRLevelMushyLayer::setEta(Real a_eta)
{
  CH_assert(m_level==0);

  // Set eta on all projection operators
  AMRLevelMushyLayer* ml = this;
  while(ml)
  {
    ml->m_projection.etaLambda(a_eta);
    ml = ml->getFinerLevel();
  }

}
void AMRLevelMushyLayer::postTimeStep()
{
  CH_TIME("AMRLevelMushyLayer::postTimeStep");

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::postTimeStep - do sync operations on level " << m_level << endl;
  }

  // Get the advection velocity as a CC variable so we can write it out
  EdgeToCell(m_advVel, *m_vectorNew[VectorVars::m_advectionVel]);

  if (m_level == 0)
  {

    // Check if the timestep failed
    // in this case - replace with old data and reset time
    // postTimestep is called from finest levels first, so once we reach level 0
    // we are done with all other levels at this timestep


    // Check every level to see if we failed
    bool timestepFailed = m_timestepFailed;
    AMRLevelMushyLayer* mlPtr = this;
    while(mlPtr)
    {
      timestepFailed = timestepFailed || mlPtr->m_timestepFailed;
      mlPtr = mlPtr->getFinerLevel();
    }

    if (timestepFailed)
    {
      AMRLevelMushyLayer* mlPtr = this;

      // First work out a minimum allowed dt

      while(mlPtr)
      {
        if (mlPtr->getFinerLevel())
        {
          mlPtr = mlPtr->getFinerLevel();
        }
        else
        {
          break;
        }
      }

      // We now have the finest level
//      Real alpha = 0.1; // Some empirically chosen multiplier (must be < 1)
//      Real maxU =  mlPtr->getMaxVelocity();
      //    m_opt.minDt = alpha*m_initial_dt_multiplier*mlPtr->m_dx/maxU;

      if (s_verbosity > 2)
      {
        pout() << "Timestep failed min dt = " << m_opt.minDt << endl;
      }

      // Check if we're going to make dt too small on the finest level, in
      // which case don't change dt anywhere else
      bool halveDt = true;
      if (m_dt/2 < m_opt.minDt)
      {
        halveDt = false;
        timestepFailed = false; // don't halve the timestep
        if (s_verbosity > 1)
        {
          pout() << "AMRLevelMushyLayer::postTimestep() - Not reducing dt as it would break the minimum dt cap" << endl;
        }

      }

      mlPtr = this;

      while(mlPtr)
      {
        mlPtr->m_timestepFailed = timestepFailed; // Store this on every level

        // Try turning this off - just halve dt
        if (m_opt.solverFailRestartMethod == m_restartResetData)
        {
          restartTimestepFromBackup(); // why was this turned off?
        }

        // Always halve dt
        if (halveDt)
        {
          Real oldDt = mlPtr->m_dt;
          Real newDt = oldDt/2;

          mlPtr->m_projection.rescalePressure(oldDt, newDt);
        }

        mlPtr = mlPtr->getFinerLevel();
        //        mlPtr->dt(-mlPtr->m_dt);

      }

      return;
    }


    // Grow m_adv_vel_centering
    m_adv_vel_centering = m_adv_vel_centering*m_opt.adv_vel_centering_growth;
    m_adv_vel_centering = min(m_adv_vel_centering, 0.5);

    // Set the same value on all levels of refinement
    setAdvVelCentering(m_adv_vel_centering);

    if (s_verbosity >= 3)
    {
      pout() << "Set adv vel centering = " << m_adv_vel_centering << endl;
    }


  } // end if level = 0

  if (hasFinerLevel())
  {
    AMRLevelMushyLayer* fineAMRMLPtr = getFinerLevel();

    // create boundary condition object
    VelBCHolder velBC(m_physBCPtr->uStarFuncBC(m_viscousBCs));

    // first do refluxing and avgDown for conservation
    // do momentum refluxing first
    Real refluxScale = -1.0 / m_dx; // petermc made negative, 7 Dec 07
    if (m_opt.reflux_momentum)
    {
      if (s_implicit_reflux)
      {
        // defer this until we're doing sync projection
      } else {
        m_vectorFluxRegisters[VectorVars::m_fluidVel]->reflux(
            *m_vectorNew[VectorVars::m_fluidVel], refluxScale);
      }
    }

    for (int scalarVar = 0; scalarVar < m_numScalarVars;
        scalarVar++)
    {

      fineAMRMLPtr->m_coarseAverageScalar.averageToCoarse(
          *m_scalarNew[scalarVar],
          *(fineAMRMLPtr->m_scalarNew[scalarVar]));

    }

    for (int vectorVar = 0; vectorVar < m_numVectorVars;
        vectorVar++)
    {
      fineAMRMLPtr->m_coarseAverageVector.averageToCoarse(
          *m_vectorNew[vectorVar],
          *(fineAMRMLPtr->m_vectorNew[vectorVar]));
    }

    // now call sync projection if necessary
    // this is a lot of code
    Real crseTime = 0.0;
    if (m_level > 0)
      crseTime = m_coarser_level_ptr->time();

    // do multilevel operations if this is the coarsest level or if
    // coarser level is not at the same time as this level
    if (m_level == 0 || (abs(crseTime - m_time) > TIME_EPS))
    {

      // first, do refluxing for scalars
      if (m_opt.reflux_enthalpy || m_opt.reflux_concentration || m_opt.reflux_lambda)
      {

        doExplicitReflux(ScalarVars::m_lambda);

        // Compute max lambda
        if (m_level==0)
        {

          Vector<int> nRef;

          AMRLevelMushyLayer* ml = this;
          while (ml)
          {
            nRef.push_back(ml->m_ref_ratio);
            ml = ml->getFinerLevel();
          }

          Vector< LevelData<FArrayBox>* > amrLambda;
          amrLambda.resize(nRef.size());
          fillAMRLambda(amrLambda);

          Real maxLambda = ::computeNorm(amrLambda,nRef,m_dx, Interval(0,0), 0, m_level) - 1;
          if (s_verbosity >= 3)
          {
            pout() << "Max(lambda-1) = " << maxLambda << endl;
          }

          //          m_diagnostics.addDiagnostic(m_diagnostics.m_lambda_err, m_time, maxLambda);

          if (m_opt.variable_eta_factor != 1.0)
          {
            Real newEta = m_projection.etaLambda();

            Real maxEta = maxAllowedEta();


            if (maxLambda > m_maxLambda)
            {
              newEta = newEta*m_opt.variable_eta_factor;
            }
            else
            {
              newEta = newEta/m_opt.variable_eta_factor;
            }

            m_maxLambda = maxLambda;

            newEta = max(newEta, m_opt.minEta);
            newEta = min(newEta, maxEta);

            if (s_verbosity >= 3)
            {
              pout() << "New eta = " << newEta << endl;
            }

            setEta(newEta);
          }
        }

        if (m_opt.reflux_enthalpy || m_opt.reflux_concentration)
        {
          Real refluxRHS = 100.0;
          refluxRHS = doHCreflux();
          if (s_verbosity >= 3)
          {
            pout() << "  Sum (reflux RHS) = " << refluxRHS << endl;
          }
        }
        else
        {
          if (s_verbosity >=3 )
          {
            pout() << "  Not doing HC reflux" << endl;
          }
        }

      } // end scalar refluxing


      // Now move onto momentum

      int finest_level = getFinestLevel();


      // now allocate container for composite velocity and lambda
      Vector<LevelData<FArrayBox>* > compVel(finest_level + 1, NULL);
      Vector<LevelData<FArrayBox>* > compLambda(finest_level + 1, NULL);
      Vector<RefCountedPtr<LevelData<FluxBox> > > compPorosityFace(finest_level + 1);
      Vector<RefCountedPtr<LevelData<FArrayBox> > > compPorosity(finest_level + 1);

      if (m_level > 0)
      {
        // coarser level will be coarse-fine BC
        AMRLevelMushyLayer* crseAMRLevelMushyLayerPtr =
            getCoarserLevel();
        DisjointBoxLayout crseGrids = crseAMRLevelMushyLayerPtr->m_grids;

        LevelData<FArrayBox>* crseVelPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim);
        crseAMRLevelMushyLayerPtr->fillVectorField(*crseVelPtr, m_time, m_fluidVel, true);
        compVel[m_level - 1] = crseVelPtr;

        LevelData<FArrayBox>* crseLambdaPtr = new LevelData<FArrayBox>(crseGrids, 1);
        crseAMRLevelMushyLayerPtr->fillScalars(*crseLambdaPtr, m_time, m_lambda, true);
        compLambda[m_level - 1] = crseLambdaPtr;

        RefCountedPtr<LevelData<FluxBox> > crsePorosityFacePtr = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(crseGrids, 1));
        crseAMRLevelMushyLayerPtr->fillScalarFace(*crsePorosityFacePtr, m_time, m_porosity, true);
        compPorosityFace[m_level - 1] = crsePorosityFacePtr;

        RefCountedPtr<LevelData<FArrayBox> > crsePorosityPtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(crseGrids, 1));
        crseAMRLevelMushyLayerPtr->fillScalars(*crsePorosityPtr, m_time, m_porosity, true);
        compPorosity[m_level - 1] = crsePorosityPtr;

      }

      AMRLevelMushyLayer* thisMLPtr = this;
      // now loop over levels and construct composite vel, lambda
      for (int lev = m_level; lev <= finest_level; lev++)
      {
        compVel[lev] = &(*thisMLPtr->m_vectorNew[VectorVars::m_fluidVel]);
        compLambda[lev] = &(*thisMLPtr->m_scalarNew[ScalarVars::m_lambda]);

        RefCountedPtr<LevelData<FluxBox> > porosityFacePtr = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(thisMLPtr->m_grids, 1));
        thisMLPtr->fillScalarFace(*porosityFacePtr, m_time, m_porosity, true);
        compPorosityFace[lev] = porosityFacePtr;

        RefCountedPtr<LevelData<FArrayBox> > porosityPtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(thisMLPtr->m_grids, 1));
        thisMLPtr->fillScalars(*porosityPtr, m_time, m_porosity, true);
        compPorosity[lev] = porosityPtr;

        const DisjointBoxLayout& levelGrids = compVel[lev]->getBoxes();
        velBC.applyBCs(*compVel[lev], levelGrids,
                       thisMLPtr->problemDomain(), thisMLPtr->m_dx, false); // not homogeneous
        thisMLPtr = thisMLPtr->getFinerLevel();
      }



      // before we do sync projection, do implicit refluxing,
      // if appropriate
      if (m_opt.reflux_momentum && s_implicit_reflux)
      {
        doMomentumReflux(compVel);
      } // end implicit reflux

      // Average down
      // do average down, starting from the second finest level

      if (m_opt.refluxAverageDown)
      {
        thisMLPtr = this;
        while (thisMLPtr->hasFinerLevel())
        {
          thisMLPtr = thisMLPtr->getFinerLevel();
        }
        thisMLPtr = thisMLPtr->getCoarserLevel();

        for (int lev = finest_level-1; lev >= m_level; --lev)
        {
          CH_assert(thisMLPtr->m_level == lev);

          LevelData<FArrayBox>& levelVel = *(compVel[lev]);

          AMRLevelMushyLayer& fineML = *(thisMLPtr->getFinerLevel());

          // quick sanity check
          CH_assert (fineML.m_level == lev+1);

          CoarseAverage& avgDown = fineML.m_coarseAverageVector;
          LevelData<FArrayBox>& fineVel = *fineML.m_vectorNew[VectorVars::m_fluidVel];
          avgDown.averageToCoarse(levelVel, fineVel);

          CoarseAverage& scalAvgDown =  fineML.m_coarseAverageScalar;

          scalAvgDown.averageToCoarse(*thisMLPtr->m_scalarNew[ScalarVars::m_enthalpy],
                                      *fineML.m_scalarNew[ScalarVars::m_enthalpy]);
          scalAvgDown.averageToCoarse(*thisMLPtr->m_scalarNew[ScalarVars::m_bulkConcentration],
                                      *fineML.m_scalarNew[ScalarVars::m_bulkConcentration]);

          thisMLPtr->updateEnthalpyVariables();

          thisMLPtr = thisMLPtr->getCoarserLevel();
        }

      }

      // re-apply physical boundary conditions
      thisMLPtr = this;
      for (int lev = m_level; lev < finest_level; lev++)
      {
        LevelData<FArrayBox>& levelVel = *thisMLPtr->m_vectorNew[VectorVars::m_fluidVel];
        const DisjointBoxLayout& levelGrids = levelVel.getBoxes();
        velBC.applyBCs(levelVel, levelGrids, thisMLPtr->problemDomain(),
                       thisMLPtr->m_dx, false); // not homogeneous
        thisMLPtr = thisMLPtr->getFinerLevel();

        // Calculate dSdt after reflux
        LevelData<FArrayBox> temporaryLD;
        thisMLPtr->compute_d_dt(ScalarVars::m_bulkConcentration, temporaryLD);
        thisMLPtr->compute_d_dt(ScalarVars::m_enthalpy, temporaryLD);
      }

      // call projection-based sync operations
      // composition projection and freestream correction computation
      if (m_opt.doSyncOperations)
      {
        m_projection.doSyncOperations(compVel, compLambda, compPorosityFace, compPorosity, m_time, m_dt);
      }

      // re-apply physical boundary conditions
      thisMLPtr = this;
      for (int lev = m_level; lev < finest_level; lev++)
      {
        LevelData<FArrayBox>& levelVel =
            *thisMLPtr->m_vectorNew[VectorVars::m_fluidVel];
        const DisjointBoxLayout& levelGrids = levelVel.getBoxes();
        velBC.applyBCs(levelVel, levelGrids, thisMLPtr->problemDomain(),
                       thisMLPtr->m_dx, false); // not homogeneous
        thisMLPtr = thisMLPtr->getFinerLevel();
      }

    } // end if doing synchronisation stuff



  } // end if there is a finer level
  else
  {
    // in case this is the only level
    // Let's still compute freestream correction
    // If we've coarsened a previous  mesh, we may have some lambda error left over
    // which we need to sort out


    if (m_level == 0 && m_opt.max_possible_level > 0)
    {
      // create boundary condition object and set physical BC's
      VelBCHolder velBC(m_physBCPtr->uStarFuncBC(m_viscousBCs));

      velBC.applyBCs(*m_vectorNew[VectorVars::m_fluidVel], m_grids, m_problem_domain,
                     m_dx, false); // not homogeneous

      Interval velComps = m_vectorNew[VectorVars::m_fluidVel]->interval();
      m_vectorNew[VectorVars::m_fluidVel]->exchange(velComps);

      if (m_opt.computeFreestreamCorrection)
      {
        Vector<LevelData<FArrayBox>* > compVel(1, NULL);
        Vector<LevelData<FArrayBox>* > compLambda(1, NULL);
        Vector<RefCountedPtr<LevelData<FluxBox> > > compPorosityFace(1);
        Vector<RefCountedPtr<LevelData<FArrayBox> > > compPorosity(1);

        fillAMRLambda(compLambda);
        fillAMRVelPorosity(compVel, compPorosityFace,compPorosity);

        m_projection.doSyncOperations(compVel, compLambda, compPorosityFace, compPorosity, m_time, m_dt);
      }

    } // end if level = 0

  } // end if there is or isn't a finer level


  // Make sure temperature, porosity etc. are consistent with the refluxed
  // and averaged down enthalpy and bulk concentration
  updateEnthalpyVariables();

  computeDiagnostics();

  // Could possibly do this in an AMR sense instead

  if (m_opt.computeVorticityStreamFunction)
  {
    computeVorticityStreamfunction();
  }

  getExtraPlotFields();

  // Do some things from level 0 but across entire hierarchy

  if (m_level == 0)
  {
    AMRLevelMushyLayer* AMRmlptr = this;
    while (AMRmlptr != NULL)
    {
      DataIterator dit = (*AMRmlptr->m_vectorNew[VectorVars::m_fluidVelErr]).dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
      {
        (*AMRmlptr->m_vectorNew[VectorVars::m_fluidVelErr])[dit].setVal(0);
        (*AMRmlptr->m_vectorNew[VectorVars::m_fluidVelErr])[dit].plus((*AMRmlptr->m_vectorNew[VectorVars::m_fluidVel])[dit]);
        (*AMRmlptr->m_vectorNew[VectorVars::m_fluidVelErr])[dit].minus((*AMRmlptr->m_vectorNew[VectorVars::m_fluidVelAnalytic])[dit]);

        (*AMRmlptr->m_scalarNew[ScalarVars::m_Terr])[dit].copy((*AMRmlptr->m_scalarNew[ScalarVars::m_temperatureAnalytic])[dit]);
        (*AMRmlptr->m_scalarNew[ScalarVars::m_Terr])[dit].minus((*AMRmlptr->m_scalarNew[ScalarVars::m_temperature])[dit]);
      }

      // Backup data from this timestep
      AMRmlptr->backupTimestep();

      Real maxU = ::computeNorm(*(AMRmlptr->m_vectorNew[VectorVars::m_fluidVel]), NULL, 1 , m_dx, Interval(0, SpaceDim-1), 0);
      if (maxU > 1e10)
      {
        pout() << "WARNING - During PostTimestep,  Max U = " << maxU << endl;
      }


      AMRmlptr = AMRmlptr->getFinerLevel();
    }
  }
}



void AMRLevelMushyLayer::doExplicitReflux(int a_var)
{
  CH_TIME("AMRLevelMushyLayer::doExplicitReflux");

  if (s_verbosity > 5)
  {
    pout () << "AMRLevelMushyLayer::doExplicitReflux, var = " << m_scalarVarNames[a_var] << ", level = " << m_level << endl;
  }
  int finest_level = getFinestLevel();

  Interval solverComps = Interval(0,0);

  Vector<Real> AmrDx(finest_level + 1);

  AMRLevelMushyLayer* thisMLPtr = this;
  int startLev = m_level;
  // if crser level exists, define it as well for BC's
  if (startLev > 0)
  {
    startLev = startLev - 1;
    thisMLPtr = thisMLPtr->getCoarserLevel();
  }
  //  AMRLevelMushyLayer* startLevelPtr = thisMLPtr;

  Vector<LevelData<FArrayBox>*> scalRefluxCorr(finest_level + 1,  NULL);
  Vector<LevelData<FArrayBox>*> scalRefluxRHS(finest_level + 1,   NULL);
  Vector<LevelData<FArrayBox>*> phiOld(finest_level + 1,   NULL);
  Vector<LevelData<FArrayBox>*> phiNew(finest_level + 1,   NULL);
  //  Vector<LevelData<FArrayBox>*> fullRHS(finest_level + 1,   NULL);

  // startLev is either m_level, or m_level-1 if that is defined
  for (int lev = startLev; lev <= finest_level; lev++)
  {
    AmrDx[lev] = thisMLPtr->m_dx;
    const DisjointBoxLayout& levelGrids = thisMLPtr->m_grids;
    // recall that AMRMultiGrid can only do one component.
    // rhs has no ghost cells
    scalRefluxRHS[lev] = new LevelData<FArrayBox>(levelGrids, 1);

    //soln has one layer of ghost cells
    IntVect ghostVect(D_DECL(1, 1, 1));
    scalRefluxCorr[lev] = new LevelData<FArrayBox>(levelGrids,  1, ghostVect);

    // Make sure we put refluxed solution back into here
    phiNew[lev] = &(*(thisMLPtr->m_scalarNew[a_var]));

    // Will later copy scalarNew into this variable
    phiOld[lev] = &(*(thisMLPtr->m_scalarOld[a_var]));
    //    fullRHS[lev] = new LevelData<FArrayBox>(phiNew[lev]->disjointBoxLayout(), 1);

    // initialize corr, RHS to 0
    DataIterator levelDit = scalRefluxRHS[lev]->dataIterator();
    LevelData<FArrayBox>& levelRefluxRHS = *(scalRefluxRHS[lev]);
    LevelData<FArrayBox>& levelRefluxCorr =  *(scalRefluxCorr[lev]);
    for (levelDit.reset(); levelDit.ok(); ++levelDit)
    {
      levelRefluxRHS[levelDit()].setVal(0.0);
      levelRefluxCorr[levelDit()].setVal(0.0);
    }
    thisMLPtr = thisMLPtr->getFinerLevel();
  }

  // loop over levels and establish RHS, corr
  thisMLPtr = this;

  for (int lev = m_level; lev <= finest_level; lev++)
  {
    if (lev < finest_level)
    {
      LevelData<FArrayBox>& levelRefluxRHS = *(scalRefluxRHS[lev]);
      // first set to 0
      DataIterator levelDit =  levelRefluxRHS.dataIterator();

      for (levelDit.reset(); levelDit.ok(); ++levelDit)
      {
        levelRefluxRHS[levelDit()].setVal(0.0);
      }
      Real refluxScale = -1.0 / AmrDx[lev]; // petermc made negative, 7 Dec 07

      thisMLPtr->m_fluxRegisters[a_var]->reflux(
          levelRefluxRHS, refluxScale);

      // initial guess for correction is RHS
      levelRefluxRHS.copyTo(solverComps,
                            *(scalRefluxCorr[lev]), solverComps);
    }

    thisMLPtr = thisMLPtr->getFinerLevel();
  } // end loop over levels to compute RHS

  // Do averaging down on this level and finer
  thisMLPtr = this;
  for (int lev=m_level; lev<finest_level; lev++)
  {

    AMRLevelMushyLayer& fineML = *(thisMLPtr->getFinerLevel());
    // quick sanity check
    CH_assert(fineML.m_level = lev + 1);
    CoarseAverage& scalAvgDown =fineML.m_coarseAverageScalar;

    scalAvgDown.averageToCoarse(*(scalRefluxRHS[lev]),
                                *(scalRefluxRHS[lev+1]));

    thisMLPtr = thisMLPtr->getFinerLevel();
  }

  // if this component is not diffusive, just add RHS
  // to scalar now increment scalars with reflux
  // correction
  // go from finest->coarsest so that we can also avgDown
  // increment NS pointer to finest level
  thisMLPtr = this;
  while (!thisMLPtr->finestLevel())
  {
    thisMLPtr = thisMLPtr->getFinerLevel();
  }
  CH_assert(thisMLPtr->m_level == finest_level);
  for (int lev = finest_level; lev >= m_level; lev--)
  {
    LevelData<FArrayBox>& levelScal =  *thisMLPtr->m_scalarNew[a_var];
    LevelData<FArrayBox>& levelCorr = *(scalRefluxRHS[lev]);
    DataIterator levelDit = levelCorr.dataIterator();
    for (levelDit.reset(); levelDit.ok(); ++levelDit)
    {
      levelScal[levelDit()] += levelCorr[levelDit()];
    }

    // if a finer level exists, do avgDown as well
    if (lev < finest_level)
    {
      AMRLevelMushyLayer& fineML =
          *(thisMLPtr->getFinerLevel());
      // quick sanity check
      CH_assert(fineML.m_level == lev + 1);
      CoarseAverage& scalAvgDown =
          fineML.m_coarseAverageScalar;

      LevelData<FArrayBox>& fineScal =
          *fineML.m_scalarNew[a_var];
      scalAvgDown.averageToCoarse(levelScal,
                                  fineScal);
    }


    thisMLPtr = thisMLPtr->getCoarserLevel();


  } // end loop over levels

  // clean up temporary scalar storage
  for (int lev = 0; lev <= finest_level; lev++)
  {
    if (scalRefluxRHS[lev] != NULL)
    {
      delete scalRefluxRHS[lev];
      scalRefluxRHS[lev] = NULL;
    }

    if (scalRefluxCorr[lev] != NULL)
    {
      delete scalRefluxCorr[lev];
      scalRefluxCorr[lev] = NULL;
    }


  }


}

Real AMRLevelMushyLayer::doHCreflux()
{
  CH_TIME("AMRLevelMushyLayer::doHCreflux");

  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::doHCreflux on level " << m_level << endl;
  }
  // Two components: enthalpy and salinity
  Interval solverComps(0, 1);
  int numComps = solverComps.size();
  int finest_level = getFinestLevel();

  Vector<Real> AmrDx(finest_level + 1);

  AMRLevelMushyLayer* thisMLPtr = this;
  int startLev = m_level;
  // if crser level exists, define it as well for BC's
  if (startLev > 0)
  {
    startLev = startLev - 1;
    thisMLPtr = thisMLPtr->getCoarserLevel();
  }
  AMRLevelMushyLayer* startLevelPtr = thisMLPtr;


  Vector<LevelData<FArrayBox>*> HCRefluxCorr(finest_level + 1, NULL);
  Vector<LevelData<FArrayBox>*> HCRefluxRHS(finest_level + 1, NULL);

  // Might need these two for nonlinear reflux
  //  Vector<LevelData<FArrayBox>*> HCOld(finest_level + 1, NULL);
  //  Vector<LevelData<FArrayBox>*> HCNew(finest_level + 1, NULL);
  //  Vector<LevelData<FArrayBox>*> fullRHS(finest_level + 1, NULL);

  // startLev is either m_level, or m_level-1 if that is defined
  thisMLPtr = startLevelPtr;
  for (int lev = startLev; lev <= finest_level; lev++)
  {
    AmrDx[lev] = thisMLPtr->m_dx;
    const DisjointBoxLayout& levelGrids = thisMLPtr->m_grids;
    // recall that AMRMultiGrid can only do one component. Bullshit! Can do more now.
    // rhs has no ghost cells. Bullshit!  ghost cell as AMRFASMultigrid creates objects from this and they need a ghost cell
    HCRefluxRHS[lev] = new LevelData<FArrayBox>(levelGrids, numComps, IntVect::Unit);

    //soln has one layer of ghost cells
    IntVect ghostVect(D_DECL(1, 1, 1));
    HCRefluxCorr[lev] = new LevelData<FArrayBox>(levelGrids,  numComps, ghostVect);

    // initialize corr, RHS to 0
    DataIterator levelDit = HCRefluxRHS[lev]->dataIterator();
    LevelData<FArrayBox>& levelRefluxRHS = *(HCRefluxRHS[lev]);
    LevelData<FArrayBox>& levelRefluxCorr =  *(HCRefluxCorr[lev]);
    for (levelDit.reset(); levelDit.ok(); ++levelDit)
    {
      levelRefluxRHS[levelDit()].setVal(0.0);
      levelRefluxCorr[levelDit()].setVal(0.0);
    }
    thisMLPtr = thisMLPtr->getFinerLevel();
  }


  // loop over levels and establish RHS, corr
  thisMLPtr = this;

  for (int lev = m_level; lev <= finest_level; lev++)
  {
    if (lev < finest_level)
    {
      LevelData<FArrayBox>& levelRefluxRHS = *(HCRefluxRHS[lev]);

      DataIterator levelDit = levelRefluxRHS.dataIterator();

      for (levelDit.reset(); levelDit.ok(); ++levelDit)
      {
        levelRefluxRHS[levelDit()].setVal(0.0);
      }
      Real refluxScale = -1.0 / AmrDx[lev]; // petermc made negative, 7 Dec 07

      thisMLPtr->m_fluxRegHC->reflux(levelRefluxRHS, refluxScale);
      // Don't think we should have been doing this setToZero.
      //      thisMLPtr->m_fluxRegHC->setToZero();
    }

    // initial guess for correction is RHS
    HCRefluxRHS[lev]->copyTo(solverComps,
                             *(HCRefluxCorr[lev]), solverComps);

    thisMLPtr = thisMLPtr->getFinerLevel();
  } // end loop over levels to compute RHS

  // Do averaging down on this level and finer
  thisMLPtr = this;
  for (int lev=m_level; lev<finest_level; lev++)
  {
    AMRLevelMushyLayer& fineML =  *(thisMLPtr->getFinerLevel());
    // quick sanity check
    CH_assert(fineML.m_level = lev + 1);

    CoarseAverage& scalAvgDown = fineML.m_coarseAverageHC;

    scalAvgDown.averageToCoarse(*(HCRefluxRHS[lev]), *(HCRefluxRHS[lev+1]));

    thisMLPtr = thisMLPtr->getFinerLevel();
  }

  // define levelOp and solver for diffusive solves;
  // need to define
  // new AMRMultiGrid for each component because they
  // may have different coefficients

  // only do all of this if  actually diffusive
  //  bool forceExplicitReflux = false;
  //  int numLevels = finest_level + 1;

  Vector<AMRLevelMushyLayer*> hierarchy;
  Vector<DisjointBoxLayout> grids;
  Vector<int> refRat;
  ProblemDomain lev0Dom;
  Real lev0Dx;
  getHierarchyAndGrids(hierarchy, grids, refRat, lev0Dom, lev0Dx);

  RelaxSolver<LevelData<FArrayBox> > bottomSolver;
  bottomSolver.m_verbosity = m_opt.verbosity_multigrid;

  Real maxRefluxRHS = computeSum(HCRefluxRHS, refRat, m_dx, Interval(0, numComps-1), m_level);

  int numLevels = finest_level + 1;

  // Need to define a two component VCAMRPoissonOp to solve for enthalpy-salinity corrections

  AMRLevelMushyLayer* finestML = this;
  while(finestML->hasFinerLevel())
  {
    finestML = finestML->getFinerLevel();
  }
  // This is a Helmholtz operator
  Real alpha = 1.0;
  Real beta = m_opt.refluxBetaSign*m_dt;

  int Hcomp = 0;
  int Ccomp = 1;

  Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(finest_level+1);
  Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(finest_level+1);
  Vector<RefCountedPtr<LevelData<FluxBox> > > porosityFace(finest_level+1);

  Vector<RefCountedPtr<LevelData<FArrayBox> > > enthalpySolidus(finest_level+1),
      enthalpyLiquidus(finest_level+1), enthalpyEutectic(finest_level+1), HC(finest_level+1);

  // Not sure if we actually need a ghost vector or not
  IntVect ivGhost = IntVect::Unit;

  AMRLevelMushyLayer* amrML = this;

  // Define these things over ALL levels
  while(amrML->getCoarserLevel())
  {
    amrML = amrML->getCoarserLevel();
  }

  // Almost certain this will always be level 0
  int baseLevel = amrML->level();

  for (int lev=baseLevel; lev<= finest_level; lev++)
  {
    bCoef[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(grids[lev], numComps, ivGhost));
    porosityFace[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(grids[lev], 1, ivGhost));
    aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], numComps, ivGhost));

    enthalpySolidus[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], 1, ivGhost));
    enthalpyLiquidus[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], 1, ivGhost));
    enthalpyEutectic[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], 1, ivGhost));
    HC[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], numComps, ivGhost));

    amrML->fillHC(*HC[lev], m_time);
    amrML->fillScalars(*enthalpySolidus[lev], m_time, m_enthalpySolidus, true, true);
    amrML->fillScalars(*enthalpyLiquidus[lev], m_time, m_enthalpyLiquidus, true, true);
    amrML->fillScalars(*enthalpyEutectic[lev], m_time, m_enthalpyEutectic, true, true);

    amrML->fillScalarFace(*porosityFace[lev], m_time, m_porosity, true, true);

    for (DataIterator dit = bCoef[lev]->dataIterator(); dit.ok(); ++dit)
    {
      (*aCoef[lev])[dit].setVal(1.0);
      (*bCoef[lev])[dit].setVal(1.0);

      // For linear reflux, we don't do this porosity stuff and just solve
      // diffusion coefficient * Lap(correction) = Reflux RHS
      if (m_opt.refluxMethod!=RefluxMethod::LinearReflux)
      {


        // for salt solve
        (*bCoef[lev])[dit].mult((*porosityFace[lev])[dit], (*porosityFace[lev])[dit].box(), 0, Ccomp);

        // bCoef for heat solve
        (*bCoef[lev])[dit].minus((*porosityFace[lev])[dit], (*porosityFace[lev])[dit].box(), 0, Hcomp);
        for (int dir=0; dir<SpaceDim; dir++)
        {
          (*bCoef[lev])[dit][dir].mult(m_parameters.heatConductivityRatio, Hcomp);
        }
        (*bCoef[lev])[dit].plus((*porosityFace[lev])[dit], (*porosityFace[lev])[dit].box(), 0, Hcomp);

      }

      for (int dir=0; dir<SpaceDim; dir++)
      {
        (*bCoef[lev])[dit][dir].mult(-m_scalarDiffusionCoeffs[ScalarVars::m_enthalpy], Hcomp);
        (*bCoef[lev])[dit][dir].mult(-m_scalarDiffusionCoeffs[ScalarVars::m_bulkConcentration], Ccomp);
      }



    }

    amrML = amrML->getFinerLevel();
  } // end loop over levels

  bool homogeneous = true;

  // Try a no flux BC
  BCHolder refluxBC = m_physBCPtr->enthalpySalinityBC(homogeneous);

  AMRMultiGrid<LevelData<FArrayBox> > *diffusionSolver;

  Vector< AMRLevelOp<LevelData<FArrayBox> >*  >  generic_ops;
  Vector< LevelTGAHelmOp<LevelData<FArrayBox>, FluxBox >*  >  ops;

  generic_ops.resize(finest_level+1);
  ops.resize(finest_level+1);

  if (m_opt.refluxMethod== RefluxMethod::LinearReflux || m_opt.refluxMethod == RefluxMethod::LinearVCReflux)
  {

    VCAMRPoissonOp2Factory diffusiveOpFactory;
    diffusiveOpFactory.define(lev0Dom, grids,
                              refRat, AmrDx[0],
                              refluxBC,
                              alpha, aCoef,
                              beta, bCoef);

    AMRLevelOpFactory<LevelData<FArrayBox> >& castFact =
        (AMRLevelOpFactory<LevelData<FArrayBox> >&) diffusiveOpFactory;

    diffusionSolver = new AMRMultiGrid<LevelData<FArrayBox> >;
    diffusionSolver->define(lev0Dom, castFact,
                            &bottomSolver, numLevels);

  }
  else if (m_opt.refluxMethod == RefluxMethod::NonlinearReflux)
  {
    // nonlinear reflux

    //      MushyLayerParams* mlParamsPtr = &m_parameters;
    BCHolder temperature_Sl_BC = m_physBCPtr->temperatureLiquidSalinityBC(homogeneous);
    //      BCHolder HC_BC  = m_physBCPtr->enthalpySalinityBC(homogeneous);
    EdgeVelBCHolder porosityEdgeBC(m_physBCPtr->porosityFaceBC());



    // Calculate modified diffusion coefficient
    /*
     * For heat equation:
     * alpha = 1 (liquid); 0.25 (mush); 0 (eutectic); 1/cp (solid)
     *
     * For salt equation:
     * alpha = 1 (liquid); 1 (mush); 0 (eutectic); 0 (solid)
     */

    for (int lev=m_level; lev<= finest_level; lev++)
    {

      // First need bounding energies
      for (DataIterator dit = HC[lev]->dataIterator(); dit.ok(); ++dit)
      {

        FArrayBox& Hs = (*enthalpySolidus[lev])[dit];
        FArrayBox& He = (*enthalpyEutectic[lev])[dit];
        FArrayBox& Hl = (*enthalpyLiquidus[lev])[dit];

        FArrayBox& HCfab = (*HC[lev])[dit];

        Box region = HCfab.box();
        region &= Hs.box();

        FORT_CALCULATE_BOUNDING_ENERGY( CHF_CONST_FRA(HCfab),
                                        CHF_FRA(Hs),
                                        CHF_FRA(He),
                                        CHF_FRA(Hl),
                                        CHF_BOX(region),
                                        CHF_CONST_REAL(m_parameters.compositionRatio),
                                        CHF_CONST_REAL(m_parameters.waterDistributionCoeff),
                                        CHF_CONST_REAL(m_parameters.specificHeatRatio),
                                        CHF_CONST_REAL(m_parameters.stefan),
                                        CHF_CONST_REAL(m_parameters.thetaEutectic),
                                        CHF_CONST_REAL(m_parameters.ThetaEutectic));


        // Update bCoef
        FluxBox& thisBCoef = (*bCoef[lev])[dit];

#if CH_SPACEDIM == 1
        FORT_MODIFY_DIFFUSION_COEFFS1D
#elif CH_SPACEDIM == 2
        FORT_MODIFY_DIFFUSION_COEFFS2D
#elif CH_SPACEDIM == 3
        FORT_MODIFY_DIFFUSION_COEFFS3D
#endif
        (CHF_CONST_FRA(HCfab),
         CHF_CONST_FRA(Hs),
         CHF_CONST_FRA(He),
         CHF_CONST_FRA(Hl),
#if CH_SPACEDIM >= 1
         CHF_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
         CHF_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
         CHF_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
         This_will_not_compile!
#endif
         CHF_BOX(region),
         CHF_CONST_REAL(m_parameters.specificHeatRatio)
        );

      } // end loop over grids

    } // end loop over levels

    VCAMRPoissonOp2Factory diffusiveOpFactory;
    diffusiveOpFactory.define(lev0Dom, grids,
                              refRat, AmrDx[0],
                              refluxBC,
                              alpha, aCoef,
                              beta, bCoef);

    AMRLevelOpFactory<LevelData<FArrayBox> >& castFact =
        (AMRLevelOpFactory<LevelData<FArrayBox> >&) diffusiveOpFactory;

    //      diffusionSolver = new AMRFASMultiGrid<LevelData<FArrayBox> >;
    //      diffusionSolver->define(lev0Dom, castFact,
    //                              &bottomSolver, numLevels);
    diffusionSolver = new AMRMultiGrid<LevelData<FArrayBox> >;
    diffusionSolver->define(lev0Dom, castFact,
                            &bottomSolver, numLevels);


  } // end loop over different reflux types
  else
  {
    MayDay::Error("Unknown reflux method");
  }

  generic_ops = diffusionSolver->getAMROperators();

  for (int ilev = 0; ilev < generic_ops.size(); ilev++)
  {
    ops[ilev] = dynamic_cast<LevelTGAHelmOp<LevelData<FArrayBox>,FluxBox>* >(generic_ops[ilev]);
    if (ops[ilev]==NULL)
    {
      MayDay::Error("dynamic cast failed---is that operator really a TGAHelmOp?");
    }
  }

  diffusionSolver->m_verbosity = m_opt.AMRMultigridVerb;
  diffusionSolver->m_eps = m_opt.AMRMultigridTolerance;
  diffusionSolver->m_normThresh = m_opt.AMRMultigridNormThresh;
  diffusionSolver->m_hang = m_opt.AMRMultigridHang;

  //    Interval solverComps(0,numComps-1);

  // now solve
  diffusionSolver->solve(HCRefluxCorr, HCRefluxRHS,
                         finest_level, m_level, false, // don't initialize to zero
                         true);  // force homogeneous


  // now increment scalars with reflux correction
  // go from finest->coarsest so that we can also avgDown
  // increment NS pointer to finest level
  thisMLPtr = this;
  while (!thisMLPtr->finestLevel())
  {
    thisMLPtr = thisMLPtr->getFinerLevel();
  }
  CH_assert(thisMLPtr->m_level == finest_level);


  // Check if we're doing averaging down
//  bool m_opt.doAvgDown = true;
//  pp.query("reflux_average_down", m_opt.doAvgDown);

  // Add correction to H-C fields
  // Also recalculate fluxes
  for (int lev = finest_level; lev >= m_level; lev--)
  {
    LevelData<FArrayBox>& levelCorr = *(HCRefluxCorr[lev]);
    DataIterator levelDit = levelCorr.dataIterator();
    for (levelDit.reset(); levelDit.ok(); ++levelDit)
    {

      //(*thisMLPtr->m_scalarNew[ScalarVars::m_enthalpy])[levelDit()].plus(levelCorr[levelDit()], 0, 0, 1);
      //(*thisMLPtr->m_scalarNew[ScalarVars::m_bulkConcentration])[levelDit()].plus(levelCorr[levelDit()], 1, 0, 1);

      if (m_opt.reflux_enthalpy)
      {
        (*thisMLPtr->m_scalarNew[ScalarVars::m_enthalpy])[levelDit()].plus(levelCorr[levelDit()],m_opt.refluxCorrSign, 0, 0, 1);
      }

      if (m_opt.reflux_concentration)
      {
        (*thisMLPtr->m_scalarNew[ScalarVars::m_bulkConcentration])[levelDit()].plus(levelCorr[levelDit()], m_opt.refluxCorrSign, 1, 0, 1);
      }
    }

    // Compute other fields on the finest level
    // for coarser levels, defer this until we've done averaging down
    if (lev == finest_level || !m_opt.refluxAverageDown)
    {
      thisMLPtr->updateEnthalpyVariables();
    }

    thisMLPtr = thisMLPtr->getCoarserLevel();

  } // end loop over levels

  // Clean up
  if (diffusionSolver != NULL)
  {
    delete diffusionSolver;
    diffusionSolver = NULL;
  }

  // clean up temporary scalar storage
  for (int lev = 0; lev <= finest_level; lev++)
  {
    if (HCRefluxRHS[lev] != NULL)
    {
      delete HCRefluxRHS[lev];
      HCRefluxRHS[lev] = NULL;
    }

    if (HCRefluxCorr[lev] != NULL)
    {
      delete HCRefluxCorr[lev];
      HCRefluxCorr[lev] = NULL;
    }

  }

  return maxRefluxRHS;

}

