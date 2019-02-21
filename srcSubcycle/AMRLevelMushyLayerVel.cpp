#include "AMRLevelMushyLayer.H"

/**
 * This source file contains methods relating to the velocity field
 */


void AMRLevelMushyLayer::calculateTimeIndAdvectionVel(Real time, LevelData<FluxBox>& a_advVel)
{

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::calculateTimeIndAdvectionVel (level " << m_level << ")"    << endl;
  }
  CH_TIME("AMRLevelMushyLayer::calculateTimeIndAdvectionVel");

  IntVect ivGhost = m_numGhost*IntVect::Unit;
  IntVect advectionGhost = m_numGhostAdvection*IntVect::Unit;

  LevelData<FArrayBox> src(m_grids, SpaceDim, ivGhost);
  LevelData<FArrayBox> vel(m_grids, SpaceDim, advectionGhost);

  DataIterator dit = m_grids.dataIterator();

  // The pressure scale is permeability for darcy flow,
  // and porosity for darcy-brinkman

  RefCountedPtr<LevelData<FluxBox> > crsePressureScaleEdgePtr, pressureScaleEdgePtr;
  RefCountedPtr<LevelData<FArrayBox> > crsePressureScalePtr, pressureScalePtr;

  RefCountedPtr<LevelData<FArrayBox> > crsePressurePtr, pressurePtr;

  AMRLevelMushyLayer* amrMLcrse = getCoarserLevel();

  LevelData<FluxBox>* velocityBCVals = new LevelData<FluxBox>(m_grids, 1, a_advVel.ghostVect());
  LevelData<FluxBox> gradP(m_grids, 1, IntVect::Zero);
  LevelData<FluxBox> Theta_l_face(m_grids, 1, IntVect::Unit);

  m_projection.gradPhi(gradP);
  fillScalarFace(Theta_l_face, time, m_liquidConcentration, true, false);

  // Construct a fluxbox containing values we want to enforce at the boundary (don't usually end up using this)
  for (DataIterator dit = velocityBCVals->dataIterator(); dit.ok(); ++dit)
  {
    FluxBox& bcVel = (*velocityBCVals)[dit];

    FArrayBox& bc_u = bcVel[0];
    FArrayBox& bc_v = bcVel[1];

    FArrayBox& Sl = Theta_l_face[dit][1];
    FArrayBox& gradP_z = gradP[dit][1];

    bc_u.setVal(0.0);

    // vertical velocity boundary cells
    Box v_box = bc_v.box();
    Box valid = m_grids[dit];
    valid &= m_problem_domain.domainBox();

    Box toRegion(valid);

    // Setting BC for vertical velocity (component = 1) in vertical direction (idir=1) on low side.
    int velComp = 1;
    Side::LoHiSide side = Side::Lo;
    int idir = 1;

    if (idir == velComp)
    {
      toRegion.surroundingNodes(idir);
      int coord = toRegion.sideEnd(side)[idir];
      toRegion.setRange(idir, coord);

    }
    else
    {
      toRegion.surroundingNodes(velComp);
      toRegion = adjCellBox(toRegion, idir, side, 1);

    }

    toRegion &= v_box;

    for (BoxIterator bit = BoxIterator(toRegion); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      bc_v(iv) = -gradP_z(iv) - m_parameters.m_buoyancySCoeff*Sl(iv);

    }

  }

  EdgeVelBCHolder edgeVelBC(m_physBCPtr->edgeVelFuncBC(m_opt.viscousBCs, velocityBCVals));

  VelBCHolder velBC(m_physBCPtr->uStarFuncBC(m_opt.viscousBCs));
  VelBCHolder velBCExtrap(m_physBCPtr->velExtrapBC());

  calculatePermeability();

  // if a coarser level exists, will need coarse-level data for proj
  if (m_level > 0)
  {
    const DisjointBoxLayout& crseGrids = amrMLcrse->m_grids;
    crsePressurePtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(crseGrids, 1));
    // coarse velocity BC data is interpolated in time

    amrMLcrse->fillScalars(*crsePressurePtr, time, ScalarVars::m_pressure, true);

    if (m_opt.scaleP_MAC)
    {
      crsePressureScaleEdgePtr = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(amrMLcrse->m_grids, 1, IntVect::Unit));
      amrMLcrse->fillScalarFace(*crsePressureScaleEdgePtr, time, m_pressureScaleVar, true);

      crsePressureScalePtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(amrMLcrse->m_grids, 1, IntVect::Unit));
      amrMLcrse->fillScalars(*crsePressureScalePtr, time, m_pressureScaleVar, true);
    }

    if (s_verbosity >= 5)
    {
      pout() << "AMRLevelMushyLayer::calculateTimeIndAdvectionVel - got CF BCs (level " << m_level << ")"    << endl;
    }

  }

  if(m_opt.scaleP_MAC)
  {
    pressureScaleEdgePtr = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids, 1, IntVect::Unit));
    fillScalarFace(*pressureScaleEdgePtr, time, m_pressureScaleVar, true);

    pressureScalePtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grids,1,IntVect::Unit));
    fillScalars(*pressureScalePtr, time, m_pressureScaleVar, true);
  }

  IntVect ghost = IntVect::Unit;

  LevelData<FluxBox> gradPhi(m_grids, 1);

  LevelData<FArrayBox> T(m_grids, 1, IntVect::Unit);
  LevelData<FArrayBox> C(m_grids, 1, IntVect::Unit);
  LevelData<FArrayBox> porosity(m_grids, 1, IntVect::Unit);

  fillScalars(porosity, time, m_porosity,            true);
  fillScalars(T,        time, m_temperature,         true);
  fillScalars(C,        time, m_liquidConcentration, true);

  m_projection.gradPhi(gradPhi);
  bool alreadyHasPressure = false;

  if (m_parameters.isViscous())
  {
    if (s_verbosity >= 5)
    {
      pout() << "AMRLevelMushyLayer::calculateTimeIndAdvectionVel - with viscosity (level " << m_level << ")"    << endl;
    }

    // Average U^* to faces for projection
    CellToEdge(*m_vectorNew[VectorVars::m_Ustar], a_advVel);
  }
  else
  {

    if (s_verbosity >= 5)
    {
      pout() << "AMRLevelMushyLayer::calculateTimeIndAdvectionVel - no viscosity (level " << m_level << ")"    << endl;
    }

    // Darcy velocity
    // Just doing Darcy. U^* = permeability * (RaT * theta - RaC*Theta)

    fillUnprojectedDarcyVelocity(a_advVel, time);

    // Subtract best guess at pressure by default for non-AMR sims
    bool useIncrementalPressure = (getMaxLevel() == 0) && !m_opt.useIncrementalPressure;

    // Don't do this on refined levels as it messes up the CF boundary condition
    if (useIncrementalPressure && m_level == 0)
    {
      // Before doing this, check that the pressure isn't getting too big
      // Only use incremental pressure if max (pressure) < 10^3
      // hopefully this is stable
      // if we don't do this, the pressure can grow out of control

      Real maxPressure = ::computeNorm(m_projection.phi(), NULL, 1, m_dx, Interval(0,0), 0);
      // This is chosen empirically and may need some refinement
      Real pressureCap = 1e3;

      // A better indicator of issues arising is that the pressure becomes negative
      // as we fix P = 0 on the bottom boundary, and it should increase up through the domain, use this condition instead
      Real minPressure = ::computeMin(m_projection.phi(), NULL, 1, Interval(0,0));
      if (minPressure > 0 && maxPressure < pressureCap)
      {

        Real phiScale = m_projection.getScale(m_opt.phiScale, m_dt);

        // Apply last calculated MAC correction
        m_projection.setPressureScaleEdgePtr(pressureScaleEdgePtr);

        m_projection.applyMacCorrection(a_advVel,
                                        NULL,
                                        phiScale);
        alreadyHasPressure = true;
      }

    }

  }

  // Apply domain BCs and project
  edgeVelBC.applyBCs(a_advVel, m_grids,
                     m_problem_domain, m_dx,
                     false); // inhomogeneous

  // Do some smoothing
  int nghost = a_advVel.ghostVect()[0];

  for (DataIterator dit = a_advVel.dataIterator(); dit.ok(); ++dit)
  {
    for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& velDir = a_advVel[dit][dir];
      Box b = velDir.box();
      b.grow(-(nghost+2));

      for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        Real lapU = -4*velDir(iv) +
            velDir(iv+BASISV(0)) + velDir(iv-BASISV(0)) + velDir(iv+BASISV(1)) + velDir(iv-BASISV(1));

        velDir(iv) = velDir(iv) - m_opt.smoothingCoeff*lapU;
      }
    }

  }

  // Allow multiple projections
  // This doesn't actually really help
  Real maxDivU = 1e10;
  int proj_i = 1;
  int maxNumProj = 1;

  Real correctScale = 1.0;

  while(maxDivU > 1e-10 && proj_i <= maxNumProj)
  {

    if (s_verbosity >= 5)
    {
      pout() << "AMRLevelMushyLayer::calculateTimeIndAdvectionVel - do projection (level " << m_level << ")"    << endl;
    }

    // Copy to U^* before projection
    EdgeToCell(a_advVel, *m_vectorNew[VectorVars::m_advUstar]);

    a_advVel.exchange();
    int exitStatus = m_projection.levelMacProject(a_advVel, m_dt, crsePressurePtr, pressureScalePtr,
                                 crsePressureScalePtr, pressureScaleEdgePtr, crsePressureScaleEdgePtr,
                                 alreadyHasPressure, correctScale);

    correctScale = correctScale/10;
//    correctScale = 0.0;

    a_advVel.exchange();
    Divergence::levelDivergenceMAC(*m_scalarNew[ScalarVars::m_divUadv], a_advVel, m_dx);
    maxDivU = ::computeNorm(*m_scalarNew[ScalarVars::m_divUadv], NULL, 1, m_dx, Interval(0,0), 0);
    Real minPressure = ::computeMin(m_projection.phi(), NULL, 1, Interval(0,0));
    pout() << "  MAC Projection (level "<< m_level << "), exit status = " << exitStatus
        << ", max(div u) = " << maxDivU << ", min pressure = " << minPressure << endl;



    proj_i++;
  }


  fillAdvVel(time, a_advVel);
  a_advVel.exchange();

  m_projection.getPhi(*m_scalarNew[ScalarVars::m_pressure]);

  // Test - overwrite with analytic advection velocity
//  bool m_opt.enforceAnalyticVel = (m_parameters.physicalProblem == PhysicalProblems::m_soluteFluxTest  );

  if (m_opt.enforceAnalyticVel)
  {
    fillAnalyticVel(a_advVel);
  } // end if enforce analytic vel

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::calculateTimeIndAdvectionVel - freestream correction (level " << m_level << ")"    << endl;
  }

  // Finally - apply freestream preservation correction from previous synchronisation step
  m_projection.applyFreestreamCorrection(a_advVel);

  // Face centred velocity field should be divergence free (to machine precision)
  Divergence::levelDivergenceMAC(*m_scalarNew[ScalarVars::m_divUadv], a_advVel, m_dx);

  // Need a m_fluidVel to be populated to some advection things the right way around
  EdgeToCell(m_advVel, *m_vectorNew[VectorVars::m_fluidVel]);
  this->fillVectorField(*m_vectorNew[VectorVars::m_fluidVel], time, VectorVars::m_fluidVel);

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::calculateTimeIndAdvectionVel - finished (level " << m_level << ")"    << endl;
  }

  // CLean up
  if (velocityBCVals != NULL)
  {
    delete velocityBCVals;
    velocityBCVals = NULL;

  }

}

void
AMRLevelMushyLayer::computeAdvectionVelocities(LevelData<FArrayBox>& advectionSourceTerm, Real advVelCentering)
{
  CH_TIME("AMRLevelMushyLayer::computeAdvectionVelocities()");

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::computeAdvectionVelocities (level " << m_level << ")" << endl;
  }

  DataIterator dit(m_grids);
  IntVect ivGhost = m_numGhostAdvection*IntVect::Unit;
  Real old_time = m_time - m_dt;

  // vel time is the time at which ths advection velocity is calculated
  // default is time centered, so at old_time + 0.5*dt
  Real vel_time = old_time + advVelCentering*m_dt;
  Real half_dt = vel_time - old_time;

  LevelData<FArrayBox> porosityAv(m_grids, 1, ivGhost);
  RefCountedPtr<LevelData<FluxBox> > porosityFaceAvPtr = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids, 1, ivGhost));

  EdgeVelBCHolder edgeVelBC(m_physBCPtr->advectionVelFuncBC(m_opt.viscousBCs));

  LevelData<FArrayBox> velOldGrown(m_grids, SpaceDim, ivGhost);
  fillVectorField(velOldGrown, old_time, m_fluidVel, true);

  // Get initial guess at advection velocity from U^n

  if (m_opt.useOldAdvVel)
  {
    // do nothing - we already have m_advVel from previous timestep
  }
  else
  {
    CellToEdge(velOldGrown, m_advVel);
  }

  edgeVelBC.applyBCs(m_advVel, m_grids, m_problem_domain, m_dx,
                     false); // inhomogeneous

  if (m_opt.doEulerPart)
  {
    if (s_verbosity >= 5)
    {
      pout() << "AMRLevelMushyLayer::computeAdvectionVelocities - doing advection of velocities" << endl;
    }

    if (m_opt.implicitAdvectionSolve)
    {

      bool doFRupdates = false;
      bool doProjection = false; // don't do CC projection
      bool compute_uDelU = true;

      // do full timestep - this should now give the same velocity as the later solve?
      advectionSourceTerm.exchange();
      computeCCvelocity(advectionSourceTerm, old_time, m_dt, doFRupdates, doProjection,
                        compute_uDelU, m_opt.usePhiForImplicitAdvectionSolve);

      m_vectorNew[VectorVars::m_viscousSolveSrc]->copyTo(*m_vectorNew[VectorVars::m_advectionImplicitSrc]);

      // New version:
      if (m_opt.usePhiForImplicitAdvectionSolve)
      {
        fillVectorField(*m_vectorNew[VectorVars::m_advUpreProjection], vel_time, m_advUpreProjection, false, true);

      }
      else
      {
        for (DataIterator dit = m_vectorNew[VectorVars::m_advUpreProjection]->dataIterator(); dit.ok(); ++dit)
        {
          (*m_vectorNew[VectorVars::m_advUpreProjection])[dit].copy((*m_vectorNew[VectorVars::m_UpreProjection])[dit]);
        }
        fillVectorField(*m_vectorNew[VectorVars::m_advUpreProjection], vel_time, m_UpreProjection, false, true);

      }


      CellToEdge(*m_vectorNew[VectorVars::m_advUpreProjection], m_advVel);


    }
    else
    {
      if (s_verbosity >= 5)
      {
        pout() << "AMRLevelMushyLayer::computeAdvectionVelocities() - explicit tracing scheme" << endl;
      }

      // Explicit advection solve (tracing)
//      int saveAdvectionMethod = m_advectionMethod;

      // Option to do a different method here
//      int velPredMethod = m_advectionMethod;
//      ppMain.query("velPredictionMethod", velPredMethod);
//      m_advectionMethod = velPredMethod; //m_noPorosity;

      fillScalarFace(*porosityFaceAvPtr, vel_time, m_porosity, true);
      fillScalars(porosityAv, vel_time, m_porosity, true);

      // Construct the thing to advect
      LevelData<FArrayBox> porosityGrown(m_grids, 1, ivGhost);
      LevelData<FArrayBox> U_to_advect(m_grids, SpaceDim, ivGhost);
      fillScalars(porosityGrown, old_time, m_porosity, true);

      LevelData<FArrayBox> ccAdvVel(m_grids, SpaceDim, ivGhost);
      fillVectorField(ccAdvVel, old_time, m_fluidVel, true);

      setVelZero(ccAdvVel, m_opt.advPorosityLimit);

      /*
       * Options:
       * m_porosityInAdvection:
       *    solve d(u/porosity)/dt + (u/porosity).grad(u/porosity) = S/porosity (-(u/porosity^2)dporosity/dt)
       *    (assume d(porosity)/dt = 0)
       * m_porosityOutsideAdvection:
       *    solve du/dt + (u/porosity).grad(u) = S + (u.grad(porosity)/porosity^2)u
       * m_noPorosity:
       *    solve du/dt + u.grad(u) = S
       */

      for (DataIterator dit = U_to_advect.dataIterator(); dit.ok(); ++dit)
      {
        U_to_advect[dit].setVal(0.0);
        U_to_advect[dit] += velOldGrown[dit];
        setVelZero(U_to_advect, m_opt.advPorosityLimit);

        // Divide each component of velocity by porosity if we're advecting u/chi
        if ( m_opt.advectionMethod == m_porosityInAdvection)
        {
          // We're advecting u/chi
          for (int dir = 0; dir < SpaceDim; dir++)
          {
            U_to_advect[dit].divide( porosityGrown[dit], porosityGrown[dit].box(), 0, dir, 1);

            // Also need to scale source term
            // this has already been done in computeAdvectionVelSourceTerm
//            advectionSourceTerm[dit].divide( porosityGrown[dit], porosityGrown[dit].box(), 0, dir, 1);
          }

        }

        // Set advection velocity = u/chi
        // we now do this for porosity in advection operator too
        if  (m_opt.advectionMethod == m_porosityOutsideAdvection
            || m_opt.advectionMethod == m_porosityInAdvection )
        {

          setVelZero(m_advVel, m_opt.advPorosityLimit);

          m_advVel[dit].divide((*porosityFaceAvPtr)[dit], m_advVel[dit].box(), 0, 0);
          for (int dir = 0; dir < SpaceDim; dir++)
          {
            ccAdvVel[dit].divide(porosityAv[dit], ccAdvVel[dit].box(), 0, dir);
          }

          // this does BCs (doesn't fill interior)
          fillVectorField(ccAdvVel, old_time, m_fluidVel, false);
        }

        if ( m_opt.advectionMethod == m_noPorosity)
        {
          // Don't need to do anything here
        }

        U_to_advect[dit].mult(m_parameters.m_advectionCoeff);

      } // end loop over grids

      // Set advection velocity to zero where porosity is low
      //        setVelZero(U_to_advect, 0.05);

      Real maxAdvectionVel = getMaxAdvVel();
      Real cfl = maxAdvectionVel*m_dt/m_dx;
      if (cfl > 1.0)
      {
        pout() << "WARNING - about to do trace adv vel with cfl = " << cfl << endl;
      }

      //        pout() << "Do trace advection vel with cfl = " << cfl << endl;

      edgeVelBC.applyBCs(m_advVel, m_grids, m_problem_domain, m_dx,
                                   false); // inhomogeneous

      // m_dt is the full timestep, and this always returns the half time velocity
      // changed this from m_dt (full timestep) to half_dt*2 incase we want to compute velocities
      // at some odd point in time
      traceAdvectionVel(m_advVel, ccAdvVel, U_to_advect,
                        advectionSourceTerm, m_patchGodVelocity, old_time, half_dt*2);

//      m_advVel.exchange();

      EdgeToCell(m_advVel, *m_vectorNew[VectorVars::m_advectionVel]);
      Real maxAdvU = ::computeNorm(*m_vectorNew[VectorVars::m_advectionVel], NULL, 1, m_dx, Interval(0,0));

      if (s_verbosity >= 1 && maxAdvU > 1e10)
      {
        pout() << "WARNING - max advection velocity (pre-projection) = " << maxAdvU << endl;
      }


      // Need to recover the advection velocity, as
      // traceAdvectionVel will replce m_advVel with the upwinded U_to_advect
      if ( m_opt.advectionMethod == m_porosityInAdvection)
      {
        for (dit.reset(); dit.ok(); ++dit)
        {
          m_advVel[dit].mult((*porosityFaceAvPtr)[dit], m_advVel[dit].box(), 0, 0);
        }
      }


//      m_advectionMethod = saveAdvectionMethod;

      // DO some post tracing smoothing in the x-direction
      // this is necessary to remove some \delta x scale instabilities which seem to arise

      m_advVel.exchange();

      horizontallySmooth(m_advVel);

      m_advVel.exchange();


    } // end if implicit/explicit solve

    setVelZero(m_advVel, m_opt.advPorosityLimit);

    //m_advVel contains predicted face centred velocity at half time step
    //Now need to project it and do lambda corrections etc
    Real project_dt = 2*half_dt; // this usually just works out at m_dt, for time centred advection velocities
    if (m_opt.implicitAdvectionSolve)
    {
      project_dt = m_dt;
    }

    correctEdgeCentredVelocity(m_advVel, project_dt);

    // Maybe replace m_advVel with time independent version in low porosity regions?

    if (m_opt.chiLimit > 0.0)
    {
      LevelData<FluxBox> mushyAdvVel(m_advVel.disjointBoxLayout(), 1, m_advVel.ghostVect());
      //        fillUnprojectedDarcyVelocity(mushyAdvVel, m_time-m_dt);
      calculateTimeIndAdvectionVel(m_time-m_dt, mushyAdvVel);

      int count = 0;
      for (DataIterator dit = m_advVel.dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& chi = (*m_scalarNew[ScalarVars::m_porosity])[dit];
        for (int dir=0; dir < SpaceDim; dir++)
        {
          Box b = m_advVel[dit][dir].box();
          FArrayBox& velDir = m_advVel[dit][dir];
          FArrayBox& newVelDir = mushyAdvVel[dit][dir];

          b = b.enclosedCells();

          b &= chi.box();

          for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
          {
            IntVect iv = bit();
            if (chi(iv) < m_opt.chiLimit)
            {
              velDir(iv) = newVelDir(iv);
              count++;
            }
          }
        }
      }

      pout() << "  ComputeAdvectionVelocities - replaced " << count << " cells with darcy flow" << endl;

      // Project again
      correctEdgeCentredVelocity(m_advVel, project_dt);
    }

    // The projection solve isn't always quite good enough, which can leave U~1e-10 where porosity ~1e-15
    // want to avoid u/chi > 1 for porosity =1e-15, so do it manually
    //    setVelZero(m_advVel, advPorosityLimit);

    LevelData<FluxBox> porosityFace(m_advVel.disjointBoxLayout(), 1, m_advVel.ghostVect());
    fillScalarFace(porosityFace, m_time-m_dt, m_porosity, true, true);
    porosityFace.exchange();

  }
  else
  {
    if (s_verbosity >= 5)
    {
      pout() << "AMRLevelMushyLayer::computeAdvectionVelocities - computing time independent velocities" << endl;
    }

    // if we're not doing advection of velocities, calculate the velocities for scalar advection by averaging our CC velocities to faces at the half time
    // and projecting them to ensure they're divergence free
    LevelData<FArrayBox> vel(m_grids, SpaceDim, ivGhost);
    Real time = m_time-m_dt;

    fillVectorField(vel, time, m_UpreProjection, true);

    CellToEdge(vel, m_advVel);

    //    fillScalarFace(m_advVel, time, m_UpreProjection, true);

    m_advVel.exchange();

    EdgeToCell(m_advVel, *m_vectorNew[VectorVars::m_advUstar]);
    correctEdgeCentredVelocity(m_advVel, m_dt);
  }

  if (m_opt.enforceAnalyticVel)
  {
    fillAnalyticVel(m_advVel);

    LevelData<FluxBox> velDiff(m_advVel.disjointBoxLayout(), 1, m_advVel.ghostVect());

    for(DataIterator dit = velDiff.dataIterator(); dit.ok(); ++dit)
    {
      velDiff[dit].copy(m_advVel[dit]);
    }


    if (m_opt.projectAnalyticVel)
    {
      correctEdgeCentredVelocity(m_advVel, m_dt);
    }


    for(DataIterator dit = velDiff.dataIterator(); dit.ok(); ++dit)
    {
      velDiff[dit].minus(m_advVel[dit], m_advVel[dit].box(), 0,0,1);
    }


  } // end if enforce analytic vel

  Divergence::levelDivergenceMAC(*m_scalarNew[ScalarVars::m_divUadv], m_advVel, m_dx);

  // Finally
  EdgeToCell(m_advVel, *m_vectorNew[VectorVars::m_advectionVel]);

}




void AMRLevelMushyLayer::computeCCvelocity(const LevelData<FArrayBox>& advectionSourceTerm,
                                           Real a_oldTime,
                                           Real a_dt,
                                           bool doFRupdates,
                                           bool a_doProjection,
                                           bool compute_uDelU,
                                           bool a_MACprojection)
{
  CH_TIME("AMRLevelMushyLayer::computeCCvelocity");

  int ustarVar = m_Ustar;
  int uvar = m_fluidVel;
  int uPreprojectionVar = m_UpreProjection;

  if (a_MACprojection)
  {
    ustarVar = m_advUstar;
    uvar = m_advectionVel;
    uPreprojectionVar = m_advUpreProjection;
  }

  VelBCHolder velBC(m_physBCPtr->uStarFuncBC(m_opt.viscousBCs));
  RefCountedPtr<LevelData<FluxBox> > crsePressureScaleEdgePtr, pressureScaleEdgePtr;
  RefCountedPtr<LevelData<FArrayBox> > crsePressureScalePtr, pressureScalePtr;
  LevelData<FArrayBox> *crseVelPtr = NULL;

  AMRLevelMushyLayer* amrMLcrse = getCoarserLevel();

  DataIterator dit = m_scalarNew[ScalarVars::m_porosity]->dataIterator();

  Real half_time = a_oldTime + a_dt/2;
  Real new_time = a_oldTime + a_dt;

  // Calculate Ustar
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::advance - computeUstar (level " << m_level << ")" << endl;
  }

  computeUstar(*m_vectorNew[VectorVars::m_UdelU], advectionSourceTerm, a_oldTime, a_dt, doFRupdates, a_MACprojection, compute_uDelU);

  Real maxUstar = ::computeNorm(*m_vectorNew[ustarVar], NULL, 1, m_dx, Interval(0,0), 0);
  if (s_verbosity >= 1 && maxUstar > 1e10)
  {
    pout() << "WARNING - max cell-centred velocity (pre-projection) = " << maxUstar << endl;
  }

  // Don't think this really helps
  Real porositySmoothing = 0.001;
  porositySmoothing = 0.0;

  if (m_opt.scaleP_CC)
  {
    // Usually geometric averaging, but try arithmetic
    pressureScaleEdgePtr= RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids, 1));

    // Need a ghost cell here for 5 point stencil
    pressureScalePtr= RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grids, 1, IntVect::Unit));

    fillScalarFace(*pressureScaleEdgePtr, half_time, m_pressureScaleVar, arithmeticAveraging, true, false, porositySmoothing);
    fillScalars(*pressureScalePtr, half_time, m_pressureScaleVar, true);

  }

  // if a coarser level exists, will need coarse-level data for proj
  if (m_level > 0)
  {
    const DisjointBoxLayout& crseGrids = amrMLcrse->m_grids;
    crseVelPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim);
    // coarse velocity BC data is interpolated in time

    if (m_addSubtractGradP)
    {
      // Coarse BC is the full projected velocity
      // we subtract off grad(P) inside the level project function

      amrMLcrse->fillVectorField(*crseVelPtr, new_time, uvar, true);
    }
    else
    {
      // If we're not doing +/- grad(P) stuff then just pass in the velocity before it was projected
      // else the divergence at coarse-fine boundaries will be crazy

      amrMLcrse->fillVectorField(*crseVelPtr, new_time, ustarVar, true);
    }

    if (m_opt.scaleP_CC)
    {
      crsePressureScaleEdgePtr = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(amrMLcrse->m_grids, 1));
      amrMLcrse->fillScalarFace(*crsePressureScaleEdgePtr, half_time, m_pressureScaleVar, true, porositySmoothing);

      // Need a ghost cell for 5 point stencil
      crsePressureScalePtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(amrMLcrse->m_grids, 1, IntVect::Unit));
      amrMLcrse->fillScalars(*crsePressureScalePtr, half_time, m_pressureScaleVar, true);
    }
  }

  //Put Ustar into U for projection
  for (dit.reset(); dit.ok(); ++dit)
  {
    (*m_vectorNew[uvar])[dit].copy((*m_vectorNew[ustarVar])[dit]);
  }

  // need to do physical boundary conditions and exchanges
  velBC.applyBCs(*m_vectorNew[uvar], m_grids, m_problem_domain,
                 m_dx, false); // inhomogeneous
  m_vectorNew[uvar]->exchange();

//  if (m_doProjection && a_doProjection)
//    {

  // Before we project, remove grad(P)]
    if (m_addSubtractGradP)
    {
      LevelData<FArrayBox> gradP(m_grids, SpaceDim);
      LevelData<FArrayBox> pressureScale(m_grids, 1);
      fillScalars(pressureScale, half_time, m_pressureScaleVar, true);
      Real old_time = a_oldTime;

      // first get gradPi
      fillPressureSrcTerm(gradP, pressureScale, old_time-a_dt/2, a_MACprojection);


      for (dit.reset(); dit.ok(); ++dit)
      {
        gradP[dit] *= a_dt ;
        (*m_vectorNew[uvar])[dit] += gradP[dit];
      }
    }

  // need to do physical boundary conditions and exchanges
  velBC.applyBCs(*m_vectorNew[uvar], m_grids, m_problem_domain,
                 m_dx, false); // inhomogeneous
  m_vectorNew[uvar]->exchange();

  for (dit.reset(); dit.ok(); ++dit)
  {
    (*m_vectorNew[uPreprojectionVar])[dit].copy((*m_vectorNew[uvar])[dit]);
  }


  if (m_opt.doProjection && a_doProjection)
    {


    CH_TIME("AMRLevelMushyLayer::levelProject");

    QuadCFInterp interp;
    Divergence::levelDivergenceCC(*m_scalarNew[ScalarVars::m_divU], *m_vectorNew[uvar], NULL, m_dx, true, interp);

    Real initMaxDivU = computeNorm(*m_scalarNew[ScalarVars::m_divU], NULL, -1, m_dx, Interval(0,0), 0);
    if (s_verbosity >= 3)
    {
      pout() << "  CCProjection: init max(div(U)) = " << initMaxDivU << endl;
    }

    Real maxDivU = initMaxDivU;
    Real relTol = 1e-5;
    int i = 0;

    // allow mulitple projections on level 0 (can't do this on refined levels as it breaks the CF boundary conditions)
    int maxNumProj = 1;
    if (m_level == 0)
    {
      maxNumProj = m_opt.maxProjBaseLevel;
    }

    //    for (int i = 0; i < numProj; i++)
    while(! (maxDivU < initMaxDivU*relTol || maxDivU < relTol || i >= maxNumProj))
    {

      setVelZero(*m_vectorNew[uvar], m_opt.ccVelPorosityLimit);

      velBC.applyBCs(*m_vectorNew[uvar], m_grids, m_problem_domain,
                           m_dx, false); // inhomogeneous

      if (i > 0 && m_level==0)
      {
        // this only works on level 0 - we don't have any CF boundaries

          m_projection.AdditionalLevelProject(*m_vectorNew[uvar], crseVelPtr,
                                         new_time, a_dt,  m_parameters.isViscous());
      }
      else
      {
        m_projection.LevelProject(*m_vectorNew[uvar], crseVelPtr,
                                  new_time, a_dt, pressureScalePtr, crsePressureScalePtr, pressureScaleEdgePtr, crsePressureScaleEdgePtr,
                                  m_parameters.isViscous());
      }
      // as things stand now, physical BC's are re-set in LevelProjection

      //Get a copy of the pressure so we can write it out
      m_projection.unscaledPi(*m_scalarNew[ScalarVars::m_pressure], a_dt);

      // need to do physical boundary conditions and exchanges
      velBC.applyBCs(*m_vectorNew[uvar], m_grids, m_problem_domain,
                     m_dx, false); // inhomogeneous
      m_vectorNew[uvar]->exchange();

      QuadCFInterp interp;
      Divergence::levelDivergenceCC(*m_scalarNew[ScalarVars::m_divU], *m_vectorNew[uvar], NULL, m_dx, true, interp);

      // Actually want norm (sum of absolute values)
      maxDivU = computeNorm(*m_scalarNew[ScalarVars::m_divU], NULL, -1, m_dx, Interval(0,0), 0);

      pout() << "  CCProjection: max(div(U)) = " << maxDivU << endl;

      if (i == 0)
      {
        initMaxDivU = maxDivU;
      }

      i++;
    }



  } // end if doing projection


  // Make sure zero porosity regions have no velocity
  setVelZero(*m_vectorNew[uvar], m_opt.ccVelPorosityLimit);

  //Testing:
//  LevelData<FArrayBox> U_chi(m_grids, SpaceDim);
//  fillVectorField(U_chi, m_time, m_U_porosity, true, true);
//  Real max = ::computeMax(U_chi, NULL, 1, Interval(0, SpaceDim-1));
//  if (max > 1e5)
//  {
//    pout() << "U_chi > 1e5!!! " << endl;
//  }


//  QuadCFInterp interp;
//  Divergence::levelDivergenceCC(*m_scalarNew[ScalarVars::m_divU], *m_vectorNew[uvar], NULL, m_dx, true, interp);
//  Real maxDivU = computeNorm(*m_scalarNew[ScalarVars::m_divU], NULL, -1, m_dx, Interval(0,0), 0);

//  pout() << " CCProjection: max(div(U)) = " << maxDivU << endl;

  // Go back to standard porosity limit now we've finished CC solve
  //  m_lowerPorosityLimit = 1e-15;


  //cleanup
  if(crseVelPtr != NULL)
  {
    delete crseVelPtr;
    crseVelPtr = NULL;
  }
}


//This function taken from the Navier Stokes code
void AMRLevelMushyLayer::computeLapVel(LevelData<FArrayBox>& a_lapVel,
                                       LevelData<FArrayBox>& a_vel, const LevelData<FArrayBox>* a_crseVelPtr)
{

  // These are set in applyOpI anyway, but do it again here so we can see what happens
  VelBCHolder velBC(m_physBCPtr->uStarFuncBC(m_opt.viscousBCs));
  velBC.applyBCs(a_vel, a_vel.disjointBoxLayout(), m_problem_domain, m_dx, false);

  // Copy velocity to new data holder before taking laplacian to avoid changing
  // original velocity
  LevelData<FArrayBox> a_vel_copy(m_grids, SpaceDim, 2*IntVect::Unit); // 2 ghost vectors for higher order laplace
  DataIterator dit(m_grids);
  for (dit.begin(); dit.ok(); ++dit)
  {
    //          a_vel_copy[dit].copy(a_vel[dit]);
    Box b = m_grids[dit];
    b&=m_problem_domain.domainBox();
//    a_vel_copy[dit].setVal(0.0);
//    a_vel_copy[dit] += a_vel[dit];
//    a_vel_copy[dit].plus(a_vel[dit], b, 0, 0, SpaceDim);
//              int temp=0;
  }

  a_vel.copyTo(Interval(0, SpaceDim-1), a_vel_copy, Interval(0, SpaceDim-1));


  // Need 2nd order BCs for laplace operator
  if (a_crseVelPtr)
  {

    m_quadCFInterpVector.coarseFineInterp(a_vel_copy, *a_crseVelPtr);
    //          extrapBC.applyBCs(a_vel_copy, a_vel_copy.disjointBoxLayout(), m_problem_domain, m_dx, false);
  }

  // Apply domain BCs again?
  velBC.applyBCs(a_vel_copy, a_vel_copy.disjointBoxLayout(), m_problem_domain, m_dx, false);


  a_vel_copy.exchange(Interval(0, SpaceDim-1));

  for (int dir = 0; dir < SpaceDim; dir++)
  {

    // todo - Nonurgent: should be possible to replace this with the Darcy Brinkman Op, to save having the viscous op lying around
    m_viscousOp[dir]->applyOpI(a_lapVel, a_vel_copy, false);

    a_lapVel.exchange(a_lapVel.interval());

    // Do some smoothing
    Box domBox = m_problem_domain.domainBox();
    for (int dir =0; dir <SpaceDim; dir ++)
    {
      if (!m_problem_domain.isPeriodic(dir))
      {
        domBox.grow(dir, -1);
      }
    }

    for (int i=0; i < m_opt.lapVelNumSmooth; i++)
    {
      for (DataIterator dit = a_lapVel.dataIterator(); dit.ok(); ++dit)
      {
        Box b = a_lapVel[dit].box();
        b &= domBox;

        FArrayBox& thisLap = a_lapVel[dit];

        for (int dir=0; dir<SpaceDim; dir++)
        {

        for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
        {
          IntVect iv = bit();

          thisLap(iv, dir) = (1-m_opt.lapVelSmoothScale)*thisLap(iv, dir) + (m_opt.lapVelSmoothScale/4)*(thisLap(iv+BASISV(0), dir) +
              thisLap(iv-BASISV(0), dir) +
              thisLap(iv+BASISV(1), dir) +
              thisLap(iv-BASISV(1), dir));
        }
        }

      }
    }
  }
  // may need to extend lapVel to cover ghost cells as well
  {

    if (m_opt.lapVelBCOrder >= 0)
    {
      BCHolder bc = m_physBCPtr->extrapolationFuncBC(m_opt.lapVelBCOrder);
      const DisjointBoxLayout& grids = a_lapVel.getBoxes();
      DataIterator dit = a_lapVel.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
      {
        bc(a_lapVel[dit], grids[dit], m_problem_domain, m_dx, false); // not homogeneous
      }
    }
  }

  // finally, do exchange
  a_lapVel.exchange(a_lapVel.interval());
}


void AMRLevelMushyLayer::computePredictedVelocities(
    LevelData<FluxBox>& U_chi_new, LevelData<FArrayBox>& a_traceVel,
    LevelData<FluxBox>& a_advVel, LevelData<FArrayBox>& U_chi,
    const LevelData<FArrayBox>& a_viscousSource, PatchGodunov& a_patchGodVelocity,
    LevelData<FluxBox>& a_grad_eLambda, LevelData<FluxBox>& a_gradPhi,
    LevelData<FluxBox>& porosityFace,
    Real a_old_time, Real a_dt, bool legacyCompute)
{

  // Get AdvectionPhysics object within the PatchGodunov object
  AdvectionPhysics* advectionPhysics =
      dynamic_cast<AdvectionPhysics*>(a_patchGodVelocity.getGodunovPhysicsPtr());
  if (advectionPhysics == NULL)
  {
    MayDay::Error(
        "AMRLevelMushyLayer::predictVelocities - unable to upcast GodunovPhysics to AdvectionPhysics");
  }

  a_patchGodVelocity.setCurrentTime(a_old_time);


  const DisjointBoxLayout& levelGrids = a_traceVel.getBoxes();
  // loop over grids and do prediction
  DataIterator dit = levelGrids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    Box gridBox = levelGrids[dit()];
    FArrayBox& thisOldVel = a_traceVel[dit()];
    FArrayBox& U_chi_old = U_chi[dit()];
    FluxBox& thisAdvVel = a_advVel[dit()];
    FluxBox& U_chi_predicted = U_chi_new[dit()];
    const FArrayBox& srcFab = a_viscousSource[dit()];

    // set up patchGodunov for this patch
    a_patchGodVelocity.setCurrentBox(gridBox);
    advectionPhysics->setCellVelPtr(&thisOldVel);
    advectionPhysics->setAdvVelPtr(&thisAdvVel);

    // compute face-centered variables using "Godunov Box"
    a_patchGodVelocity.computeWHalf(U_chi_predicted, U_chi_old, srcFab,
                                    a_dt, gridBox);

    // Make sure U_chi_predicted is indeed U/chi
    if (m_opt.advectionMethod == m_porosityOutsideAdvection)
    {
      // we advected U so need to divide by porosity

      for (int velComp = 0; velComp < SpaceDim; velComp++)
      {
        for (int dir = 0; dir < SpaceDim; dir++)
        {
          FArrayBox& this_U_chi_Dir = U_chi_predicted[dir];
          FArrayBox& porosityDir = porosityFace[dit][dir];


          this_U_chi_Dir.divide(porosityDir, 0, velComp, 1);
        }
      }
    }


    // now loop over directions -- normal direction velocities
    // are copied from advection velocities then dividied by porosity.
    // for tangential direction we use face-centered predicted u/porosity
    // with the MAC correction calculated previously

    // u_chi_predicted (and this_U_chi_dir) should contain whatever we have advected (usually u/chi).
    // meanwhile, thisAdvVel (and thisAdvVelDir) contains whatever we were advecting with (usually u).

      for (int velComp = 0; velComp < SpaceDim; velComp++)
      {
        for (int dir = 0; dir < SpaceDim; dir++)
        {
          FArrayBox& this_U_chi_Dir = U_chi_predicted[dir];
          FArrayBox& thisAdvVelDir = thisAdvVel[dir];
          FArrayBox& porosityDir = porosityFace[dit][dir];

          if (dir == velComp)
          {
            // normal direction -- copy from advVel->uHalf
            // srcComp is always 0, since advVel only has
            // one component
            int srcComp = 0;
            int destComp = velComp;

            // Copy whatever advection velocity we just used to do advection
            // if m_porosityInAdvection, thisAdvVelDir = U/chi
            // if m_porosityOutsideAdvection, thisAdvVelDir = U/chi
            this_U_chi_Dir.copy(thisAdvVelDir, srcComp, destComp, 1);

            // now need to subtract off grad(eLambda)
            // this is the due to the difference between
            // _advecting_ and _advected_ velocities
            FluxBox& thisCorr = a_grad_eLambda[dit()];
            FArrayBox& thisCorrDir = thisCorr[dir];
            this_U_chi_Dir.minus(thisCorrDir, 0, dir, 1);

          } // end if dir is velcomp
          else
          {
            // add MAC correction to traced velocities so they are (roughly) divergence free
            // if we are advecting u/chi, need to divide correction by chi too

            // The MAC correction computed earlier in this timestep
            FArrayBox& gradPhiDir = a_gradPhi[dit()][dir];

            // Scale correction with chi too
            // as the correction was computed for U, and we're in fact correcting U/chi
            gradPhiDir.divide(porosityDir, 0, dir, 1);

            // subtract correction
            this_U_chi_Dir.minus(gradPhiDir, velComp,
                                 velComp, 1);

          } // end if tangential direction

        } // end loop over face directions

      } // end loop over velocity components

  } // end loop over grids
}



void AMRLevelMushyLayer::computeUstar(LevelData<FArrayBox>& a_UdelU,
                                      const LevelData<FArrayBox>& advectionSourceTerm,
                                      Real a_oldTime, Real a_dt, bool doFRupdates,
                                      bool a_MACprojection, bool compute_uDelU)
{
  CH_TIME("AMRLevelMushyLayer::computeUstar");


  // old u should be m_fluidVel, whether we're doing MAC projection or not
  int ustarVar = m_Ustar;

  if (a_MACprojection)
  {
    ustarVar = m_advUstar;
  }


  // Setup and solve poisson equation for Ustar
  Real old_time = a_oldTime; // m_time - m_dt;
  //  Real half_time = m_time - m_dt / 2;
  //  Real new_time = m_time; // define this just to be explicit in the code

  Real old_crseTime=-1, new_crseTime=-1;

  if (m_level>0)
  {
    AMRLevelMushyLayer* amrMLcrse = getCoarserLevel();
    new_crseTime = amrMLcrse->m_time;
    old_crseTime = new_crseTime - amrMLcrse->m_dt;
  }

  int srcGhost = 0;
  IntVect ivGhost = srcGhost * IntVect::Unit;
  LevelData<FArrayBox> src(m_grids, SpaceDim, ivGhost);

  DataIterator dit = m_grids.dataIterator();

  Real src_time = a_oldTime + m_opt.CCVelSrcTermCentering * a_dt;

  computeUstarSrc(src, advectionSourceTerm, src_time, a_MACprojection, compute_uDelU);

  if (!m_parameters.isViscous())
  {
    // explicit update.
    for (dit.reset(); dit.ok(); ++dit)
    {

      // Put old velocity into U^*
      fillVectorField(*m_vectorNew[ustarVar], old_time, m_fluidVel, true);

      src[dit].mult(a_dt);
      (*m_vectorNew[ustarVar])[dit].copy(src[dit]);


    }
  }  // end if inviscid
  else
  {



    /*
     * multi component solve not implemented yet. Currently we solve in each direction separately,
     * but it might be faster if the solver could do both simultaneously, to avoid lots of extra multigrid
     * coarsening and refinement. One to look at in the future.
     * TODO - Future: implement multi-directional U^* solver
     */
    if (m_opt.multiCompUStarSolve)
    {
      MayDay::Error("multiCompUStarSolve not implemented yet.");
      // Start by defining everything in a multi-comp way
//      LevelBackwardEuler UstarBE;
//      LevelTGA UstarTGA;
//
//      Vector<AMRLevelMushyLayer*> hierarchy;
//        Vector<DisjointBoxLayout> allGrids;
//        Vector<int> refRat;
//        ProblemDomain lev0Dom;
//        Real lev0Dx;
//        getHierarchyAndGrids(hierarchy, allGrids, refRat, lev0Dom, lev0Dx);
//
//        // Make these thingsmulticomp
//        defineUstarMultigrid();
//        UstarBE.define(allGrids, refRat, lev0Dom, s_uStarOpFactMultiComp, s_uStarAMRMGMultiComp);
//        UstarTGA.define(allGrids, refRat, lev0Dom, s_uStarOpFactMultiComp, s_uStarAMRMGMultiComp);

    }
    else
    {

      Vector<RefCountedPtr<LevelBackwardEuler> > UstarBE;
          Vector<RefCountedPtr<LevelTGA> > UstarTGA;
          defineUstarSolver(UstarBE, UstarTGA);

      // Do the solve for each component separately (as the solver is hardwired for one component only)
      LevelData<FArrayBox> UstarCrse, UoldCrse, *UstarCrseComp, *UoldCrseComp, Uold;
      LevelData<FluxBox> diffusiveFlux(m_grids, SpaceDim);
      LevelFluxRegister *fineFluxRegPtr=NULL, *crseFluxRegPtr=NULL;

      if (m_vectorFluxRegisters[VectorVars::m_fluidVel])
      {
        fineFluxRegPtr = &(*m_vectorFluxRegisters[VectorVars::m_fluidVel]);
      }

      Uold.define(m_grids, SpaceDim, m_numGhost*IntVect::Unit);
      fillVectorField(Uold, old_time, m_fluidVel, true);

      UstarCrseComp = NULL;
      UoldCrseComp = NULL;

      if (m_level > 0)
      {
        AMRLevelMushyLayer* amrMLcrse = getCoarserLevel();
        UstarCrse.define(amrMLcrse->m_grids, SpaceDim,
                         m_numGhost * IntVect::Unit);
        UstarCrseComp = new LevelData<FArrayBox>(amrMLcrse->m_grids, 1,
                                                 m_numGhost * IntVect::Unit);

        UoldCrse.define(amrMLcrse->m_grids, SpaceDim,
                        m_numGhost * IntVect::Unit);
        UoldCrseComp = new LevelData<FArrayBox>(amrMLcrse->m_grids, 1,
                                                m_numGhost * IntVect::Unit);

        // NB these used to be m_fluidVel
        //      amrMLcrse->fillVectorField(UstarCrse, new_crseTime, m_Ustar, true);
        //      amrMLcrse->fillVectorField(UoldCrse, old_crseTime, m_Ustar, true);

        // Coarse boundary conditions should be either
        //    a) the fully projected velocity on the coarse level in we're including grad(p) in the source
        // or b) U^*, if we're not including grad(p)

        int uvar = m_fluidVel;

        if (!m_addSubtractGradP)
        {
          uvar = m_Ustar;
        }

        if (a_MACprojection)
        {
          uvar = m_advectionVel;
        }
        amrMLcrse->fillVectorField(UstarCrse, new_crseTime, uvar, true);
        amrMLcrse->fillVectorField(UoldCrse, old_crseTime, uvar, true);

        crseFluxRegPtr = &(*(amrMLcrse->m_vectorFluxRegisters[VectorVars::m_fluidVel]));
      }


      // Set U^* to zero if it apears to be uninitialised
      Real maxUstar = 0;
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
      {
        for (int dir=0; dir<SpaceDim; dir++)
        {
          maxUstar = max(maxUstar, abs((*m_vectorNew[ustarVar])[dit].max(dir)));
        }
      }

      if (maxUstar > 1e100)
      {
        setValLevel(*m_vectorNew[ustarVar], 0.0);
      }


      // Solve each component (x, y, z) separately.
      pout() << "  U* solve ";

      for (int comp = 0; comp < SpaceDim; comp++)
      {
        Interval intvl(comp, comp);

        LevelData<FArrayBox> UstarComp;
        aliasLevelData(UstarComp, &(*m_vectorNew[ustarVar]), intvl);

        LevelData<FArrayBox> UoldComp;
        aliasLevelData(UoldComp, &Uold, intvl);

        LevelData<FArrayBox> compSrc;
        aliasLevelData(compSrc, &src, intvl);


        Vector<LevelData<FArrayBox>*> UstarVectComp;
        Vector<LevelData<FArrayBox>*> rhsVectComp;

        if (m_level > 0)
        {
          UstarCrse.copyTo(intvl, *UstarCrseComp, Interval(0, 0));
          UoldCrse.copyTo(intvl, *UoldCrseComp, Interval(0, 0));

          UstarVectComp.push_back(UstarCrseComp);
          rhsVectComp.push_back(NULL); //this isn't used
        }

        UstarVectComp.push_back(&UstarComp);
        rhsVectComp.push_back(&compSrc);

        if (!doFRupdates)
        {
          crseFluxRegPtr = NULL;
          fineFluxRegPtr = NULL;
        }

        int exitStatus = -1;
        Real resid = 0;

        if (m_opt.timeIntegrationOrder == 1)
        {
          UstarBE[comp]->updateSoln(UstarComp, UoldComp, compSrc,
                                    &(*fineFluxRegPtr), &(*crseFluxRegPtr),
                                    UoldCrseComp, UstarCrseComp,
                                    old_time,  old_crseTime, new_crseTime,
                                    a_dt, m_level, false, comp); // False - don't zero phi

#ifdef CH_FORK
          exitStatus = UstarBE[comp]->exitStatus();
          resid = UstarBE[comp]->finalResidual();
#endif

        }
        else
        {
          UstarTGA[comp]->updateSoln(UstarComp, UoldComp, compSrc,
                                     &(*fineFluxRegPtr), &(*crseFluxRegPtr),
                                     UoldCrseComp, UstarCrseComp,
                                     old_time, old_crseTime, new_crseTime,
                                     a_dt, m_level, false, comp);  // False - don't zero phi
#ifdef CH_FORK
          exitStatus = UstarTGA[comp]->exitStatus();
          resid = UstarTGA[comp]->finalResidual();
#endif
        }

        pout() << " Component " << comp << ": residual = " << resid;



      } // end loop over components

      pout () << endl;

      //Clean up
      if (UstarCrseComp != NULL)
      {
        delete UstarCrseComp;
        UstarCrseComp = NULL;
      }
      if (UoldCrseComp != NULL)
      {
        delete UoldCrseComp;
        UoldCrseComp = NULL;
      }


    } // end separate component solve


  } // end if viscous

}


void AMRLevelMushyLayer::computeAdvectionVelSourceTerm(LevelData<FArrayBox>& a_src)
{
  CH_TIME("AMRLevelMushyLayer::computeAdvectionVelSourceTerm");

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::computeAdvectionVelSourceTerm"            << endl;
  }


  IntVect ivGhost = a_src.ghostVect();
  LevelData<FArrayBox> darcy(m_grids, SpaceDim, ivGhost);
  LevelData<FArrayBox> viscous(m_grids, SpaceDim, ivGhost);
  LevelData<FArrayBox> buoyancy(m_grids, SpaceDim, ivGhost);
  LevelData<FArrayBox> pressure(m_grids, SpaceDim, ivGhost);
  LevelData<FArrayBox> temperature(m_grids, 1, ivGhost);
  LevelData<FArrayBox> liquidConc(m_grids, 1, ivGhost);
  LevelData<FArrayBox> porosity(m_grids, 1, ivGhost);
  LevelData<FArrayBox> pressureScale(m_grids, 1, ivGhost);

  //source term is at the old_time
  Real src_time = m_time - m_dt;

  LevelData<FArrayBox> velOld(m_grids, SpaceDim, m_numGhostAdvection * IntVect::Unit);
  LevelData<FArrayBox> permeability(m_grids, 1, m_numGhostAdvection * IntVect::Unit);

  calculatePermeability();
  fillScalars(permeability, src_time, m_permeability, true);

  fillScalars(temperature, src_time, m_temperature, true);
  fillScalars(liquidConc, src_time, m_liquidConcentration, true);
  fillScalars(porosity, src_time, m_porosity, true);
  fillScalars(pressureScale, src_time, m_pressureScaleVar, true);

  LevelData<FArrayBox> *crseVelPtr = NULL;
  AMRLevelMushyLayer* amrMushyLayerCoarserPtr;

  fillVectorField(velOld, src_time, m_fluidVel, true, true);

  if (m_level > 0)
  {
    amrMushyLayerCoarserPtr = getCoarserLevel();
    crseVelPtr = new LevelData<FArrayBox>(amrMushyLayerCoarserPtr->m_grids, SpaceDim, IntVect::Unit);
    amrMushyLayerCoarserPtr->fillVectorField(*crseVelPtr, src_time, m_fluidVel, true);
  }

  setVelZero(velOld, m_opt.advVelsrcChiLimit);

  //Apply BCs to grad(P) term - extrapolation
  DataIterator dit = m_grids.dataIterator();
  BCHolder bcExtrap = m_physBCPtr->extrapFuncBC();

  //  m_projection.gradPiBCs(pressure, true);
  fillPressureSrcTerm(pressure, pressureScale, src_time-m_dt/2,
                      false); // false - make sure we use Pi (true means use phi)


  // Make copy of these options in case we want to change them for just this timestep
  bool darcySrc = m_opt.advVelDarcySrc ;
  bool viscousSrc = m_opt.advVelViscousSrc;

  if (m_time <= m_opt.skipTrickySourceTerm)
  {
    pout() << "AMRLevelMushyLayer::computeAdvectionVelSourceTerm - time = " << m_time << " < " << m_opt.skipTrickySourceTerm;
    pout() << ", skipping darcy and viscous src terms";
    darcySrc = false;
    viscousSrc = false;
  }

  //Calculate Laplacian(U)
  if (m_parameters.isViscous())
  {
    bool recomputeLapVel = true;

    if (m_level > 0 && m_opt.advSrcAllowLaggedLapVel)
    {
      amrMushyLayerCoarserPtr = getCoarserLevel();
      Real crseOldTime = amrMushyLayerCoarserPtr->time() - amrMushyLayerCoarserPtr->dt();
      Real thisOldTime = m_time-m_dt ;

      // Don't recompute lap(u) if this is the second subcycled step - issues with CF BCs
      if (abs(thisOldTime- crseOldTime) > TIME_EPS)
      {
        recomputeLapVel = false;
        pout() << "Using old lap(U) in advection vel source term" << endl;
      }
    }

    if (recomputeLapVel)
    {
      computeLapVel(viscous, velOld, crseVelPtr);
      for (dit.reset(); dit.ok(); ++dit)
      {
        (*m_vectorNew[VectorVars::m_advSrcLapU])[dit].copy(viscous[dit]);
      }

    }
    else
    {
      for (dit.reset(); dit.ok(); ++dit)
      {
        viscous[dit].copy((*m_vectorNew[VectorVars::m_advSrcLapU])[dit]);
      }
    }

  }
  else
  {
    for (dit.reset(); dit.ok(); ++dit)
    {
      viscous[dit].setVal(0.0);
    }
  }

  fillBuoyancy(buoyancy, temperature, liquidConc, porosity);
  //  fillBuoyancy(buoyancy, src_time);



  for (dit.reset(); dit.ok(); ++dit)
  {
    darcy[dit].setVal(0.0);
    darcy[dit] += velOld[dit];

    // Divide each component of this term by the (scalar) permeability
    for (int dir = 0; dir < SpaceDim; dir++)
    {
      darcy[dit].divide(permeability[dit],
                        darcy[dit].box(), 0, dir);

      darcy[dit].mult(porosity[dit],
                      darcy[dit].box(), 0, dir);


    }

    darcy[dit].mult(m_parameters.m_darcyCoeff);
    viscous[dit].mult(m_parameters.m_viscosityCoeff);

  }

  // Finally, set darcy src term=0 near CF interface
  setVelZero(darcy, m_opt.advVelsrcChiLimit, 2); // set to zero up to 2 cells from interface

  for (dit.reset(); dit.ok(); ++dit)
  {
    a_src[dit].setVal(0.0);

    if (viscousSrc)
    {
      a_src[dit] += viscous[dit];
    }

    if (m_opt.advVelPressureSrc)
    {
      a_src[dit] -= pressure[dit];
    }

    if (m_opt.advVelBuoyancySrc)
    {
      a_src[dit] += buoyancy[dit];
    }

    if (darcySrc)
    {
      a_src[dit] -= darcy[dit];
    }



  } // end dataiterator

  // Add in extra (u.grad(porosity)/porosity^2)u term if needed
  if (m_opt.advectionMethod == m_porosityOutsideAdvection)
  {
    LevelData<FArrayBox> extraSrc(m_grids, SpaceDim, IntVect::Unit);

    Gradient::levelGradientCC(extraSrc, porosity, m_dx, extraSrc.ghostVect()[0]);

    for (DataIterator dit = extraSrc.dataIterator(); dit.ok(); ++dit)
    {
      for (int dir=0; dir < SpaceDim; dir++)
      {
        // Make u.grad(porosity)
        extraSrc[dit].mult(velOld[dit], dir, dir);

        // Make u.grad(porosity)/porosity^2
        extraSrc[dit].divide(porosity[dit], 0, dir);
        extraSrc[dit].divide(porosity[dit], 0, dir);

        // Make (u.grad(porosity)/porosity^2)u
        extraSrc[dit].mult(velOld[dit], dir, dir);
      }

      // Add this extra source term to the full source term
      a_src[dit].plus(extraSrc[dit], 0, 0, SpaceDim);
    }
  }
  else if (m_opt.advectionMethod == m_porosityInAdvection)
  {
    // Need a different source term for this problem

    LevelData<FArrayBox> extraSrc(m_grids, SpaceDim, IntVect::Unit);

    for (DataIterator dit = extraSrc.dataIterator(); dit.ok(); ++dit)
    {

      extraSrc[dit].setVal(0.0);

      for (int dir=0; dir < SpaceDim; dir++)
      {
        extraSrc[dit].plus(m_dPorosity_dt[dit], 0, dir);

        // Make u (dchi/dt)
        extraSrc[dit].mult(velOld[dit], dir, dir);

        // Make u (dchi/dt)/porosity^2
        extraSrc[dit].divide(porosity[dit], 0, dir);
        extraSrc[dit].divide(porosity[dit], 0, dir);

        // Also need to scale the original source term
        a_src[dit].divide(porosity[dit], 0, dir);
      }

      // Add this extra source term to the full source term
      a_src[dit].minus(extraSrc[dit], 0, 0, SpaceDim);

    }
  }

  Real maxSrcVal = ::computeNorm(a_src, NULL, 1, m_dx);

  if (maxSrcVal > 1e10)
  {
    pout() << "WARNING - norm of advection velocity src term = " << maxSrcVal << endl;
  }

  setVelZero(a_src, m_opt.advVelsrcChiLimit);

  a_src.exchange();

  for (dit.reset(); dit.ok(); ++dit)
  {
    // This source term is centered at step n, so store it in vectorOld.
    (*m_vectorOld[VectorVars::m_advectionSrc])[dit].copy(a_src[dit], 0, 0, SpaceDim);

    // Also need to store in vectorNew so that we can interpolate in time when subcycling
    (*m_vectorNew[VectorVars::m_advectionSrc])[dit].copy(a_src[dit], 0, 0, SpaceDim);
  }

  if (crseVelPtr != NULL)
  {
    delete crseVelPtr;
    crseVelPtr = NULL;
  }

}



void AMRLevelMushyLayer::computeUDelU(LevelData<FArrayBox>& U_adv_src, const LevelData<FArrayBox>& advectionSourceTerm, Real a_oldTime, Real a_dt)
{
  // May want to use pre existing value for this, rather than re computing.


  LevelData<FArrayBox> UdelU_porosity(m_grids, SpaceDim, m_numGhostAdvection*IntVect::Unit);
  DataIterator dit = UdelU_porosity.dataIterator();

  if (m_opt.doEulerPart)
  {

    // option 0 - use advection velocity to upwind old velocities to cell faces, and use these in calculation
    // option 1 - average advection velocities to cell centers, and use this along with gradient of advection velocities to computer u.del(u)
    if (m_opt.uDeluMethod == 0)
    {
      bool doVelFRupdates = true;
      predictVelocities(UdelU_porosity, m_advVel, advectionSourceTerm, a_oldTime, a_dt, doVelFRupdates);

    }
    else if (m_opt.uDeluMethod == 1)
    {


      LevelData<FArrayBox> ccVel(m_grids, SpaceDim, IntVect::Unit);
      LevelData<FluxBox> vel_chi(m_grids, 1, IntVect::Unit);
      LevelData<FluxBox> porosityFace(m_grids, 1, IntVect::Unit);

      fillScalarFace(porosityFace, a_oldTime, m_porosity, true);
      for (DataIterator dit = UdelU_porosity.dataIterator(); dit.ok(); ++dit)
      {
        vel_chi[dit].copy(m_advVel[dit]);
        for (int dir=0; dir<SpaceDim; dir++)
        {
          vel_chi[dit].divide(porosityFace[dit], porosityFace[dit].box(), 0,0);
        }
      }

      m_advVel.exchange();
      EdgeToCell(m_advVel, ccVel);
      ccVel.exchange();

      Gradient::levelGradientCC(UdelU_porosity, vel_chi, m_dx);
      for (DataIterator dit = UdelU_porosity.dataIterator(); dit.ok(); ++dit)
      {
        UdelU_porosity[dit].mult(ccVel[dit]);
      }

      Real maxUdelU = ::computeMax(UdelU_porosity, NULL,1, Interval(0, SpaceDim-1));

      if (maxUdelU > 1e3)
      {
        pout() << "max(uDelU) = " << maxUdelU << endl;
      }


    }


  }
  else
  {
    // Don't want this term
    for (dit.reset(); dit.ok(); ++dit)
    {
      UdelU_porosity[dit].setVal(0.0);
    }
  }


  setVelZero(UdelU_porosity, m_opt.uDelU_porosityLimit, m_opt.uDelU_grow); // grow by 1

  UdelU_porosity.copyTo(U_adv_src);

  for (dit.reset(); dit.ok(); ++dit)
  {
    (*m_vectorNew[VectorVars::m_UdelU])[dit].copy(UdelU_porosity[dit], 0, 0, SpaceDim);
  }

}


void AMRLevelMushyLayer::computeUstarSrc(LevelData<FArrayBox>& src,
                                         const LevelData<FArrayBox>& advectionSourceTerm,
                                         Real src_time,
                                         bool a_MACprojection, bool compute_uDelU)
{
  // Get coarse pointer if necessary
  LevelData<FArrayBox> crseOldVel;
  AMRLevelMushyLayer* amrMLcrse;
  int srcGhost = src.ghostVect()[0];

  bool pressureSrc = m_addSubtractGradP;
  bool advSrc = m_opt.CCAdvSrc;

  if (m_opt.CCPressureSrcOverride)
  {
    pressureSrc = m_opt.CCPressureSrc;
  }

  if (m_time <= m_opt.skipTrickySourceTerm)
  {
    pout() << "AMRLevelMushyLayer::computeUstarSrc - time = " << m_time << " < " << m_opt.skipTrickySourceTerm;
    pout() << ", skipping advection src term";
    advSrc = false;
  }

  Real oldTime = m_time - m_dt;

  if (m_level > 0)
  {
    amrMLcrse = getCoarserLevel();

    crseOldVel.define(amrMLcrse->m_grids, SpaceDim,
                      srcGhost * IntVect::Unit);
    amrMLcrse->fillVectorField(crseOldVel, oldTime, m_fluidVel, true);

    //          nRefCrse = amrMLcrse->m_ref_ratio;
  }
  DataIterator dit = m_grids.dataIterator();

  // Create data structures
  IntVect ivGhost = srcGhost * IntVect::Unit;
  LevelData<FArrayBox> U_src(m_grids, SpaceDim, ivGhost); // = (1 - dt pi_0 chi / (2*pi)  +(Da dt/2) laplacian) u^n
  LevelData<FArrayBox> U_nonviscous_src(m_grids, SpaceDim, ivGhost); //
  LevelData<FArrayBox> P_src(m_grids, SpaceDim, ivGhost); // = chi * grad(P^{n-1/2})
  LevelData<FArrayBox> buoyancy_src(m_grids, SpaceDim, ivGhost); // = chi*Ra*(theta^n+theta^{n+1})/2 e_z
  LevelData<FArrayBox> darcy_src(m_grids, SpaceDim, ivGhost); // = (pr/Da) * chi*U/Pi
  LevelData<FArrayBox> U_correction(m_grids, SpaceDim, ivGhost); //= (1-acoef)*U_old to correct for how TGA and backward euler work

  // Create time centred versions of scalar fields
  LevelData<FArrayBox> porosity_centred(m_grids, 1, ivGhost);
  LevelData<FArrayBox> permeability_centred(m_grids, 1, ivGhost);
  LevelData<FArrayBox> temperature_centred(m_grids, 1, ivGhost);
  LevelData<FArrayBox> liquidConcentration_centred(m_grids, 1, ivGhost);
  LevelData<FArrayBox> pressureScale(m_grids, 1, ivGhost);

  LevelData<FArrayBox> U_half(m_grids, SpaceDim, ivGhost);
  fillVectorField(U_half, src_time, m_fluidVel, true);

  LevelData<FArrayBox> U_old(m_grids, SpaceDim, ivGhost);
  fillVectorField(U_old, oldTime, m_fluidVel, true);
  //  LevelData<FArrayBox> U_scaled(m_grids, SpaceDim, ivGhost);
  //  fillVectorField(U_scaled, oldTime, m_fluidVel, true);

  LevelData<FArrayBox> U_adv_src(m_grids, SpaceDim, ivGhost);
  setValLevel(U_adv_src, 0.0);
  if (compute_uDelU)
  {
    computeUDelU(U_adv_src, advectionSourceTerm, oldTime, m_dt);
  }


  fillScalars(temperature_centred,         src_time, m_temperature, true);
  fillScalars(porosity_centred,            src_time, m_porosity, true);
  fillScalars(permeability_centred,        src_time, m_permeability, true);
  fillScalars(liquidConcentration_centred, src_time, m_liquidConcentration, true);
  fillScalars(pressureScale, src_time, m_pressureScaleVar, true);

  // Now the Pressure bit
  fillPressureSrcTerm(P_src, pressureScale, oldTime-m_dt/2, a_MACprojection);

  // Buoyancy
  fillBuoyancy(buoyancy_src, temperature_centred, liquidConcentration_centred,porosity_centred);
  //  fillBuoyancy(buoyancy_src, src_time);

  // Finally the buoyancy bit
  for (dit.reset(); dit.ok(); ++dit)
  {

    //    fillBuoyancy(buoyancy_src[dit], temperature_centred[dit], liquidConcentration_centred[dit],
    //                 porosity_centred[dit],  (*m_vectorNew[VectorVars::m_bodyForce])[dit]);


    darcy_src[dit].copy(U_old[dit]);

    for (int dir = 0; dir < SpaceDim; dir++)
    {
      darcy_src[dit].mult(porosity_centred[dit], 0, dir, 1);
      darcy_src[dit].divide(permeability_centred[dit], darcy_src[dit].box(), 0, dir, 1);
    }
    darcy_src[dit].mult(m_parameters.m_darcyCoeff);

  }

  //Put it all together
  for (dit.reset(); dit.ok(); ++dit)
  {

    src[dit].setVal(0.0);

    if (m_opt.CCBuoyancySrc)
    {
      src[dit] += buoyancy_src[dit];
    }

    if (advSrc)
    {
      src[dit] -= U_adv_src[dit];
    }

    if (m_addSubtractGradP || pressureSrc)
    {
      src[dit] -= P_src[dit];
    }

    if (m_opt.explicitDarcyTerm || m_opt.CCDarcySrc)
    {
      src[dit] -= darcy_src[dit];
    }

  }

  src.exchange();

  src.copyTo(*m_vectorNew[VectorVars::m_viscousSolveSrc]);

}


void AMRLevelMushyLayer::predictVelocities(LevelData<FArrayBox>& a_uDelU,
                                           LevelData<FluxBox>& a_advVel,
                                           const LevelData<FArrayBox>& a_src,
                                           Real old_time,
                                           Real a_dt,
                                           bool doFRupdates)
{
  CH_TIME("AMRLevelMushyLayer::predictVelocities()");

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::predictVelocities: " << m_level << endl;
  }

  const DisjointBoxLayout& levelGrids = a_advVel.getBoxes();

  // for tracing, will need to fill in boundary values for
  // grown copy of velocity
  IntVect advectionGhostVect = m_numGhostAdvection * IntVect::Unit;

  LevelData<FArrayBox> traceVel(levelGrids, SpaceDim, advectionGhostVect);
  LevelData<FArrayBox> UtoAdvect_old(m_grids, SpaceDim, advectionGhostVect);

  // Advection velocity is either U or U/chi depending on how we're doing things
  LevelData<FluxBox> advectionVelocity(levelGrids, 1, advectionGhostVect);

  LevelData<FArrayBox> porosity(levelGrids, 1, advectionGhostVect);
  LevelData<FluxBox> porosityFace(levelGrids, 1, advectionGhostVect);

  fillScalars(porosity, old_time, m_porosity, true);
  fillScalarFace(porosityFace, old_time, m_porosity, true);

  // now fill in velocity field
  // m_time is new time
  fillVectorField(traceVel, old_time, m_fluidVel, true);

  // tracing uses inviscid BC's
  {
    VelBCHolder velBC(m_physBCPtr->tracingVelFuncBC());
    velBC.applyBCs(traceVel, levelGrids, m_problem_domain, m_dx, false); // not homogeneous
  }

  for (DataIterator dit = a_advVel.dataIterator(); dit.ok(); ++dit)
  {
    advectionVelocity[dit].copy(a_advVel[dit]);

    if ( m_opt.advectionMethod == m_noPorosity )
    {
      // do nothing
    }
    else if ( m_opt.advectionMethod == m_porosityOutsideAdvection
        || m_opt.advectionMethod == m_porosityInAdvection)
    {

      advectionVelocity[dit].divide(porosityFace[dit], porosityFace[dit].box(), 0, 0);
      for (int dir=0; dir<SpaceDim; dir++)
      {
        traceVel[dit].divide(porosity[dit], porosity[dit].box(), 0, dir);
      }
    }

    for (int dir=0; dir<SpaceDim; dir++)
    {
      advectionVelocity[dit][dir].mult(m_parameters.m_advectionCoeff);
    }
  }

  if (m_opt.advectionMethod == m_porosityInAdvection)
  {
    fillVectorField(UtoAdvect_old, old_time, m_U_porosity, true);
  }
  else if (m_opt.advectionMethod == m_porosityOutsideAdvection ||  m_opt.advectionMethod == m_noPorosity)
  {
    fillVectorField(UtoAdvect_old, old_time, m_fluidVel, true);
  }


  setVelZero(UtoAdvect_old, m_opt.advVelChiLimit);
  setVelZero(advectionVelocity, m_opt.advVelChiLimit);

  // will need edge-centered storage for all velocity components
  // This is the half time predicted value of whatever we were advecting (u or u/chi)
  LevelData<FluxBox> U_chi_advected(levelGrids, SpaceDim);

  { //this is to allow this fluxbox to go out of scope
    // will also need grad(eLambda) and grad(phi_mac)
    LevelData<FluxBox>& grad_eLambda = m_projection.grad_eLambda();
    LevelData<FluxBox> gradPhi(levelGrids, SpaceDim);
    m_projection.gradPhi(gradPhi);

    // Need to multiply gradPhi by porosity to get correct correction field
    LevelData<FluxBox> pressureScale(levelGrids, 1);
    fillScalarFace(pressureScale, old_time, m_pressureScaleVar, true);


    // Get the true pressure correction (scaled with chi if appropriate)
    if (m_opt.scaleP_MAC)
    {
      for (DataIterator dit=gradPhi.dataIterator(); dit.ok(); ++dit)
      {
        for (int dir=0; dir<SpaceDim; dir++)
        {
          gradPhi[dit].mult(pressureScale[dit], pressureScale[dit].box(), 0, dir);
        }
      }
    }

    traceVel.exchange();
    advectionVelocity.exchange();
    UtoAdvect_old.exchange();
    grad_eLambda.exchange();
    gradPhi.exchange();
    pressureScale.exchange();

    // This should U_chi_advected - if we've only advected U, then it should divide by chi before returning
    computePredictedVelocities(U_chi_advected, traceVel, advectionVelocity,
                               UtoAdvect_old, a_src, m_patchGodVelocity,
                               grad_eLambda, gradPhi, pressureScale, old_time, a_dt,
                               m_opt.legacyComputePredictVel);

    U_chi_advected.exchange();

    for (DataIterator dit = U_chi_advected.dataIterator(); dit.ok(); ++dit)
    {
      for (int dir=0; dir<SpaceDim; dir++)
      {
        EdgeToCell(U_chi_advected[dit], dir, (*m_vectorNew[VectorVars::m_U_porosity])[dit], 0, dir);
      }
    }

  } // FluxBox goes out of scope, memory reclaimed.

  // now loop over grids to compute advection and update flux registers
  DataIterator dit = a_uDelU.dataIterator();

  // for nonconservative update, also need a cell-centered
  // advection velocity -- get this by averaging edge-centered
  // advection velocity to cell centers.  this one shouldn't need
  // any ghost cells.

  // reduce, reuse, recycle...
  // cellAdvVel is U_advection averaged to cells. Never any factors of chi in here.
  LevelData<FArrayBox> cellAdvVel(m_grids, SpaceDim);
  EdgeToCell(a_advVel, cellAdvVel);

  // Compute flux for flux registers - U*(U/chi)
  // Also might want this to compute div(flux) as the source term
  LevelData<FluxBox> momentumFlux(U_chi_advected.disjointBoxLayout(), SpaceDim, U_chi_advected.ghostVect());
  U_chi_advected.copyTo(Interval(0, SpaceDim-1), momentumFlux, Interval(0, SpaceDim-1));

  for (dit.reset(); dit.ok(); ++dit)
  {
    for (int dir = 0; dir < SpaceDim; dir++)
    {
      FArrayBox& thisFlux = momentumFlux[dit][dir];
      FArrayBox& thisAdvVel = a_advVel[dit][dir];

      // now need to loop over traced velocity components
      for (int velComp = 0; velComp < SpaceDim; ++velComp)
      {
        thisFlux.mult(thisAdvVel, 0, velComp, 1);
      }
    }
  }

  if (m_opt.uDelUConservativeForm)
  {
    Divergence::levelDivergenceMACMultiComp(a_uDelU, momentumFlux, m_dx);
  }
  else
  {

    for (dit.reset(); dit.ok(); ++dit)
    {
      FluxBox& thisU_chiHalf = U_chi_advected[dit];
      FArrayBox& thisCellAdvVel = cellAdvVel[dit];
      const Box& thisBox = levelGrids[dit];
      FArrayBox& this_uDelU = a_uDelU[dit];
      this_uDelU.setVal(0.0);

      // to do this in a dimensionality independent way,
      // loop over directions.
      for (int dir = 0; dir < SpaceDim; dir++)
      {
        FArrayBox& U_chi_HalfDir = thisU_chiHalf[dir];

        // this will increment this_uDelU with
        // cellVel*d(uHalfDir)/d[xyz] for all velocity components
        FORT_UDELS(CHF_FRA(this_uDelU), CHF_FRA(U_chi_HalfDir),
                   CHF_FRA(thisCellAdvVel), CHF_BOX(thisBox), CHF_REAL(m_dx),
                   CHF_INT(dir));


      } // end loop over directions
    } // end loop over dataiterator
  } // end if conservative/non conservative calculation


  for (dit.reset(); dit.ok(); ++dit)
  {


    // now need to do flux register stuff
    if (doFRupdates)
    {

      Interval velComps(0, SpaceDim - 1);
      int reflux_var = m_fluidVel;

      Real fluxMult = m_dt; // petermc made positive, 7 Dec 07

      // was m_hasFiner, but only care about finer levels if they're defined,
      if (!finestLevel())
      {
        CH_assert(m_vectorFluxRegisters[reflux_var]->isDefined());
        for (int dir = 0; dir < SpaceDim; dir++)
        {
          // because we eventually want -D_R(U-del-U), signs
          // on fluxMult are opposite of what you might expect...
          if (m_opt.reflux_normal_momentum)
          {
            m_vectorFluxRegisters[reflux_var]->incrementCoarse(
                momentumFlux[dit][dir], fluxMult, dit(), velComps,
                velComps, dir);
          } else {
            // if we're not refluxing normal momentum component,
            // then do this component by component
            for (int velComp = 0; velComp < SpaceDim; velComp++)
            {
              if (velComp != dir)
              {
                Interval velInt(velComp, velComp);
                m_vectorFluxRegisters[reflux_var]->incrementCoarse(
                    momentumFlux[dit][dir], fluxMult, dit(),
                    velInt, velInt, dir);
              } // end if a tangential velocity component
            } // end loop over velocity components
          } // end if we're not refluxing normal momentum
        } // end loop over directions
      } // end if a finer level exists

      // If we have a coarser level...
      if (m_level > 0)
      {
        LevelFluxRegister& crseFR =
            *getCoarserLevel()->m_vectorFluxRegisters[reflux_var];
        CH_assert(crseFR.isDefined());
        // because we will eventually want -D_R(U-del-U) signs
        // on fluxMult are opposite of what you might expect...

        for (int dir = 0; dir < SpaceDim; dir++)
        {
          if (m_opt.reflux_normal_momentum)
          {
            crseFR.incrementFine(momentumFlux[dit][dir], fluxMult, dit(),
                                 velComps, velComps, dir, Side::Lo);

            crseFR.incrementFine(momentumFlux[dit][dir], fluxMult, dit(),
                                 velComps, velComps, dir, Side::Hi);
          } else {
            // for case where we only want to reflux tangential
            // momentum components, do this component by component
            // so that we can isolate the normal component.
            for (int velComp = 0; velComp < SpaceDim; velComp++)
            {
              if (velComp != dir)
              {
                Interval velInt(velComp, velComp);

                crseFR.incrementFine(momentumFlux[dit][dir], fluxMult,
                                     dit(), velInt, velInt, dir, Side::Lo);

                crseFR.incrementFine(momentumFlux[dit][dir], fluxMult,
                                     dit(), velInt, velInt, dir, Side::Hi);
              } // end if this is a tangential component
            } // end loop over velocity components
          } // end if we're not refluxing normal momentum
        } // end loop over directions
      } // end if level > 0
    } // end loop over grids
    // at this point, we've computed the advective term, and updated
    // the flux registers with the momentum flux

  } // end if do flux register updates

}

void AMRLevelMushyLayer::traceAdvectionVel(LevelData<FluxBox>& a_advVel,
                                           LevelData<FArrayBox>& a_old_vel, LevelData<FArrayBox>& U_chi,
                                           const LevelData<FArrayBox>& a_viscousSource,
                                           PatchGodunov& a_patchGodVelocity, Real a_old_time, Real a_dt)
{
  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::traceAdvectionVel()" << endl;
  }

  CH_TIME("AMRLevelMushyLayer::traceAdvectionVel()");

  // Get AdvectionPhysics object within the PatchGodunov object
  AdvectionPhysics* advectionPhysics =
      dynamic_cast<AdvectionPhysics*>(a_patchGodVelocity.getGodunovPhysicsPtr());

  if (advectionPhysics == NULL)
  {
    MayDay::Error(
        "AMRLevelMushyLayer::traceAdvectionVel - unable to upcast GodunovPhysics to AdvectionPhysics");
  }

  a_patchGodVelocity.setCurrentTime(a_old_time);

  // this should now be unnecessary
  // also need to fill grown advection velocity
  // by cell-to-edge averaging old-time velocity
  //    CellToEdge(a_old_vel, a_advVel);
  CornerCopier cornerCopy(a_advVel.disjointBoxLayout(), a_advVel.disjointBoxLayout(), m_problem_domain, a_advVel.ghostVect(), true);
  a_advVel.exchange(cornerCopy); // definitely need this

  CornerCopier cornerCopy2(U_chi.disjointBoxLayout(), U_chi.disjointBoxLayout(), m_problem_domain, U_chi.ghostVect(), true);
  U_chi.exchange(cornerCopy2); // may need this

  // loop over grids and predict face-centered velocities at the half
  // time using the patchGodunov infrastructure
  const DisjointBoxLayout& levelGrids = a_advVel.getBoxes();
  DataIterator dit = a_advVel.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    //    Box gridBox = levelGrids[dit()];

    FArrayBox& thisOldVel = a_old_vel[dit()];
    FluxBox& thisAdvVel = a_advVel[dit()];
    FArrayBox& U_chi_old = U_chi[dit];
    const FArrayBox& srcFab = a_viscousSource[dit()];

    Box gridBox(levelGrids[dit()]);
    //    gridBox.grow(1); // should we be doing this?

    // temporary storage for the advected velocities
    FluxBox U_chi_new(gridBox, SpaceDim); //srcFab.box()

    // set up PatchGodunov and AdvectionPhysics for this patch
    a_patchGodVelocity.setCurrentBox(gridBox);
    advectionPhysics->setCellVelPtr(&thisOldVel); // u^n
    advectionPhysics->setAdvVelPtr(&thisAdvVel); // advection velocity
    //      advectionPhysics->setVelocities(&thisOldVel, &thisAdvVel);

    // predict face-centered U/chi
    a_patchGodVelocity.computeWHalf(U_chi_new, U_chi_old, srcFab, a_dt,
                                    gridBox);

    // copy the normal component of the U_chi_new into "thisAdvVel"
    // (the rest is discarded)
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      int srcComp = idir;
      int destComp = 0;
      int numComp = 1;

      thisAdvVel[idir].copy(U_chi_new[idir], srcComp, destComp, numComp);
    }

    EdgeToCell(a_advVel, *m_vectorNew[VectorVars::m_advectionVel]);
    Real maxAdvU = ::computeNorm(*m_vectorNew[VectorVars::m_advectionVel], NULL, 1, m_dx, Interval(0,0), 0);

    if (s_verbosity >= 1 && maxAdvU > 1e10)
    {
      pout() << "WARNING - max advection velocity (post tracing, pre-projection) = " << maxAdvU << endl;
    }

  } // end loop over grids

  //  int temp = 0;
}

void AMRLevelMushyLayer::computeInflowOutflowAdvVel()
{
  for (DataIterator dit = m_totalAdvVel.dataIterator(); dit.ok(); ++dit)
  {
    for (int idir = 0; idir<SpaceDim; idir++)
    {

      m_totalAdvVel[dit][idir].copy(m_advVel[dit][idir]);

      m_totalAdvVel[dit][idir].plus(m_frameAdvVel[dit][idir]);
    }
  }
}



void AMRLevelMushyLayer::correctEdgeCentredVelocity(LevelData<FluxBox>& a_advVel, Real a_dt)
{
  CH_TIME("AMRLevelMushyLayer::correctEdgeCentredVelocity()");
  if (s_verbosity >= 4)
  {
    pout() << "AMRLevelMushyLayer::correctEdgeCentredVelocity()" << endl;
  }


  Real old_time = m_time - m_dt;
  Real half_time = old_time + 0.5*m_dt;

  RefCountedPtr<LevelData<FluxBox> > crsePressureScaleEdgePtr, pressureScaleEdgePtrOneGhost;
  RefCountedPtr<LevelData<FArrayBox> > pressureScalePtr, crsePressureScalePtr;

  EdgeVelBCHolder edgeVelBC(m_physBCPtr->advectionVelFuncBC(m_opt.viscousBCs));
  edgeVelBC.applyBCs(a_advVel, m_grids, m_problem_domain, m_dx,
                     false); // inhomogeneous
  if (s_verbosity >= 6)
  {
    pout() << "  AMRLevelMushyLayer::correctEdgeCentredVelocity() - applied init BCs" << endl;
  }

  if(m_opt.scaleP_MAC)
  {
    pressureScaleEdgePtrOneGhost = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids, 1, IntVect::Unit));
    fillScalarFace(*pressureScaleEdgePtrOneGhost, half_time, m_pressureScaleVar, true);

    horizontallySmooth(*pressureScaleEdgePtrOneGhost);

    pressureScalePtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grids,1,IntVect::Unit));
    fillScalars(*pressureScalePtr, half_time, m_pressureScaleVar, true);

    if (m_level > 0)
    {
      if (s_verbosity >= 6)
      {
        pout() << "  AMRLevelMushyLayer::correctEdgeCentredVelocity() - get CF boundary conditions" << endl;
      }

      AMRLevelMushyLayer* mlCrse = getCoarserLevel();

      if (s_verbosity >= 6)
      {
        pout() << "  AMRLevelMushyLayer::correctEdgeCentredVelocity() - got coarse mushy layer object" << endl;
      }

      crsePressureScaleEdgePtr = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(mlCrse->m_grids, 1, IntVect::Unit));
      mlCrse->fillScalarFace(*crsePressureScaleEdgePtr, half_time, m_pressureScaleVar, true);

      crsePressureScalePtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(mlCrse->m_grids, 1, IntVect::Unit));
      mlCrse->fillScalars(*crsePressureScalePtr, half_time, m_pressureScaleVar, true);
    }

    if (s_verbosity >= 6)
    {
      pout() << "  AMRLevelMushyLayer::correctEdgeCentredVelocity() - got pressure scale vars" << endl;
    }
  }

  a_advVel.exchange();
  EdgeToCell(a_advVel, *m_vectorNew[VectorVars::m_advUstar]);

  if (m_opt.doProjection)
  {
    if (s_verbosity >= 6)
    {
      pout() << "  AMRLevelMushyLayer::correctEdgeCentredVelocity() - do projection" << endl;
    }

    int projNum = 0;


    //    Divergence::levelDivergenceMAC(*m_scalarNew[ScalarVars::m_divUadv], a_advVel, m_dx);
    Real maxDivU = 10*m_opt.maxDivUFace ; //::computeNorm(*m_scalarNew[ScalarVars::m_divUadv], NULL, 1, m_dx, Interval(0,0));

    while ( (maxDivU > m_opt.maxDivUFace && projNum < m_opt.maxNumMACProj) ||
        (m_opt.enforceAnalyticVel && projNum < 1) ) // do at least one projection of analytic vel
    {
      projNum++;

      int exitStatus = m_projection.levelMacProject(a_advVel, old_time, a_dt, pressureScalePtr, crsePressureScalePtr,
                                   pressureScaleEdgePtrOneGhost, crsePressureScaleEdgePtr);


      Divergence::levelDivergenceMAC(*m_scalarNew[ScalarVars::m_divUadv], a_advVel, m_dx);

      maxDivU = ::computeNorm(*m_scalarNew[ScalarVars::m_divUadv], NULL, 1, m_dx, Interval(0,0), 0);
      pout() << "  MAC Projection (#" << projNum << " on level "<< m_level << "), exit status = " << exitStatus << ", max(div u) = " << maxDivU << endl;


    }

  }

  //  Divergence::levelDivergenceMAC(*m_scalarNew[ScalarVars::m_divUadv], a_advVel, m_dx);

  // Fill CF boundaries
  if (m_level > 0)
  {
    fillAdvVel(old_time, a_advVel);
  }

  //  m_projection.applyFreestreamCorrection(a_advVel, m_dt);
  //  if (!freestreamBeforeProjection)
  //  {
  if (currentCFLIsSafe())
  {
    m_projection.applyFreestreamCorrection(a_advVel);
  }
  //  }

  a_advVel.exchange();

}
