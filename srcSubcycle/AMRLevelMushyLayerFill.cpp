#include "AMRLevelMushyLayer.H"

void AMRLevelMushyLayer::computeViscosity()
{
  LevelData<FArrayBox> &a_viscosity = *m_scalarNew[m_viscosity];
  const LevelData<FArrayBox> &a_liquid_concentration =
      *m_scalarNew[m_liquidConcentration];

  switch (m_parameters.m_viscosityFunction)
  {
  case ViscosityFunction::linearViscosity:

    // Viscosity varies linearly with increasing solute concentration from 1.0
    // to max_viscosity

    for (DataIterator dit = a_viscosity.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox &viscosity = a_viscosity[dit];

      // Liquid concentration is between -composition ratio and 0
      viscosity.copy(a_liquid_concentration[dit]);

      // This will make viscosity between 0 and 1.0
      viscosity.plus(m_parameters.compositionRatio);
      viscosity.mult(1.0 / m_parameters.compositionRatio);

      // This makes viscosity between 1.0 and max_viscosity
      viscosity.mult(m_parameters.max_viscosity - 1.0);
      viscosity.plus(1.0);
    }

    break;

  default:

    // Default is uniform viscosity = 1
    setValLevel(a_viscosity, 1.0);

    break;
  }
}

void AMRLevelMushyLayer::fillAdvVel(Real time, LevelData<FluxBox> &a_advVel)
{
  // Fill interior ghost cells
  if (m_level > 0)
  {
    int interpRadius = a_advVel.ghostVect()[0];

    AMRLevelMushyLayer *amrMLcrse = getCoarserLevel();

    // Refinement ratio between this and coarser level
    int crseRefRatio = amrMLcrse->m_ref_ratio;

    ProblemDomain crseDomain = amrMLcrse->m_problem_domain;
    DisjointBoxLayout crseGrids(amrMLcrse->m_grids);

    PiecewiseLinearFillPatchFace patcher(m_grids, crseGrids, 1, crseDomain,
                                         crseRefRatio, interpRadius);

    Real timeInterpCoeff = 0;

    // As advection velocities are stored at the half time step, we can't
    // necessarilly interpolate at exactly the right point in time, so accept
    // that we're making an O(dt) error here.
    patcher.fillInterp(a_advVel, amrMLcrse->m_advVel, amrMLcrse->m_advVel,
                       timeInterpCoeff, 0, 0, 1);
  }
}
void AMRLevelMushyLayer::fillAnalyticVel(FArrayBox &velDir, int dir, int comp,
                                         bool project)
{
  const Box &b = velDir.box();
  for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
  {
    Real x, y, offsetx, offsety;
    IntVect iv = bit();
    RealVect loc, offset;

    IndexType index = b.ixType();
    if (dir == 0)
    {
      offsety = 0.5;
      offsetx = index.cellCentered() ? 0.5 : 0;
    }
    else
    {
      offsetx = 0.5;
      offsety = index.cellCentered() ? 0.5 : 0;
    }

    offset[0] = offsetx;
    offset[1] = offsety;

    ::getLocation(iv, loc, m_dx, offset);

    x = loc[0];
    y = loc[1];

    Real freqx = (m_problem_domain.isPeriodic(0)) ? 2 * M_PI / m_opt.domainWidth
                                                  : M_PI / m_opt.domainWidth;
    //    Real freqy = (m_problem_domain.isPeriodic(1)) ? 2*M_PI/m_domainWidth :
    //    M_PI/m_domainWidth ;

    Real Uscale;
    if (m_parameters.nonDimVel > 0)
    {
      Uscale = m_parameters.nonDimVel * 100;
    }
    else if (m_parameters.darcy > 0)
    {
      Uscale =
          max(m_parameters.rayleighTemp, m_parameters.rayleighComposition) *
          m_parameters.darcy;
    }
    else
    {
      Uscale = sqrt(
          max(m_parameters.rayleighTemp, m_parameters.rayleighComposition));
    }

    Real magU = Uscale; //*0.15;

    if (m_opt.analyticVelType == -1)
    {

      // This velocity should have zero divergence
      if (dir == 0)
      {
        //        velDir(iv, comp) =
        //        magU*sin(freqx*(x+shift))*sin(freqx*(x+shift))*sin(freqx*(y+shift))*cos(freqx*(y+shift));
        velDir(iv, comp) = 0;
      }
      else
      {
        velDir(iv, comp) = magU * sin(freqx * x) * (m_domainHeight - y);
      }
    }
    else if (m_opt.analyticVelType == PhysicalProblems::m_soluteFluxTest)
    {

      if (dir == 0)
      {
        //              y = y + m_dx/2;
        velDir(iv, comp) =
            magU * (cos(freqx * x) + 0.5 * freqx * x) * exp(-freqx * y);
      }
      else
      {
        //              x = x + m_dx/2;
        velDir(iv, comp) = magU * (sin(freqx * x)) * sin(freqx * y);

        //              fab(iv) = -magU*cos(freq*x)*sin(freq*y);

        // velDir(iv) = sin(freq*x);
      }
    }
    else if (m_opt.analyticVelType == PhysicalProblems::m_mushyLayer)
    {
      if (project)
      {
        // This version is divergent, and needs to be projected
        if (dir == 0)
        {
          //              y = y + m_dx/2;
          velDir(iv, comp) =
              -magU * sin(freqx * x) * sin(freqx * x) * cos(freqx * y);
        }
        else
        {
          //              x = x + m_dx/2;
          velDir(iv, comp) =
              magU * cos(freqx * x) * sin(freqx * y) * sin(freqx * y);
        }
      }
      else
      {
        Real shift = 0.5 * m_dx;
        shift = 0;
        // This velocity should have zero divergence
        if (dir == 0)
        {
          //          velDir(iv, comp) = magU*sin(freqx*x)*cos(freqx*y);
          //          velDir(iv, comp) = magU*sin(freqx*x)*cos(freqx*y);

          velDir(iv, comp) =
              magU * sin(freqx * (x + shift)) * sin(freqx * (x + shift)) *
              sin(freqx * (y + shift)) * cos(freqx * (y + shift));
        }
        else
        {

          velDir(iv, comp) =
              -magU * sin(freqx * (x + shift)) * cos(freqx * (x + shift)) *
              sin(freqx * (y + shift)) * sin(freqx * (y + shift));
          //          velDir(iv, comp) = -magU*cos(freqx*x)*sin(freqx*y);
        }
      }
    }
    else
    {
      if (dir == 0)
      {
        velDir(iv, comp) = -magU * sin(freqx * x) * cos(freqx * y);
      }
      else
      {
        velDir(iv, comp) = magU * cos(freqx * x) * sin(freqx * y);
      }
    }
  }
}

void AMRLevelMushyLayer::fillAnalyticVel(LevelData<FArrayBox> &a_advVel)
{
  DataIterator dit = m_grids.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    FArrayBox &vel = a_advVel[dit];
    for (int dir = 0; dir < SpaceDim; dir++)
    {
      fillAnalyticVel(vel, dir, dir, m_opt.projectAnalyticVel);
    }
  }
}

void AMRLevelMushyLayer::fillAnalyticVel(LevelData<FluxBox> &a_advVel)
{
  DataIterator dit = m_grids.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    for (int dir = 0; dir < SpaceDim; dir++)
    {
      FArrayBox &velDir = a_advVel[dit][dir];
      fillAnalyticVel(velDir, dir, 0, m_opt.projectAnalyticVel);
    }
  }
}

void AMRLevelMushyLayer::fillFixedPorosity(LevelData<FArrayBox> &a_porosity)
{
  DataIterator dit = a_porosity.dataIterator();

  Real innerRadius =
      m_opt.fixedPorosityFractionalInnerRadius * m_opt.domainWidth;

  for (dit.reset(); dit.ok(); ++dit)
  {

    BoxIterator bit(a_porosity[dit].box());

    for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      Real x = loc[0];
      Real y = loc[1];

      // Porosity varies linearly from 0.4 at boundary to
      // 1.0 near the middle, c.f. Le Bars and Worster (2006)

      if (m_parameters.m_porosityFunction ==
          ParamsPorosityFunctions::m_porosityLinear)
      {

        Real boundaryPorosity = m_parameters.bcValPorosityLo[0];

        Real xDistFromMiddle = abs(x - m_opt.domainWidth / 2);
        Real yDistFromMiddle = abs(y - m_opt.domainWidth / 2);

        Real gradient =
            (1 - boundaryPorosity) / (0.5 * m_opt.domainWidth - innerRadius);

        Real maxDistFromMiddle = max(xDistFromMiddle, yDistFromMiddle);

        //          Real porosity = maxChi + boundaryPorosity -
        //          2*maxDistFromMiddle;
        Real porosity = 1.0;
        if (maxDistFromMiddle >= innerRadius)
        {
          porosity = 1 - gradient * (maxDistFromMiddle - innerRadius);
        }

        Real maxPorosity = 1.0;
        porosity = min(porosity, maxPorosity);

        // initialDataPoiseuille(x, y, valNew,  m_porosity, -1);
        a_porosity[dit](iv) = porosity;
        //            (*m_scalarOld[ScalarVars::m_porosity])[dit](iv) =
        //            porosity;
      }
      else if (m_parameters.m_porosityFunction ==
               ParamsPorosityFunctions::m_porosityGaussian)
      {

        Real boundaryPorosity = m_parameters.bcValPorosityLo[0];

        Real xd = abs(x - m_opt.domainWidth / 2);
        Real yd = abs(y - m_opt.domainWidth / 2);

        // Real maxDistFromMiddle = max(xDistFromMiddle, yDistFromMiddle);

        // Real porosity = 1 + boundaryPorosity - 2*maxDistFromMiddle;

        Real porosity = boundaryPorosity +
                        m_opt.fixedPorosityMaxChi *
                            exp(-(xd * xd + yd * yd) /
                                (m_opt.FixedPorositySTD * m_opt.domainWidth));

        Real maxPorosity = 1.0;
        porosity = min(porosity, maxPorosity);
        porosity = max(porosity, 0.0);

        // initialDataPoiseuille(x, y, valNew,  m_porosity, -1);
        a_porosity[dit](iv) = porosity;
        //            (*m_scalarOld[ScalarVars::m_porosity])[dit](iv) =
        //            porosity;
      }
      else if (m_parameters.m_porosityFunction ==
               ParamsPorosityFunctions::m_porosityEdge)
      {
        Real boundaryPorosity = m_parameters.bcValPorosityLo[0];

        Real xd = abs(x - m_opt.domainWidth / 2);
        Real yd = abs(y - m_opt.domainWidth / 2);

        // Real maxDistFromMiddle = max(xDistFromMiddle, yDistFromMiddle);

        // Real porosity = 1 + boundaryPorosity - 2*maxDistFromMiddle;

        Real porosity =
            boundaryPorosity *
            (1 - m_opt.fixedPorosityMaxChi *
                     exp(-(xd * xd + yd * yd) /
                         (m_opt.FixedPorositySTD * m_opt.domainWidth)));

        Real maxPorosity = 1.0;
        porosity = min(porosity, maxPorosity);
        porosity = max(porosity, 0.0);

        // initialDataPoiseuille(x, y, valNew,  m_porosity, -1);
        a_porosity[dit](iv) = porosity;
        //            (*m_scalarOld[ScalarVars::m_porosity])[dit](iv) =
        //            porosity;
      }
      else if (m_parameters.m_porosityFunction ==
               ParamsPorosityFunctions::m_porosityConstant)
      {
        // Assume all boundary values are equal, and just pick one
        a_porosity[dit](iv) = m_parameters.bcValPorosityLo[0];
        //            (*m_scalarOld[ScalarVars::m_porosity])[dit](iv) =
        //            m_parameters.bcValPorosityLo[0];
      }
      else if (m_parameters.m_porosityFunction ==
               ParamsPorosityFunctions::m_porosityTimeDependent)
      {
        Real time = m_time;
        if (m_opt.fixedPorosityEndTime >= 0 &&
            m_time >= m_opt.fixedPorosityEndTime)
        {
          time = m_opt.fixedPorosityEndTime;
        }
        //        Real scale = min((1-m_parameters.bcValPorosityHi[1]),
        //        time/m_opt.porosityTimescale);
        Real scale = time / m_opt.porosityTimescale;

        scale = min(scale, 4.0);
        //        Real scale = min(0.9,
        Real chi = 1 - scale * pow(y / m_domainHeight, 2) * 0.3 *
                           sin(M_PI * x / m_opt.domainWidth * 2);
        chi = min(chi, 1.0);
        chi = max(chi, m_opt.lowerPorosityLimit);
        a_porosity[dit](iv) = chi;
      }
    }
  }
}

void AMRLevelMushyLayer::fillTCl(LevelData<FArrayBox> &a_phi, Real a_time,
                                 bool doInterior, bool quadInterp)
{

  bool apply_bcs = false;
  fillMultiComp(a_phi, a_time, m_temperature, m_liquidConcentration, doInterior,
                quadInterp, apply_bcs);

  BCHolder bc = m_physBCPtr->temperatureLiquidSalinityBC();
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
  {
    bc(a_phi[dit], m_grids[dit], m_problem_domain, m_dx, false);
  }
}

void AMRLevelMushyLayer::fillMultiComp(LevelData<FArrayBox> &a_phi, Real a_time,
                                       int scal1, int scal2, bool doInterior,
                                       bool quadInterp, bool apply_bcs)
{
  CH_assert(a_phi.nComp() == 2);

  if (apply_bcs)
  {
    MayDay::Warning("trying to apply bcs in multi comp");
  }

  LevelData<FArrayBox> temp1(a_phi.disjointBoxLayout(), 1, a_phi.ghostVect());
  LevelData<FArrayBox> temp2(a_phi.disjointBoxLayout(), 1, a_phi.ghostVect());
  fillScalars(temp1, a_time, scal1, doInterior, quadInterp);
  fillScalars(temp2, a_time, scal2, doInterior, quadInterp);

  // Need to do copying this way to transfer ghost cells
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
  {
    a_phi[dit].copy(temp1[dit], 0, 0, 1);
    a_phi[dit].copy(temp2[dit], 0, 1, 1);
  }
  //
  //  fillScalars(a_phi, a_time, scal1, doInterior, quadInterp, 0, apply_bcs);
  //  // 0th component fillScalars(a_phi, a_time, scal2, doInterior, quadInterp,
  //  1, apply_bcs); // 1st component
}

void AMRLevelMushyLayer::fillHC(LevelData<FArrayBox> &a_phi, Real a_time,
                                bool doInterior, bool quadInterp)
{
  // Don't apply BCs for HC as they are complicated and we need to do it later
  bool apply_bcs = false;
  fillMultiComp(a_phi, a_time, ScalarVars::m_enthalpy,
                ScalarVars::m_bulkConcentration, doInterior, quadInterp,
                apply_bcs);

  BCHolder bc = m_physBCPtr->enthalpySalinityBC();
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
  {
    bc(a_phi[dit], m_grids[dit], m_problem_domain, m_dx, false);
  }

  doSanityChecks(a_phi, "HC");
}

void AMRLevelMushyLayer::fillScalarFace(LevelData<FluxBox> &a_scal, Real a_time,
                                        const int a_var, bool doInterior,
                                        bool quadInterp)
{
  CH_TIME("AMRLevelMushyLayer::fillScalarFace");
  LOG_DEBUG("FillScalarFace - field: " << m_scalarVarNames[a_var]
                                       << ", level: " << m_level);

  CellToEdgeAveragingMethod method = arithmeticAveraging; // Arithmetic mean
  // Geometric mean for permeability was making pressure projection impossible
  // || a_var == m_pressureScaleVar
  if (a_var == m_porosity ||
      a_var == m_permeability) // || a_var == m_permeability
  {
    method = geometricAveraging; // Geometric mean
    //    method = arithmeticAveraging;
  }

  fillScalarFace(a_scal, a_time, a_var, method, doInterior, quadInterp);
}

void AMRLevelMushyLayer::fillScalarFace(LevelData<FluxBox> &a_scal, Real a_time,
                                        const int a_var,
                                        CellToEdgeAveragingMethod method,
                                        bool doInterior, bool quadInterp,
                                        Real smoothing)
{
  CH_TIME("AMRLevelMushyLayer::fillScalarFace");

  // Need a ghost vector here to get domain faces correct
  LevelData<FArrayBox> temp(m_grids, a_scal.nComp(), IntVect::Unit);
  fillScalars(temp, a_time, a_var, doInterior,
              quadInterp); // NB - this includes exchanges (and corner copiers)

  // To make debugging easier
  DataIterator dit = temp.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    a_scal[dit].setVal(0.0);
  }

  if (smoothing > 0.0)
  {
    smoothScalarField(temp, a_var, smoothing);
  }

  CellToEdge2(temp, a_scal, method);

  // Redo BCs. For now just use extrap. Should eventually do this properly,

  // I think porosity is the only thing we average to faces like this?
  if (a_var == m_porosity)
  {
    EdgeVelBCHolder bc(m_physBCPtr->porosityFaceBC());
    bc.applyBCs(a_scal, a_scal.disjointBoxLayout(), m_problem_domain, m_dx,
                false);
  }
  else if (a_var == m_permeability)
  {
    EdgeVelBCHolder bc(m_physBCPtr->permeabilityFaceBC());
    bc.applyBCs(a_scal, a_scal.disjointBoxLayout(), m_problem_domain, m_dx,
                false);
  }

  // This is important
  a_scal.exchange();

  doRegularisationOps(a_scal, a_var);

  // Try a corner copier?
  if (a_scal.ghostVect()[0] > 0)
  {
    CornerCopier cornerCopy(a_scal.disjointBoxLayout(),
                            a_scal.disjointBoxLayout(), m_problem_domain,
                            a_scal.ghostVect(), true);
    a_scal.exchange(a_scal.interval(), cornerCopy);
  }
}

Real AMRLevelMushyLayer::getCoarseTimeInterpCoeff(Real a_time)
{

  Real crse_new_time = getCoarserLevel()->m_time;
  Real crse_dt = getCoarserLevel()->dt();
  Real crse_old_time = crse_new_time - crse_dt;
  Real crse_time_interp_coeff;

  // check for "essentially 0 or 1"
  if (abs(a_time - crse_old_time) < TIME_EPS)
  {
    crse_time_interp_coeff = 0.0;
  }
  else if (abs(a_time - crse_new_time) < TIME_EPS)
  {
    crse_time_interp_coeff = 1.0;
  }
  else
  {
    crse_time_interp_coeff = (a_time - crse_old_time) / crse_dt;
  }

  if (crse_time_interp_coeff < 0. || crse_time_interp_coeff > 1.)
  {
    LOG_INFO("ERROR crse_time_interp_coeff = " << crse_time_interp_coeff);
    LOG_INFO("Computing time interp coeff on level " << m_level);
    LOG_INFO("Level " << m_level << ":     a_time=" << a_time
                      << ",   m_time=" << m_time << ",m_dt=" << m_dt);
    LOG_INFO("Level " << getCoarserLevel()->m_level
                      << ": old_time=" << crse_old_time
                      << ", new_time=" << crse_new_time << ",m_dt=" << crse_dt);
  }

  CH_assert(crse_time_interp_coeff >= 0.);
  CH_assert(crse_time_interp_coeff <= 1.);
  return crse_time_interp_coeff;
}
// Fill a single component of a scalar field
void AMRLevelMushyLayer::fillScalars(LevelData<FArrayBox> &a_scal, Real a_time,
                                     const int a_var, bool doInterior,
                                     bool quadInterp, int a_comp,
                                     bool apply_bcs)
{

  CH_TIME("AMRLevelMushyLayer::fillScalars");
  LOG_DEBUG_EXTRA("Fill scalar field: " << m_scalarVarNames[a_var]
                                        << ", level: " << m_level);

  Interval scalComps = Interval(a_comp, a_comp);
  Interval srcComps = Interval(0, 0);

  Real old_time = m_time - m_dt;

  if (doInterior)
  {
    CH_TIME("AMRLevelMushyLayer::fillScalars::interior");

    if (abs(a_time - old_time) < TIME_EPS)
    {
      m_scalarOld[a_var]->copyTo(srcComps, a_scal, scalComps);
    }
    else if (abs(a_time - m_time) < TIME_EPS)
    {
      m_scalarNew[a_var]->copyTo(srcComps, a_scal, scalComps);
    }
    else
    {
      // do linear interpolation in time
      timeInterp(a_scal, a_time, *m_scalarOld[a_var], old_time,
                 *m_scalarNew[a_var], m_time, srcComps, scalComps);
    }
  }

  // Only do CF interp if we have ghost vectors
  if (m_level > 0 && a_scal.ghostVect() >= IntVect::Unit)
  {
    CH_TIME("AMRLevelMushyLayer::fillScalars::coarseFineInterp");

    // Make sure fill objects are defined
    if (!m_piecewiseLinearFillPatchScalarOne.isDefined())
    {
      defineCFInterp();
    }

    // fill in coarse-fine BC data by conservative linear interp
    AMRLevelMushyLayer &crseLevel = *getCoarserLevel();

    LevelData<FArrayBox> &oldCrseScal = *(crseLevel.m_scalarOld[a_var]);
    LevelData<FArrayBox> &newCrseScal = *(crseLevel.m_scalarNew[a_var]);

    const DisjointBoxLayout &crseGrids = oldCrseScal.getBoxes();
    //    const ProblemDomain& crseDomain = crseLevel.problemDomain();
    //    int nRefCrse = crseLevel.refRatio();

    Real crse_time_interp_coeff = getCoarseTimeInterpCoeff(a_time);

    const IntVect &scalGrowVect = a_scal.ghostVect();
    int scalGrow = scalGrowVect[0];

    Interval subcomps(a_comp, a_comp);
    LevelData<FArrayBox> alias;
    aliasLevelData(alias, &a_scal, subcomps);

    {
      CH_TIME("AMRLevelMushyLayer::fillScalars::linearFillPatch");

      // Now we're using an alias:
      int srcComp = 0; // the data to interpolate from always has 1 comp
      int destComp = 0;
      int numComp = 1;

      if (scalGrow == 1)
      {
        m_piecewiseLinearFillPatchScalarOne.fillInterp(
            alias, oldCrseScal, newCrseScal, crse_time_interp_coeff, srcComp,
            destComp, numComp);
      }
      else if (scalGrow == 2)
      {
        m_piecewiseLinearFillPatchScalarTwo.fillInterp(
            alias, oldCrseScal, newCrseScal, crse_time_interp_coeff, srcComp,
            destComp, numComp);
      }
      else if (scalGrow == 3)
      {
        m_piecewiseLinearFillPatchScalarThree.fillInterp(
            alias, oldCrseScal, newCrseScal, crse_time_interp_coeff, srcComp,
            destComp, numComp);
      }
      else if (scalGrow == 4)
      {
        m_piecewiseLinearFillPatchScalarFour.fillInterp(
            alias, oldCrseScal, newCrseScal, crse_time_interp_coeff, srcComp,
            destComp, numComp);
      }
      else
      {
        LOG_INFO("ERROR: No fill patch object for " << scalGrow
                                                    << "ghost cells");
        MayDay::Error("No fill patch object for this number of ghost cells");
      }
    }

    if (quadInterp && m_quadCFInterpScalar.isDefined())
    {
      CH_TIME("AMRLevelMushyLayer::fillScalars::quadCFinterp");

      LevelData<FArrayBox> avCrseScal(crseGrids, 1, oldCrseScal.ghostVect());
      Real crse_new_time = getCoarserLevel()->m_time;
      Real crse_dt = getCoarserLevel()->dt();
      Real crse_old_time = crse_new_time - crse_dt;
      ::timeInterp(avCrseScal, a_time, oldCrseScal, crse_old_time, newCrseScal,
                   crse_new_time, Interval(0, 0));

      m_quadCFInterpScalar.coarseFineInterp(alias, avCrseScal);
    }
  }

  const ProblemDomain &physDomain = problemDomain();

  int numGhost = a_scal.ghostVect()[0];

  // Do domain BCs if we have ghost cells
  if (numGhost > 0 && apply_bcs)
  {
    CH_TIME("AMRLevelMushyLayer::fillScalars::domainBCs");

    BCHolder thisBC;
    getScalarBCs(thisBC, a_var, false); // inhomogeneous

    // loop over boxes
    DataIterator dit = a_scal.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      thisBC(a_scal[dit], m_grids[dit], physDomain, m_dx,
             false); // inhomogeneous
    }
  }

  doRegularisationOps(a_scal, a_var, a_comp);

  {
    CH_TIME("AMRLevelMushyLayer::fillScalars::exchange");

    a_scal.exchange();

    if (a_scal.ghostVect()[0] > 0 && m_opt.scalarExchangeCorners)
    {
      CornerCopier cornerCopy(m_grids, m_grids, m_problem_domain,
                              a_scal.ghostVect(), true);
      a_scal.exchange(a_scal.interval(), cornerCopy);
    }
  }

  LOG_FUNCTION_EXIT();
}

void AMRLevelMushyLayer::fillUnprojectedDarcyVelocity(
    LevelData<FluxBox> &a_advVel, Real time)
{
  IntVect ghost = a_advVel.ghostVect();
  DataIterator dit = a_advVel.dataIterator();
  LevelData<FluxBox> permeability_face(m_grids, 1, ghost);
  LevelData<FluxBox> T_face(m_grids, 1, ghost);
  LevelData<FluxBox> C_face(m_grids, 1, ghost);
  LevelData<FluxBox> viscosity_face(m_grids, 1, ghost);
  LevelData<FluxBox> viscous_velocity_face(m_grids, 1, ghost);
  fillScalarFace(permeability_face, time, ScalarVars::m_permeability, true,
                 false);
  fillScalarFace(T_face, time, ScalarVars::m_temperature, true, false);
  fillScalarFace(C_face, time, ScalarVars::m_liquidConcentration, true, false);
  fillScalarFace(viscosity_face, time, ScalarVars::m_viscosity, true, false);

  if (m_opt.explicitViscousVelSolve)
  {
    LevelData<FArrayBox> viscous(m_grids, SpaceDim, ghost);
    LevelData<FArrayBox> *crseVelPtr = nullptr;

    computeLapVel(viscous, *m_vectorOld[m_fluidVel], crseVelPtr);

    CellToEdge(viscous, viscous_velocity_face);
  }

  //  int gravityDir = SpaceDim-1;

  for (dit.reset(); dit.ok(); ++dit)
  {
    a_advVel[dit].setVal(0.0);

    FArrayBox &fabVelz = a_advVel[dit][SpaceDim - 1];

    fabVelz.plus(T_face[dit][SpaceDim - 1], m_parameters.m_buoyancyTCoeff);
    fabVelz.plus(C_face[dit][SpaceDim - 1], -m_parameters.m_buoyancySCoeff);
    fabVelz.plus(m_parameters.body_force);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_advVel[dit][idir].mult(permeability_face[dit][idir]);

      if (m_parameters.m_viscosityFunction !=
          ViscosityFunction::uniformViscosity)
      {
        a_advVel[dit][idir].divide(viscosity_face[dit][idir]);
      }

      if (m_opt.explicitViscousVelSolve)
      {
        a_advVel[dit][idir].plus(viscous_velocity_face[dit][idir]);
      }
    }
  }

  EdgeToCell(a_advVel, *m_vectorNew[VectorVars::m_advUstar]);
  Real maxVel = ::computeNorm(*m_vectorNew[VectorVars::m_advUstar], nullptr, 1,
                              m_dx, Interval(0, SpaceDim - 1), 0);
  if (maxVel > 1e290)
  {
    LOG_INFO("WARNING: max(unprojected advection velocity) = " << maxVel);
  }
  else
  {
    LOG_INFO("Max unprojected vel=" << maxVel);
  }
}

void AMRLevelMushyLayer::fillPressureSrcTerm(
    LevelData<FArrayBox> &gradP, LevelData<FArrayBox> &pressureScale,
    Real a_time, bool a_MACprojection)
{

  computeGradP(gradP, a_time, a_MACprojection);

  DataIterator dit = gradP.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    if ((a_MACprojection && m_opt.scaleP_MAC) ||
        (!a_MACprojection && m_opt.scaleP_CC))
    {
      for (int dir = 0; dir < SpaceDim; dir++)
      {
        gradP[dit].mult(pressureScale[dit], 0, dir);
      }
    }
  }
}

void AMRLevelMushyLayer::fillBuoyancy(LevelData<FArrayBox> &a_buoyancy,
                                      Real a_time)
{
  IntVect ghost = a_buoyancy.ghostVect();

  if (m_level == 0)
  {
    DisjointBoxLayout dbl = a_buoyancy.disjointBoxLayout();

    LevelData<FArrayBox> temperature(dbl, 1, ghost);
    LevelData<FArrayBox> liquidConc(dbl, 1, ghost);
    LevelData<FArrayBox> porosity(dbl, 1, ghost);

    fillScalars(temperature, a_time, m_temperature, true, true);
    fillScalars(liquidConc, a_time, m_liquidConcentration, true, true);
    fillScalars(porosity, a_time, m_porosity, true, true);

    fillBuoyancy(a_buoyancy, temperature, liquidConc, porosity);
  }
  else
  {
    AMRLevelMushyLayer *mlCrse = getCoarserLevel();
    LevelData<FArrayBox> crseBuoyancy(mlCrse->m_grids, SpaceDim, ghost);
    mlCrse->fillBuoyancy(crseBuoyancy, a_time);

    // Refine back to this level
    m_fineInterpVector.interpToFine(a_buoyancy, crseBuoyancy);

    // Fill ghost cells
    PiecewiseLinearFillPatch filpatcher(
        m_grids, mlCrse->m_grids, a_buoyancy.nComp(), mlCrse->m_problem_domain,
        mlCrse->m_ref_ratio, ghost[0], false);

    Interval scalComps = a_buoyancy.interval();

    filpatcher.fillInterp(a_buoyancy, crseBuoyancy, crseBuoyancy, 0, 0, 0,
                          scalComps.size());
  }
}

void AMRLevelMushyLayer::fillBuoyancy(LevelData<FArrayBox> &a_buoyancy,
                                      LevelData<FArrayBox> &a_temperature,
                                      LevelData<FArrayBox> &a_liquidConc,
                                      LevelData<FArrayBox> &a_porosity)
{

  for (DataIterator dit = a_buoyancy.dataIterator(); dit.ok(); ++dit)
  {
    if (m_opt.buoyancy_zero_time >= 0 && m_time > m_opt.buoyancy_zero_time)
    {
      a_buoyancy[dit].setVal(0);
    }
    else
    {
      fillBuoyancy(a_buoyancy[dit], a_temperature[dit], a_liquidConc[dit],
                   a_porosity[dit],
                   (*m_vectorNew[VectorVars::m_bodyForce])[dit]);
    }
  }
}

void AMRLevelMushyLayer::fillBuoyancy(FArrayBox &buoyancy,
                                      FArrayBox &temperature,
                                      FArrayBox &liquidConc,
                                      FArrayBox &porosity,
                                      FArrayBox &bodyForce)
{
  int gravityDir = SpaceDim - 1;

  // buoyancy = porosity * (Ra_T * theta - Ra_c*Theta_l) e_z
  buoyancy.setVal(0.0);
  buoyancy.plus(temperature, m_parameters.m_buoyancyTCoeff, 0, gravityDir, 1);
  buoyancy.plus(liquidConc, -m_parameters.m_buoyancySCoeff, 0, gravityDir, 1);
  buoyancy.plus(bodyForce, 1, gravityDir, gravityDir, 1);

  buoyancy.mult(porosity, buoyancy.box(), 0, gravityDir, 1);
}

void AMRLevelMushyLayer::ensureGreaterThan(LevelData<FArrayBox> &a_phi,
                                           Real min_val, string field_name,
                                           int comp)
{
  LOG_INFO("Check " << field_name << " > " << min_val);
  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    // Would probably be quicker to write this in Fortran,
    // but as we only do this after regridding it's probably not worth it
    FArrayBox &phiFAB = a_phi[dit];
    for (BoxIterator bit = BoxIterator(phiFAB.box()); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      Real val = phiFAB.get(iv, comp);
      if (val < min_val)
      {
        LOG_WARNING(field_name << " (" << iv << ")=" << val << " which is < "
                               << min_val);
        phiFAB.set(iv, 0, min_val);
      }
    }
  }
}
