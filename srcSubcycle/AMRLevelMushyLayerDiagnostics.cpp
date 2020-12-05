#include "AMRLevelMushyLayer.H"



void AMRLevelMushyLayer::computeDiagnostics()
{
  if (!m_opt.computeDiagnostics)
  {
    return;
  }

  // Only compute diagnostics on level 0
  if (m_level > 0)
  {
    return;
  }

  // If a diagnostic period has been declared, check this time has passed since we last produced diagnostics
  if (m_opt.diagnostics_period > 0 &&
      m_time - m_prev_diag_output < m_opt.diagnostics_period)
  {
    return;
  }

  CH_TIME("AMRLevelMushyLayer::computeDiagnostics");

  bool calcDiagnostics = (m_level == 0); // do this on all processors!

  Vector<AMRLevelMushyLayer*> mlVect;
  Vector<DisjointBoxLayout> grids;
  Vector<int> refRat;
  ProblemDomain lev0Dom;
  Real lev0Dx;
  getHierarchyAndGrids(mlVect, grids, refRat, lev0Dom, lev0Dx);
  int nLevels = mlVect.size();

  // Compute Nusselt number
  {
    CH_TIME("AMRLevelMushyLayer::computeDiagnostics::computeNusselt");

    if ((m_parameters.physicalProblem == PhysicalProblems::m_HRL
        || m_parameters.physicalProblem == PhysicalProblems::m_rayleighBenard)
        && m_level == 0)
    {
      fillScalars(*m_scalarNew[ScalarVars::m_temperature], m_time, m_temperature);

      Real Nu = ::computeNusselt(*m_scalarNew[ScalarVars::m_temperature], *m_vectorNew[VectorVars::m_fluidVel],
                                 m_dx, m_parameters,
                                 m_opt.domainWidth, m_domainHeight);

      m_diagnostics.addDiagnostic(DiagnosticNames::diag_Nu, m_time, Nu);
    }
    else if (m_parameters.physicalProblem == PhysicalProblems::m_convectionMixedPorous
        && m_level == 0)
    {

      // Calculate average dT/dx at left and right boundaries

      //    Vector<LevelData<FArrayBox>* > T(nLevels);
      Vector<LevelData<FArrayBox>* > gradT(nLevels);
      Vector<LevelData<FArrayBox>* > AMRNu(nLevels);
      Vector<LevelData<FluxBox>* > gradEdgeT(nLevels);
      LevelData<FArrayBox>* crseT = nullptr;
      //    LevelData<FArrayBox>* fineT = nullptr;

      // Use domain flux registers to calculate Nu
      Vector<LevelDomainFluxRegister*> Tflux(nLevels);

      LevelDomainFluxRegister* fineFR = nullptr;
      LevelDomainFluxRegister* coarseFR = nullptr;

      for (int lev=0; lev<nLevels; lev++)
      {
        Tflux[lev] = new LevelDomainFluxRegister();
      }

      Real scale = 1.0;

      for (int lev=0; lev<nLevels; lev++)
      {
        AMRLevelMushyLayer* ml = mlVect[lev];

        gradT[lev] = new LevelData<FArrayBox>(ml->m_grids, SpaceDim);
        gradEdgeT[lev] = new LevelData<FluxBox>(ml->m_grids, SpaceDim);
        AMRNu[lev] = new LevelData<FArrayBox>(ml->m_grids, 1);

        if (lev > 0)
        {
          coarseFR = Tflux[lev-1];
        }
        else
        {
          coarseFR = nullptr;
        }

        if (lev < (nLevels-1))
        {
          fineFR = Tflux[lev+1];
        }
        else
        {
          fineFR = nullptr;
        }


        Tflux[lev]->define(ml->m_problem_domain, ml->m_grids,refRat[lev], ml->m_dx, fineFR, coarseFR, SpaceDim);

        AMRLevelMushyLayer* mlCrse = ml->getCoarserLevel();
        if (mlCrse)
        {
          crseT = mlCrse->m_scalarNew[ScalarVars::m_temperature];
        }

        Box domainBox = ml->m_problem_domain.domainBox();

        // Ensure BCs are set
        ml->fillScalars(*ml->m_scalarNew[ScalarVars::m_temperature], ml->m_time, m_temperature, false, true);

        Gradient::levelGradientMAC(*gradEdgeT[lev], (*ml->m_scalarNew[ScalarVars::m_temperature]),
                                   ml->m_dx);

        scale = 1.0/ml->m_dx;
        Tflux[lev]->incrFlux(*gradEdgeT[lev], scale);

      } // end loop over levels


      // Scale with coarsest dx
      scale = -1.0;

      // This makes sure we get the most accurate fluxes (using the finest level available)
      Real Nuleft  = Tflux[0]->getFluxHierarchy(0, Side::Lo, scale);
      Real Nuright = Tflux[0]->getFluxHierarchy(0, Side::Hi, scale);

      Real Nu = 0.5*abs(Nuleft + Nuright);

      m_diagnostics.addDiagnostic(DiagnosticNames::diag_Nu, m_time, Nu);
      m_diagnostics.addDiagnostic(DiagnosticNames::diag_NuLeft, m_time, Nuleft);
      m_diagnostics.addDiagnostic(DiagnosticNames::diag_NuRight, m_time, Nuright);

      // Also compute the flux in the middle of the domain on level 0

      LevelData<FluxBox> fluidFlux(m_grids, 2, IntVect::Unit);
      computeTotalAdvectiveFluxes(fluidFlux);

      Real Nu_sum = 0;
      Box middleBox(m_problem_domain.domainBox());
      middleBox.surroundingNodes(0);

      int width = middleBox.hiVect()[0]-middleBox.loVect()[0];
      int shrink_cells = int(width/2);
      middleBox.growLo(0, -shrink_cells);
      middleBox.growHi(0, -shrink_cells);
//      middleBox.shift(0, 0.5);

      // we already have dt/dx, now add u*T to it
      int numBoxes = 0;
      for (DataIterator dit = gradEdgeT[0]->dataIterator(); dit.ok(); ++dit)
      {
        FluxBox& fb = (*gradEdgeT[0])[dit];
        FArrayBox& horizFlux = fb[0];


        Box b(horizFlux.box());
        b &= middleBox;
        int numBoxCells = b.numPts();

        Real thisNu = 0;

        if (numBoxCells > 0)
        {
          FArrayBox& fluidFluxFAB = fluidFlux[dit][0];
          horizFlux.minus(fluidFluxFAB);
          thisNu =  horizFlux.sum(b, 0);
        }
        else
        {
          // Skip to next box
          continue;
        }

        thisNu =  thisNu/numBoxCells;

        Nu_sum = Nu_sum + thisNu;
        numBoxes++;
      }

      Nu_sum = abs(Nu_sum) / numBoxes;
      m_diagnostics.addDiagnostic(DiagnosticNames::diag_NuMiddle, m_time, Nu_sum);

      // Cleanup after nusselt calculation
      for (int lev=0; lev<nLevels; lev++)
      {
        if (Tflux[lev] != nullptr)
        {
          delete Tflux[lev];
          Tflux[lev] = nullptr;
        }

        if (gradT[lev] != nullptr)
        {
          delete gradT[lev];
          gradT[lev] = nullptr;
        }

        if ( gradEdgeT[lev] != nullptr)
        {
          delete gradEdgeT[lev] ;
          gradEdgeT[lev]  = nullptr;
        }

        if ( AMRNu[lev]  != nullptr)
        {
          delete  AMRNu[lev] ;
          AMRNu[lev]  = nullptr;
        }

      }

    }

    if (m_diagnostics.diagnosticIsIncluded(DiagnosticNames::diag_maxUhalf))
    {
      CH_TIME("AMRLevelMushyLayer::computeDiagnostics::computeMaxVel");

      // Compute max velocities at midpoints
      // Have to do this ourselves, rather than using Chombo's computeMax(), as we
      // want to compute max over a single strip of cells

      Box horizBox = ::adjCellHi(m_problem_domain.domainBox(), 1, 1);
      horizBox.shift(1, -int(m_numCells[1]/2));

      Real maxVertVel = 0.0;

      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
      {
        Box b = m_grids[dit];
        b &= horizBox;

        if (b.size() >= IntVect::Unit)
        {
          FArrayBox& vel = (*m_vectorNew[VectorVars::m_fluidVel])[dit];
          maxVertVel = max(maxVertVel, vel.max(b, 1));
        }
      }

      // Do broadcast/gather on this
      int srcProc = 0;
      Vector<Real> allMax(numProc(), 0.0);
      gather(allMax, maxVertVel, srcProc);
      Real globalMaxVertVel = 0;
      if (procID() == srcProc)
      {
        for (int ivec = 0; ivec<numProc(); ivec++)
        {
          globalMaxVertVel = max(globalMaxVertVel, allMax[ivec]);
        }
      }
      broadcast(globalMaxVertVel, srcProc);


      Box vertBox = ::adjCellHi(m_problem_domain.domainBox(), 0, 1);
      vertBox.shift(0, -int(m_numCells[0]/2));

      Real maxHorizVel = 0.0;
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
      {
        Box b = m_grids[dit];
        b &= vertBox;

        if (b.size() >= IntVect::Unit)
        {
          FArrayBox& vel = (*m_vectorNew[VectorVars::m_fluidVel])[dit];
          maxHorizVel = max(maxHorizVel, vel.max(b, 0));
        }
      }

      gather(allMax, maxHorizVel, srcProc);
      Real globalMaxHorizVel = 0;
      if (procID() == srcProc)
      {
        for (int ivec = 0; ivec<numProc(); ivec++)
        {
          globalMaxHorizVel = max(globalMaxHorizVel, allMax[ivec]);
        }
      }
      broadcast(globalMaxHorizVel, srcProc);

      m_diagnostics.addDiagnostic(DiagnosticNames::diag_maxVhalf, m_time, globalMaxHorizVel);
      m_diagnostics.addDiagnostic(DiagnosticNames::diag_maxUhalf, m_time, globalMaxVertVel);

    }
  }


  if (m_diagnostics.diagnosticIsIncluded(DiagnosticNames::diag_soluteFluxTop))
  {
    CH_TIME("AMRLevelMushyLayer::computeDiagnostics::computeSoluteFluxes");

    // Need to calculate:
    // 1. Sum of new bulk salinity over valid regions
    // 2. Sum of old bulk salinity over valid regions
    // 3. Salt flux at top of domain over valid regions
    // 4. Salt flux at bottom of domain over valid regions

    Vector<LevelData<FArrayBox>* > S_new, S_old, dSdt, H_new, horizAvFs, Fs_vert_diffusion, Fs_vert_fluid, Fs_vert_frame;
    Vector<LevelData<FluxBox>* > soluteFluxTop, soluteFluxBottom;

    S_new.resize(nLevels);
    H_new.resize(nLevels);
    horizAvFs.resize(nLevels);
    Fs_vert_diffusion.resize(nLevels);
    Fs_vert_fluid.resize(nLevels);
    Fs_vert_frame.resize(nLevels);

    // Remember we are on level 0 here
    AMRLevelMushyLayer* ml = this;

    // Compute fluxes and setup diagnostics on all levels
    for (int lev = 0; lev < nLevels; lev++)
    {
      ml->computeDiagnosticSoluteFluxes();
      ml->computeTotalFlux();

      S_new[lev] = ml->m_scalarNew[ScalarVars::m_bulkConcentration];
      H_new[lev] = ml->m_scalarNew[ScalarVars::m_enthalpy];
      horizAvFs[lev] = new LevelData<FArrayBox>(ml->m_grids, 1, IntVect::Zero); // no ghost vector
      Fs_vert_diffusion[lev] = new LevelData<FArrayBox>(ml->m_grids, 1, IntVect::Zero); // no ghost vector
      Fs_vert_fluid[lev] = new LevelData<FArrayBox>(ml->m_grids, 1, IntVect::Zero); // no ghost vector
      Fs_vert_frame[lev] = new LevelData<FArrayBox>(ml->m_grids, 1, IntVect::Zero); // no ghost vector

      ml->m_scalarNew[ScalarVars::m_averageVerticalFlux]->copyTo(*horizAvFs[lev]);
      ml->m_scalarNew[ScalarVars::m_FsVertDiffusion]->copyTo(*Fs_vert_diffusion[lev]);
      ml->m_scalarNew[ScalarVars::m_FsVertFluid]->copyTo(*Fs_vert_fluid[lev]);
      ml->m_scalarNew[ScalarVars::m_FsVertFrame]->copyTo(*Fs_vert_frame[lev]);

      ml = ml->getFinerLevel();
    }

    // Get horizontally averaged salt flux at different levels in the domain
    // 10%, 20% and 30% above the bottom

    Vector<Real> averageVertFs;
    horizontallyAverage(averageVertFs, *m_scalarNew[ScalarVars::m_averageVerticalFlux]);
    int j_top = m_problem_domain.domainBox().bigEnd()[SpaceDim-1];
    int j_bottom = m_problem_domain.domainBox().smallEnd()[SpaceDim-1];
    int domHeight = j_top-j_bottom;
    int j_10 = j_bottom + int(0.1*domHeight);
    int j_20 = j_bottom + int(0.2*domHeight);
    int j_30 = j_bottom + int(0.3*domHeight);
    int j_40 = j_bottom + int(0.4*domHeight);
    int j_50 = j_bottom + int(0.5*domHeight);

    m_diagnostics.addDiagnostic(DiagnosticNames::diag_Fs10, m_time, averageVertFs[j_10]);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_Fs20, m_time, averageVertFs[j_20]);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_Fs30, m_time, averageVertFs[j_30]);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_Fs40, m_time, averageVertFs[j_40]);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_Fs50, m_time, averageVertFs[j_50]);

    // Another diagnostic - average vertical solute flux across the whole domain
    Real volume = 0.0;
    Real domainAverageFs = computeSum(volume, horizAvFs, refRat, lev0Dx, Interval(0,0), 0);
    domainAverageFs = domainAverageFs/volume;
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_averageVerticalSaltFlux, m_time, domainAverageFs);

    Real L2FsVertDiffusion = computeNorm(Fs_vert_diffusion, refRat, lev0Dx, Interval(0,0), 2);
    Real L2FsVertFluid = computeNorm(Fs_vert_fluid, refRat, lev0Dx, Interval(0,0), 2);
    Real L2FsVertFrame = computeNorm(Fs_vert_frame, refRat, lev0Dx, Interval(0,0), 2);

    m_diagnostics.addDiagnostic(DiagnosticNames::diag_L2FsVertDiffusion, m_time, L2FsVertDiffusion);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_L2FsVertFluid, m_time, L2FsVertFluid);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_L2FsVertFrame, m_time, L2FsVertFrame);


    Real L1FsVertDiffusion = computeNorm(Fs_vert_diffusion, refRat, lev0Dx, Interval(0,0), 1);
    Real L1FsVertFluid = computeNorm(Fs_vert_fluid, refRat, lev0Dx, Interval(0,0), 1);
    Real L1FsVertFrame = computeNorm(Fs_vert_frame, refRat, lev0Dx, Interval(0,0), 1);

    m_diagnostics.addDiagnostic(DiagnosticNames::diag_L1FsVertDiffusion, m_time, L1FsVertDiffusion);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_L1FsVertFluid, m_time, L1FsVertFluid);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_L1FsVertFrame, m_time, L1FsVertFrame);

    Real L0FsVertDiffusion = computeNorm(Fs_vert_diffusion, refRat, lev0Dx, Interval(0,0), 0);
    Real L0FsVertFluid = computeNorm(Fs_vert_fluid, refRat, lev0Dx, Interval(0,0), 0);
    Real L0FsVertFrame = computeNorm(Fs_vert_frame, refRat, lev0Dx, Interval(0,0), 0);

    m_diagnostics.addDiagnostic(DiagnosticNames::diag_L0FsVertDiffusion, m_time, L0FsVertDiffusion);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_L0FsVertFluid, m_time, L0FsVertFluid);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_L0FsVertFrame, m_time, L0FsVertFrame);

    //Clean up
    for(int lev = 0; lev < nLevels; lev++)
    {
      delete horizAvFs[lev];
      horizAvFs[lev] = nullptr;

      delete Fs_vert_diffusion[lev];
      Fs_vert_diffusion[lev] = nullptr;

      delete Fs_vert_fluid[lev];
      Fs_vert_fluid[lev]= nullptr;

      delete Fs_vert_frame[lev];
      Fs_vert_frame[lev]= nullptr;
    }

    // Check for flux conservation

    // Calculate sums over valid regions

    // Volume averaged sum i.e. sum[H*dx^(SpaceDim)]
    AMRSaltSum_new = ::computeSum(S_new, refRat, lev0Dx, Interval(0,0), 0);
    AMREnthalpySum_new = ::computeSum(H_new, refRat, lev0Dx, Interval(0,0), 0);

    Real Fs_bottom = m_saltDomainFluxRegister.getFluxHierarchy(1, Side::Lo, 1.0);
    Real Fs_top = m_saltDomainFluxRegister.getFluxHierarchy(1, Side::Hi, 1.0);

    Real Fh_bottom = m_heatDomainFluxRegister.getFluxHierarchy(1, Side::Lo, 1.0);
    Real Fh_top =m_heatDomainFluxRegister.getFluxHierarchy(1, Side::Hi, 1.0);

    Real scale = 1/(m_dt*m_opt.domainWidth);
    Real totalF_bottom = Fs_bottom*scale;
    Real totalF_top = Fs_top*scale;

    Real totalFh_top = Fh_top*scale;
    Real totalFh_bottom = Fh_bottom*scale;

    Real ds_diff = (AMRSaltSum_new - AMRSaltSum_old);
    Real H_diff = (AMREnthalpySum_new - AMREnthalpySum_old);

    // Calculate flux differences by considering all sides/directions
    Real salt_flux_diff = 0;
    Real heat_flux_diff = 0;
    for (int dir=0; dir < SpaceDim; dir++)
    {
      Real Fs_hi = m_saltDomainFluxRegister.getFluxHierarchy(dir, Side::Hi, 1.0);
      Real Fs_lo = m_saltDomainFluxRegister.getFluxHierarchy(dir, Side::Lo, 1.0);

      Real Fh_hi = m_heatDomainFluxRegister.getFluxHierarchy(dir, Side::Hi, 1.0);
      Real Fh_lo = m_heatDomainFluxRegister.getFluxHierarchy(dir, Side::Lo, 1.0);

      salt_flux_diff += Fs_hi - Fs_lo;
      heat_flux_diff += Fh_hi - Fh_lo;
    }

    Real salt_mismatch = (ds_diff + salt_flux_diff);
    Real heat_mismatch = (H_diff + heat_flux_diff);

    Real rel_salt_mismatch = salt_mismatch/salt_flux_diff;
    Real rel_heat_mismatch = heat_mismatch/heat_flux_diff;

    m_diagnostics.addDiagnostic(DiagnosticNames::diag_soluteFluxBottom, m_time, totalF_bottom);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_soluteFluxTop, m_time, totalF_top);

    m_diagnostics.addDiagnostic(DiagnosticNames::diag_heatFluxAbsMismatch, m_time, heat_mismatch);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_saltFluxAbsMismatch, m_time, salt_mismatch);

    m_diagnostics.addDiagnostic(DiagnosticNames::diag_heatFluxRelMismatch, m_time, rel_heat_mismatch);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_saltFluxRelMismatch, m_time, rel_salt_mismatch);

    m_diagnostics.addDiagnostic(DiagnosticNames::diag_heatFluxBottom, m_time, totalFh_bottom);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_heatFluxTop, m_time, totalFh_top);

    // Clean up ready for next full timestep
    AMRSaltSum_old = AMRSaltSum_new;
    AMREnthalpySum_old = AMREnthalpySum_new;
    ml = this;

    while(ml)
    {
      setValLevel(ml->m_saltFluxTop, 0.0);
      setValLevel(ml->m_saltFluxBottom, 0.0);

      (ml->m_heatDomainFluxRegister).setToZero();
      (ml->m_saltDomainFluxRegister).setToZero();

      ml = ml->getFinerLevel();

    }

  } // end if we care about solute fluxes


  if (calcDiagnostics
      && m_diagnostics.diagnosticIsIncluded(DiagnosticNames::diag_maxLambda))
  {
    CH_TIME("AMRLevelMushyLayer::computeDiagnostics::computeMaxLambda");

    Real maxVel = getMaxVelocity();
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_maxVel, m_time, maxVel);


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

    ml = this->getCoarsestLevel();
    for (int lev = ml->m_level; lev < amrLambda.size(); lev++)
    {
      LevelData<FArrayBox>* levelLambda = amrLambda[lev];
      for (DataIterator dit = levelLambda->dataIterator(); dit.ok(); ++dit)
      {
        (*levelLambda)[dit].plus(-1.0);
      }
    }

    Real sumLambda = ::computeNorm(amrLambda,nRef, m_dx, Interval(0,0), 2, m_level) ;
    Real maxLambda = ::computeNorm(amrLambda, nRef, m_dx, Interval(0,0), 0, m_level); //::computeSum(amrLambda,nRef,m_dx, Interval(0,0), m_level);

    Real maxFreestreamCorr = ::computeNorm(*m_vectorNew[VectorVars::m_freestreamCorrection], nullptr, 1, m_dx, Interval(0, SpaceDim-1), 0);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_maxFreestreamCorrection, m_time, maxFreestreamCorr);

    m_diagnostics.addDiagnostic(DiagnosticNames::diag_maxLambda, m_time, maxLambda);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_sumLambda, m_time, sumLambda);

    ml = this->getCoarsestLevel();
    for (int lev = ml->m_level; lev < amrLambda.size(); lev++)
    {
      LevelData<FArrayBox>* levelLambda = amrLambda[lev];
      for (DataIterator dit = levelLambda->dataIterator(); dit.ok(); ++dit)
      {
        (*levelLambda)[dit].plus(1);
      }
    }
  }

  // Work out average liquid salinity at various vertical points in domain
  if (calcDiagnostics
      && m_diagnostics.diagnosticIsIncluded(DiagnosticNames::diag_avSalinity))
  {
    Vector<Real> averagedSalinity;
    horizontallyAverage(averagedSalinity, *m_scalarNew[ScalarVars::m_liquidConcentration]);

    // Also do average of liquid salinity over liquid regions
    Real averageSl = averageOverLiquidRegion(ScalarVars::m_liquidConcentration);


    int y_0 = 0.0;
    int y_1 = round(averagedSalinity.size()*0.2*(2-1));
    int y_2 = round(averagedSalinity.size()*0.2*(3-1));
    int y_3 = round(averagedSalinity.size()*0.2*(4-1));

    m_diagnostics.addDiagnostic(DiagnosticNames::diag_HorizAvSalinity0, m_time, averagedSalinity[y_0]);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_HorizAvSalinity20, m_time, averagedSalinity[y_1]);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_HorizAvSalinity40, m_time, averagedSalinity[y_2]);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_HorizAvSalinity60, m_time, averagedSalinity[y_3]);

    m_diagnostics.addDiagnostic(DiagnosticNames::diag_avSalinity, m_time, averageSl);


  }

  // Work out mushy layer depth
  Real depth = computeMushDepth();

  if (calcDiagnostics
          && m_diagnostics.diagnosticIsIncluded(DiagnosticNames::diag_mushDepth))
  {
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_mushDepth, m_time, depth);
  }

  if (calcDiagnostics
        && m_diagnostics.diagnosticIsIncluded(DiagnosticNames::diag_mushyAverageBulkConc))
    {

    Real mushVol;
    Real mushAvBulkC = averageOverMushyRegion(ScalarVars::m_bulkConcentration, mushVol);
//    Real mushAvBulkC = 0.0;
    Real mushAvPorosity = averageOverMushyRegion(ScalarVars::m_porosity, mushVol);

    // Convert mush volume in number of cells to physical volume (in dimensionless units)
    mushVol *= pow(m_dx, SpaceDim);

    m_diagnostics.addDiagnostic(DiagnosticNames::diag_mushyAverageBulkConc, m_time, mushAvBulkC);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_mushyAveragePorosity, m_time, mushAvPorosity);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_mushyVol, m_time, mushVol);


    }

  // Now lets work out some chimney geometry:
  // - how big
  // - what spacing
  bool doChimneyDiagnostics = false; // turn this off for now as it seems to be very slow
  if (m_level == 0 && doChimneyDiagnostics)
  {
    computeChimneyDiagnostics();
  }


}

Real AMRLevelMushyLayer::computeMushDepth(Real a_porosity_criteria)
{
  Vector<Real> averagedPorosity;
  horizontallyAverage(averagedPorosity, *m_scalarNew[ScalarVars::m_porosity]);

  int depth_i = 0;
  Real depth = -1.0;
  for (int i = 0; i < averagedPorosity.size() ; i++)
  {
    depth_i++;
    if (averagedPorosity[i] < a_porosity_criteria)
    {
      depth = (averagedPorosity.size()-depth_i)*m_dx;
      break;
    }
  }

  return depth;
}
