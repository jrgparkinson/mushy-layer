#include "computeSum.H"
#include "AMRLevelMushyLayer.H"
#include "analyticSolns.H"
#include "SetValLevel.H"

/**
 * This source file contains methods for initialising and defining objects
 */

void AMRLevelMushyLayer::setDefaults()
{
  // Need to make sure these are initialised,
  // else  Some conditional statements will depends on uninitialised value(s)
  m_newGrids_different = false;
  m_prev_diag_output = m_time;

}


void AMRLevelMushyLayer::define(AMRLevel* a_coarserLevelPtr,
                                const ProblemDomain& a_problemDomain, int a_level, int a_refRatio)
{

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::define (level " << m_level << ")" << endl;
  }

  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr, a_problemDomain, a_level, a_refRatio);

  if (a_coarserLevelPtr != NULL)
  {
    AMRLevelMushyLayer* amrMLcrse =
        dynamic_cast<AMRLevelMushyLayer*>(a_coarserLevelPtr);

    if (amrMLcrse != NULL)
    {
      define(amrMLcrse->m_opt, amrMLcrse->m_parameters);
    }
    else
    {
      MayDay::Error(
          "AMRLevelMushyLayer::define: a_coarserLevelPtr is not castable to AMRLevelMushyLayer*");
    }
  }

  setDefaults();

  // Why wasn't this line here before?
  m_problem_domain = a_problemDomain;

  m_numCells = m_problem_domain.domainBox().size();

  m_domainHeight = m_opt.domainWidth*m_numCells[SpaceDim-1]/m_numCells[0];

  // Compute the grid spacing
  m_dx = m_opt.domainWidth / m_numCells[0];

  m_numGhost = 1;
  m_numGhostAdvection = 4;

  m_scalarVarNames = Vector<string>(m_numScalarVars, string("scalar"));
  m_scalarVarNames[ScalarVars::m_enthalpy] = string("Enthalpy");
  m_scalarVarNames[ScalarVars::m_bulkConcentration] = string("Bulk concentration");
  m_scalarVarNames[ScalarVars::m_temperature] = string("Temperature");
  m_scalarVarNames[ScalarVars::m_porosity] = string("Porosity");
  m_scalarVarNames[ScalarVars::m_liquidConcentration] = string("Liquid concentration");
  m_scalarVarNames[ScalarVars::m_solidConcentration] = string("Solid concentration");
  m_scalarVarNames[ScalarVars::m_pressure] = string("Pressure");
  m_scalarVarNames[ScalarVars::m_permeability] = string("Permeability");
  m_scalarVarNames[ScalarVars::m_viscosity] = string("Viscosity");
  m_scalarVarNames[ScalarVars::m_lambda] = string("lambda");
  m_scalarVarNames[ScalarVars::m_lambda_porosity] = string("lambda_porosity");
  m_scalarVarNames[ScalarVars::m_enthalpySolidus] = string("Enthalpy solidus");
  m_scalarVarNames[ScalarVars::m_enthalpyLiquidus] = string("Enthalpy liquidus");
  m_scalarVarNames[ScalarVars::m_enthalpyEutectic] = string("Enthalpy eutectic");
  m_scalarVarNames[ScalarVars::m_temperatureAnalytic] = string("T analytic");
  m_scalarVarNames[ScalarVars::m_porosityAnalytic] = string("Porosity analytic");
  m_scalarVarNames[ScalarVars::m_saltEqnSrcGodunov] = string("Salt src term godunov");

  m_scalarVarNames[ScalarVars::m_Terr] = string("T err");
  m_scalarVarNames[ScalarVars::m_enthalpySrc] = string("H src");
  m_scalarVarNames[ScalarVars::m_divUadv] = string("div U face");
  m_scalarVarNames[ScalarVars::m_dHdt] = string("dHdt");
  m_scalarVarNames[ScalarVars::m_dSdt] = string("dSdt");
  m_scalarVarNames[ScalarVars::m_averageVerticalFlux] = string("average vertical solute flux");
  m_scalarVarNames[ScalarVars::m_soluteFluxAnalytic] = string("analytic vertical solute flux");
  m_scalarVarNames[ScalarVars::m_verticalFlux] = string("vertical solute flux");
  m_scalarVarNames[ScalarVars::m_divU] = string("div U cell");
  m_scalarVarNames[ScalarVars::m_averageHeatFlux] = string("average vertical heat flux");
  m_scalarVarNames[ScalarVars::m_streamfunction] = string("streamfunction");
  m_scalarVarNames[ScalarVars::m_vorticity] = string("vorticity");
  m_scalarVarNames[ScalarVars::m_FsVertFluid] = string("FsVertFluid");
  m_scalarVarNames[ScalarVars::m_FsVertFrame] = string("FsVertFrame");
  m_scalarVarNames[ScalarVars::m_FsVertDiffusion] = string("FsVertDiffusion");
  m_scalarVarNames[ScalarVars::m_divUcorr] = string("Div U correction");
  m_scalarVarNames[ScalarVars::m_pi] = string("pi");
  m_scalarVarNames[ScalarVars::m_phi] = string("phi");
  m_scalarVarNames[ScalarVars::m_MACBC] = string("MAC projection BC");
  m_scalarVarNames[ScalarVars::m_MACrhs] = string("MAC projection rhs");
  m_scalarVarNames[ScalarVars::m_CCrhs] = string("CC projection rhs");

  m_vectorVarNames = Vector<string>(m_numVectorVars, string("vector"));
  m_vectorVarNames[VectorVars::m_fluidVel] = string("Darcy velocity");
  m_vectorVarNames[VectorVars::m_U_porosity] = string("U divided by porosity");
  m_vectorVarNames[VectorVars::m_Ustar] = string("U star");
  m_vectorVarNames[VectorVars::m_advUstar] = string("adv U star");
  m_vectorVarNames[VectorVars::m_advectionVel] = string("Advection velocity");
  m_vectorVarNames[VectorVars::m_viscousSolveSrc] = string("Viscous solve src");
  m_vectorVarNames[VectorVars::m_UdelU] = string("u del U");
  m_vectorVarNames[VectorVars::m_advectionSrc] = string("Advection Source");
  m_vectorVarNames[VectorVars::m_fluidVelAnalytic] = string("Darcy Vel Analytic");
  m_vectorVarNames[VectorVars::m_fluidVelErr] = string("Darcy Vel Error");
  //  m_vectorVarNames[m_fluidRefluxCorr] = string("Vel reflux correction");
  m_vectorVarNames[VectorVars::m_dUdt] = string("dUdt");
  m_vectorVarNames[VectorVars::m_FsDiffusion] = string("FsDiffusion");
  m_vectorVarNames[VectorVars::m_FsFluid] = string("FsFluid");
  m_vectorVarNames[VectorVars::m_Fs] = string("Fs");
  m_vectorVarNames[VectorVars::m_freestreamCorrection] = string("Freestream correction");
  m_vectorVarNames[VectorVars::m_advUpreProjection] = string("Unprojected advection vel");
  m_vectorVarNames[VectorVars::m_UpreProjection] = string("Unprojected CC vel");
  m_vectorVarNames[VectorVars::m_advectionImplicitSrc] = string("Implicit Advection Src");
  m_vectorVarNames[VectorVars::m_MACcorrection] = string("MAC correction");
  m_vectorVarNames[VectorVars::CCcorrection] = string("CC correction");
  //  m_vectorVarNames[m_lapUinit] = string("Lap U init");
  //  m_vectorVarNames[m_lapUfinal] = string("Lap U final");
  m_vectorVarNames[VectorVars::m_advSrcLapU] = string("AdvSrc Lap U");

  m_vectorVarNames[VectorVars::m_bodyForce] = string("Body force");
  m_vectorVarNames[VectorVars::m_advVelCorr] = string("Explicit advVel corr");


  if (m_opt.minimalOutput)
  {
    m_outputScalarVars.push_back(ScalarVars::m_enthalpy);
    m_outputScalarVars.push_back(ScalarVars::m_bulkConcentration);
  }
  else
  {
    // These are the default output vars
    m_outputScalarVars.push_back(ScalarVars::m_enthalpy);
    m_outputScalarVars.push_back(ScalarVars::m_bulkConcentration);
    m_outputScalarVars.push_back(ScalarVars::m_temperature);
    m_outputScalarVars.push_back(ScalarVars::m_porosity);
    m_outputScalarVars.push_back(ScalarVars::m_liquidConcentration);

    m_outputScalarVars.push_back(ScalarVars::m_streamfunction);
    m_outputScalarVars.push_back(ScalarVars::m_permeability);
    m_outputScalarVars.push_back(ScalarVars::m_lambda);
    //    m_outputScalarVars.push_back(ScalarVars::m_lambda_porosity);

    m_outputScalarVars.push_back(ScalarVars::m_pressure);

    if (m_parameters.m_viscosityFunction != ViscosityFunction::uniformViscosity)
    {
      m_outputScalarVars.push_back(ScalarVars::m_viscosity);
    }


    if (m_opt.debug)
    {


      m_outputScalarVars.push_back(ScalarVars::m_solidConcentration);

      m_outputScalarVars.push_back(ScalarVars::m_saltEqnSrcGodunov);
      if (m_parameters.physicalProblem == PhysicalProblems::m_solidificationNoFlow)
      {
        m_outputScalarVars.push_back(ScalarVars::m_temperatureAnalytic);
        m_outputScalarVars.push_back(ScalarVars::m_porosityAnalytic);
        m_outputScalarVars.push_back(ScalarVars::m_Terr);
      }

      m_outputScalarVars.push_back(ScalarVars::m_enthalpySrc);
      m_outputScalarVars.push_back(ScalarVars::m_divUadv);
      m_outputScalarVars.push_back(ScalarVars::m_dHdt);
      m_outputScalarVars.push_back(ScalarVars::m_averageVerticalFlux);
      m_outputScalarVars.push_back(ScalarVars::m_dSdt);
      m_outputScalarVars.push_back(ScalarVars::m_soluteFluxAnalytic);
      m_outputScalarVars.push_back(ScalarVars::m_verticalFlux);
      m_outputScalarVars.push_back(ScalarVars::m_divU);
      m_outputScalarVars.push_back(ScalarVars::m_averageHeatFlux);

      m_outputScalarVars.push_back(ScalarVars::m_vorticity);
      m_outputScalarVars.push_back(ScalarVars::m_FsVertFrame);
      m_outputScalarVars.push_back(ScalarVars::m_FsVertFluid);
      m_outputScalarVars.push_back(ScalarVars::m_FsVertDiffusion);
      //  m_outputScalarVars.push_back(ScalarVars::m_divUcorr);
      m_outputScalarVars.push_back(ScalarVars::m_pi);
      m_outputScalarVars.push_back(ScalarVars::m_phi);
      m_outputScalarVars.push_back(ScalarVars::m_MACBC);

      m_outputScalarVars.push_back(ScalarVars::m_MACrhs);
      m_outputScalarVars.push_back(ScalarVars::m_CCrhs);

    }

  }

  if (m_opt.minimalOutput)
  {
    m_outputVectorVars.push_back(VectorVars::m_advectionVel);
  }
  else
  {
    // These are the default output vars

    m_outputVectorVars.push_back(VectorVars::m_fluidVel);
    m_outputVectorVars.push_back(VectorVars::m_advectionVel);

    //    m_outputVectorVars.push_back(VectorVars::m_FsDiffusion);
    //    m_outputVectorVars.push_back(VectorVars::m_FsFluid);
    m_outputVectorVars.push_back(VectorVars::m_Fs);


    if (m_opt.debug)
    {
      m_outputVectorVars.push_back(VectorVars::m_FsDiffusion);
      m_outputVectorVars.push_back(VectorVars::m_FsFluid);


      //  m_outputVectorVars.push_back(VectorVars::m_U_porosity);
      m_outputVectorVars.push_back(VectorVars::m_Ustar);
      m_outputVectorVars.push_back(VectorVars::m_advUstar);

      //  m_outputVectorVars.push_back(VectorVars::m_fluidVelAnalytic);
      //  m_outputVectorVars.push_back(VectorVars::m_fluidVelErr);
      //  m_outputVectorVars.push_back(m_fluidRefluxCorr);
      m_outputVectorVars.push_back(VectorVars::m_advectionSrc);
      m_outputVectorVars.push_back(VectorVars::m_viscousSolveSrc);
      //  m_outputVectorVars.push_back(VectorVars::m_dUdt);
      m_outputVectorVars.push_back(VectorVars::m_freestreamCorrection);
      m_outputVectorVars.push_back(VectorVars::m_advUpreProjection);
      m_outputVectorVars.push_back(VectorVars::m_UpreProjection);
      m_outputVectorVars.push_back(VectorVars::m_advectionImplicitSrc);
      m_outputVectorVars.push_back(VectorVars::CCcorrection);
      m_outputVectorVars.push_back(VectorVars::m_MACcorrection);

      m_outputVectorVars.push_back(VectorVars::m_UdelU);
      m_outputVectorVars.push_back(VectorVars::m_U_porosity);

      m_outputVectorVars.push_back(VectorVars::m_advVelCorr);

      //        m_outputVectorVars.push_back(VectorVars::m_advSrcLapU);
      //  m_outputVectorVars.push_back(m_lapUfinal);


    }

  }

  m_numOutputComps = m_outputScalarVars.size() + m_outputVectorVars.size()*SpaceDim;

  // Now sort out what we need for checkpoint files
  // This should be a fairly small list - trying to make these files as small as possible
  // to save storage space on the disk. Should only output the bare minimum needed to
  // restart simulations
  m_chkVectorVars.push_back(VectorVars::m_fluidVel);

  m_chkScalarVars.push_back(ScalarVars::m_enthalpy);
  m_chkScalarVars.push_back(ScalarVars::m_bulkConcentration);
  m_chkScalarVars.push_back(ScalarVars::m_pressure);
  m_chkScalarVars.push_back(ScalarVars::m_lambda);

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::define - finished (level " << m_level << ")" << endl;
  }
}


void AMRLevelMushyLayer::defineCFInterp()
{
  AMRLevelMushyLayer* mlCrse = getCoarserLevel();

  if (mlCrse)
  {
    DisjointBoxLayout* crseGridsPtr = &(mlCrse->m_grids);
    int nRefCrse = mlCrse->m_ref_ratio;

    m_quadCFInterpScalar.define(m_grids, crseGridsPtr, m_dx, nRefCrse, 1,
                                m_problem_domain);
    m_quadCFInterpVector.define(m_grids, crseGridsPtr, m_dx, nRefCrse, SpaceDim,
                                m_problem_domain);


    //    AMRLevelMushyLayer* crseML = getCoarserLevel();
    const ProblemDomain& crseDomain = mlCrse->problemDomain();

    m_piecewiseLinearFillPatchScalarOne.define(m_grids, *crseGridsPtr,  1, crseDomain, nRefCrse,
                                               1, false);
    m_piecewiseLinearFillPatchScalarTwo.define(m_grids, *crseGridsPtr,  1, crseDomain, nRefCrse,
                                               2, false);
    m_piecewiseLinearFillPatchScalarThree.define(m_grids, *crseGridsPtr,  1, crseDomain, nRefCrse,
                                                 3, false);
    m_piecewiseLinearFillPatchScalarFour.define(m_grids, *crseGridsPtr,  1, crseDomain, nRefCrse,
                                                4, false);

  }
}

void AMRLevelMushyLayer::levelSetup()
{
  CH_TIME("AMRLevelMushyLayer::levelSetup");

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::levelSetup (level " << m_level << ")" << endl;
  }

  // In case these have changed
  defineIBCs();

  m_physBCPtr->Dx(m_dx);

  AMRLevelMushyLayer* amrMLCoarserPtr = getCoarserLevel();
  AMRLevelMushyLayer* amrMLFinerPtr = getFinerLevel();

  m_hasCoarser = (amrMLCoarserPtr != NULL);
  m_hasFiner = (amrMLFinerPtr != NULL);

  int nRefCrse = -1;

  DisjointBoxLayout* crseGridsPtr = NULL;

  Projector *crsProj = NULL;
  Projector *fineProj = NULL;
  Projector *crsProjBackup = NULL;
  Projector *fineProjBackup = NULL;
  //    DisjointBoxLayout *crseGrids = NULL;

  bool scaleFineFluxes = true;

  if (m_hasFiner)
  {
    fineProj = &(amrMLFinerPtr->m_projection);
    fineProjBackup = &(amrMLFinerPtr->m_projectionBackup);
  }
  if (m_hasCoarser)
  {
    crsProj = &(amrMLCoarserPtr->m_projection);
    crsProjBackup = &(amrMLCoarserPtr->m_projectionBackup);
    crseGridsPtr = &(amrMLCoarserPtr->m_grids);
    nRefCrse = amrMLCoarserPtr->m_ref_ratio;
  }

  LevelDomainFluxRegister* fineDomainFRheat = NULL;
  LevelDomainFluxRegister* coarseDomainFRheat = NULL;
  LevelDomainFluxRegister* fineDomainFRsalt = NULL;
  LevelDomainFluxRegister* coarseDomainFRsalt = NULL;

  if (m_hasCoarser)
  {
    coarseDomainFRheat = &(getCoarserLevel()->m_heatDomainFluxRegister);
    coarseDomainFRsalt = &(getCoarserLevel()->m_saltDomainFluxRegister);
  }

  if (m_hasFiner)
  {
    fineDomainFRheat = &(getFinerLevel()->m_heatDomainFluxRegister);
    fineDomainFRsalt = &(getFinerLevel()->m_saltDomainFluxRegister);
  }

  m_saltDomainFluxRegister.define(m_problem_domain, m_grids, m_ref_ratio, m_dx,
                                  fineDomainFRsalt, coarseDomainFRsalt);

  m_heatDomainFluxRegister.define(m_problem_domain, m_grids, m_ref_ratio, m_dx,
                                  fineDomainFRheat, coarseDomainFRheat);

  if (m_hasCoarser)
  {
    nRefCrse = m_coarser_level_ptr->refRatio();

    //    const DisjointBoxLayout& coarserLevelDomain = amrMLCoarserPtr->m_grids;

    m_coarseAverageScalar.define(m_grids, 1, nRefCrse);
    m_coarseAverageHC.define(m_grids, 2, nRefCrse);

    m_coarseAverageVector.define(m_grids, SpaceDim, nRefCrse);

    m_fineInterpScalar.define(m_grids, 1, nRefCrse, m_problem_domain);

    m_fineInterpVector.define(m_grids, SpaceDim, nRefCrse,
                              m_problem_domain);



    // This may look twisted but you have to do this this way because the
    // coarser levels get setup before the finer levels so, since a flux
    // register lives between this level and the next FINER level, the finer
    // level has to do the setup because it is the only one with the
    // information at the time of construction.

    // Maintain flux registers
    for (int a_scalarVar = 0; a_scalarVar < m_numScalarVars;
        a_scalarVar++)
    {
      if (m_makeFluxRegForScalarVar[a_scalarVar])
      {
        amrMLCoarserPtr->m_fluxRegisters[a_scalarVar] = RefCountedPtr<
            LevelFluxRegister>(
                new LevelFluxRegister(m_grids, amrMLCoarserPtr->m_grids,
                                      m_problem_domain, amrMLCoarserPtr->m_ref_ratio,
                                      1, scaleFineFluxes));
        amrMLCoarserPtr->m_fluxRegisters[a_scalarVar]->setToZero();
      }
    }

    amrMLCoarserPtr->m_fluxRegHC = RefCountedPtr<
        LevelFluxRegister>(
            new LevelFluxRegister(m_grids, amrMLCoarserPtr->m_grids,
                                  m_problem_domain, amrMLCoarserPtr->m_ref_ratio,
                                  2, scaleFineFluxes));
    amrMLCoarserPtr->m_fluxRegHC->setToZero();

    for (int a_vectorVar = 0; a_vectorVar < m_numVectorVars;
        a_vectorVar++)
    {
      if (m_makeFluxRegForVectorVar[a_vectorVar])
      {
        amrMLCoarserPtr->m_vectorFluxRegisters[a_vectorVar] =
            RefCountedPtr<LevelFluxRegister>(
                new LevelFluxRegister(m_grids,
                                      amrMLCoarserPtr->m_grids,
                                      m_problem_domain,
                                      amrMLCoarserPtr->m_ref_ratio,
                                      SpaceDim));
        amrMLCoarserPtr->m_vectorFluxRegisters[a_vectorVar]->setToZero();
      }
    }



    // Some of our level structures can't be setup until we've created grids on this level
    Vector<Box> b = m_grids.boxArray();
    if (m_grids.boxArray().size() > 0)
    {
      m_projection.verbosity(m_opt.projection_verbosity);
      m_projection.define(m_grids, crseGridsPtr, m_problem_domain, m_dx,
                          fineProj, crsProj, nRefCrse, m_level, *m_physBCPtr, m_opt.usePiAdvectionBCs);

      m_projectionBackup.verbosity(m_opt.projection_verbosity);
      m_projectionBackup.define(m_grids, crseGridsPtr, m_problem_domain, m_dx,
                                fineProjBackup, crsProjBackup, nRefCrse, m_level, *m_physBCPtr, m_opt.usePiAdvectionBCs);

    }
    else
    {
      // Still need to set refratio on finer projection
      //                        m_projection.m_nRefCrse = nRefCrse;
    }


  } // end if has coarser
  else
  {
    // on the coarsest level here

    m_projection.define(m_grids, crseGridsPtr, m_problem_domain, m_dx,
                        fineProj, crsProj, m_ref_ratio, m_level, *m_physBCPtr, m_opt.usePiAdvectionBCs);
    m_projection.verbosity(m_opt.projection_verbosity);

    m_projectionBackup.define(m_grids, crseGridsPtr, m_problem_domain, m_dx,
                              fineProj, crsProj, m_ref_ratio, m_level, *m_physBCPtr, m_opt.usePiAdvectionBCs);
    m_projectionBackup.verbosity(m_opt.projection_verbosity);

  }

  if (!hasFinerLevel())
  {
    m_projection.isFinestLevel(true);
    m_projectionBackup.isFinestLevel(true);


  }

  defineCFInterp();

  //Some things don't care about coarser/finer levels, define them here:


  m_advectionPhysicsVelocity.define(m_problem_domain, m_dx);
  m_advectionPhysicsVelocity.setNComp(SpaceDim);

  m_advPhysHC.define(m_problem_domain, m_dx);
  m_advPhysHC.setNComp(2);


  m_advPhysTSl.define(m_problem_domain, m_dx);
  m_advPhysTSl.setNComp(2);

  m_advPhysH.define(m_problem_domain, m_dx);
  m_advPhysH.setNComp(1);
  m_advPhysC.define(m_problem_domain, m_dx);
  m_advPhysC.setNComp(1);

  m_advPhysT.define(m_problem_domain, m_dx);
  m_advPhysT.setNComp(1);
  m_advPhysSl.define(m_problem_domain, m_dx);
  m_advPhysSl.setNComp(1);

  // Collect all the BC stuff here
  setAdvectionBCs();

  // Don't do slope limiting. Not possible to do more than one of these at the same time

  // Primitive limiting stops 2nd order convergence?
  // However definitely need it to do advection properly. Without we get porosities > 1 !
  bool usePrimLimiting = true && m_opt.useLimiting;
  bool useCharLimiting = false && m_opt.useLimiting;
  bool useFlattening = false; // Can't do this

  m_patchGodVelocity.define(m_problem_domain, m_dx, &m_advectionPhysicsVelocity,
                            m_opt.velAdvNormalPredOrder, m_opt.velAdvUseFourthOrderSlopes, usePrimLimiting,
                            useCharLimiting, useFlattening, m_opt.velAdvUseArtVisc, m_opt.velAdvArtVisc);

  // Trying to fix issue with slopes at internal box boundaries
  m_patchGodVelocity.highOrderLimiter(m_opt.velAdvHigherOrderLimiter);

  // For scalars
  bool usePrimLimitingHC = usePrimLimiting;
  bool useCharLimitingHC = useCharLimiting;
  bool useFlatteningHC = useFlattening;

  // this is crucial for getting convergence.
  usePrimLimitingHC = m_opt.useLimiting;

  m_patchGodHC.define(m_problem_domain, m_dx, &m_advPhysHC,
                      m_opt.HCNormalPredOrder, m_opt.HCUseFourthOrderSlopes, usePrimLimitingHC,
                      useCharLimitingHC, useFlatteningHC, m_opt.HCUseArtVisc, m_opt.HCArtVisc);
  m_patchGodHC.highOrderLimiter(m_opt.HCHigherOrderLimiter);

  m_patchGodTSl.define(m_problem_domain, m_dx, &m_advPhysTSl,
                       m_opt.HCNormalPredOrder, m_opt.HCUseFourthOrderSlopes, usePrimLimitingHC,
                       useCharLimitingHC, useFlatteningHC, m_opt.HCUseArtVisc,  m_opt.HCArtVisc);
  m_patchGodTSl.highOrderLimiter(m_opt.HCHigherOrderLimiter);

  // Single component solve
  m_patchGodH.define(m_problem_domain, m_dx, &m_advPhysH,
                     m_opt.HCNormalPredOrder, m_opt.HCUseFourthOrderSlopes, usePrimLimitingHC,
                     useCharLimitingHC, useFlatteningHC, m_opt.HCUseArtVisc,  m_opt.HCArtVisc);
  m_patchGodC.define(m_problem_domain, m_dx, &m_advPhysC,
                     m_opt.HCNormalPredOrder, m_opt.HCUseFourthOrderSlopes, usePrimLimitingHC,
                     useCharLimitingHC, useFlatteningHC, m_opt.HCUseArtVisc,  m_opt.HCArtVisc);

  m_patchGodT.define(m_problem_domain, m_dx, &m_advPhysT,
                     m_opt.HCNormalPredOrder, m_opt.HCUseFourthOrderSlopes, usePrimLimitingHC,
                     useCharLimitingHC, useFlatteningHC, m_opt.HCUseArtVisc,  m_opt.HCArtVisc);
  m_patchGodSl.define(m_problem_domain, m_dx, &m_advPhysSl,
                      m_opt.HCNormalPredOrder, m_opt.HCUseFourthOrderSlopes, usePrimLimitingHC,
                      useCharLimitingHC, useFlatteningHC, m_opt.HCUseArtVisc,  m_opt.HCArtVisc);


  for (int var = 0; var < m_numScalarVars; var++)
  {
    if (m_makeFluxRegForScalarVar[var])
    {

      AdvectionPhysics advectionPhysicsScalar;

      advectionPhysicsScalar.define(m_problem_domain,m_dx);
      advectionPhysicsScalar.setNComp(1);

      advectionPhysicsScalar.setPhysIBC(m_scalarIBC[var]);

      m_patchGodScalars[var] = RefCountedPtr<PatchGodunov>(new PatchGodunov());
      m_patchGodScalars[var]->define(m_problem_domain,
                                     m_dx,
                                     &advectionPhysicsScalar,
                                     m_opt.velAdvNormalPredOrder,
                                     m_opt.velAdvUseFourthOrderSlopes,
                                     usePrimLimiting,
                                     useCharLimiting,
                                     useFlattening,
                                     m_opt.velAdvUseArtVisc,
                                     m_opt.velAdvArtVisc);
    }
  }


}


void AMRLevelMushyLayer::defineUstarMultigrid()
{
  CH_TIME("AMRLevelMushyLayer::defineUstarMultigrid");

  // Define multigrid solver for this level and coarser level if one exists

  Vector<AMRLevelMushyLayer*> hierarchy;
  Vector<DisjointBoxLayout> allGrids, solverGrids;
  Vector<int> refRat;
  ProblemDomain lev0Dom;
  Real lev0Dx;
  IntVect ivGhost = IntVect::Unit;
  getHierarchyAndGrids(hierarchy, allGrids, refRat, lev0Dom, lev0Dx);

  Real old_time = m_time-m_dt;

  int nlevels = allGrids.size();
  int coarsestLevel  = 0;

  Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef;
  Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > cCoef;

  aCoef.resize(nlevels);
  bCoef.resize(nlevels);
  cCoef.resize(nlevels);

  // I'm not sure we actually need this
  solverGrids.resize(nlevels);

//  bool is_time_dependent = true;
//  if (!m_opt.doEulerPart)
//  {
//    is_time_dependent = false;
//    is_time_dependent = true;
//  }

  // If we're doing a multi-component solve, need SpaceDim components, otherwise just one component.
  int num_comp = m_opt.multiCompUStarSolve ? SpaceDim : 1;

  for (int lev = coarsestLevel; lev < nlevels; lev++)
  {

    int relativeLev = lev - coarsestLevel;

    AMRLevelMushyLayer* thisLevel =
        (AMRLevelMushyLayer*) (hierarchy[lev]);
    DisjointBoxLayout& levelGrids = thisLevel->m_grids;

    solverGrids[relativeLev] = levelGrids;

    bCoef[relativeLev] = RefCountedPtr<LevelData<FluxBox> >(
        new LevelData<FluxBox>(levelGrids, num_comp, ivGhost)); // = 1
    aCoef[relativeLev] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(levelGrids, num_comp, ivGhost)); // = 1
    cCoef[relativeLev] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(levelGrids, num_comp, ivGhost)); // = porosity.pi_0/(pi)


    // Only actually fill levels if they're at m_level or coarser
    if (lev <= m_level)
    {
      LevelData<FArrayBox> porosity(levelGrids, 1, ivGhost);
      LevelData<FArrayBox> permeability(levelGrids, 1, ivGhost);

      // Try and lag this for stability (was previously new_time)
      Real coeff_time = old_time;
      thisLevel->fillScalars(porosity, coeff_time, m_porosity, true);
      thisLevel->fillScalars(permeability, coeff_time, m_permeability, true);

      //Make aCoef
      DataIterator dit = aCoef[relativeLev]->dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
      {
        // do this or darcy as a source term

        // this is what multiplies the U term (darcy's law)
//        if (m_opt.explicitDarcyTerm)
//        {
//          (*cCoef[relativeLev])[dit].setVal(0.0);
//        }
//        else
//        {
          for (int comp_i = 0; comp_i < num_comp; comp_i++)
          {
            (*cCoef[relativeLev])[dit].copy(porosity[dit], 0, comp_i);
            (*cCoef[relativeLev])[dit].divide(permeability[dit], 0, comp_i);
          }

          (*cCoef[relativeLev])[dit].mult(m_parameters.m_darcyCoeff);
//        }

        //                                                      (*cCoef[relativeLev])[dit].setVal(0.0); // testing

        // This is what multiplies du/dt
          if (isVelocityTimeDependent())
          {
            (*aCoef[relativeLev])[dit].setVal(1.0);
          }
          else
          {
            // no time dependence
            (*aCoef[relativeLev])[dit].setVal(0.0);
          }

        // this is what multiplies laplacian(u). Note the minus sign.
        (*bCoef[relativeLev])[dit].setVal(- m_parameters.m_viscosityCoeff);
      } // end loop over boxes

    } // end if level <= m_level
  } // end loop over levels

  if (m_opt.multiCompUStarSolve)
  {
    MayDay::Error("Multcomponent U start solve not implemented yet.");
//    BCHolder viscousBC = m_physBCPtr->enthalpySalinityBC()  uStarFuncBC(m_viscousBCs); //  ->velFuncBC(idir, m_viscousBCs);
//
//    RefCountedPtr<DarcyBrinkmanOpFactory> vcamrpop = RefCountedPtr<DarcyBrinkmanOpFactory>(new DarcyBrinkmanOpFactory());
//    vcamrpop->define(lev0Dom, allGrids, refRat, lev0Dx, viscousBC,
//                     0.0, aCoef, -1.0, bCoef, cCoef); // Note that we should set m_dt*etc in bCoef, not beta!
//
//    s_uStarOpFactMultiComp = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(vcamrpop); // m_UstarVCAMRPOp[idir]);
//
//    s_uStarAMRMGMultiComp->define(lev0Dom, *s_uStarOpFactMultiComp,
//                                  &s_botSolverUStar, nlevels);
//
//    s_uStarAMRMGMultiComp->setSolverParameters(numSmooth, numSmooth, numSmooth,
//                                               numMG, maxIter, tolerance, hang, normThresh);
  }
  else
  {

    for (int idir = 0; idir < SpaceDim; idir++)
    {

      BCHolder viscousBC = m_physBCPtr->velFuncBC(idir, m_opt.viscousBCs);

      RefCountedPtr<DarcyBrinkmanOpFactory> vcamrpop = RefCountedPtr<DarcyBrinkmanOpFactory>(new DarcyBrinkmanOpFactory());
      vcamrpop->define(lev0Dom, allGrids, refRat, lev0Dx, viscousBC,
                       0.0, aCoef, -1.0, bCoef, cCoef); // Note that we should set m_dt*etc in bCoef, not beta!

      m_uStarOpFact[idir] = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(vcamrpop); // m_UstarVCAMRPOp[idir]);

      m_uStarAMRMG[idir]->define(lev0Dom, *m_uStarOpFact[idir],
                                 &s_botSolverUStar, nlevels);

      m_uStarAMRMG[idir]->setSolverParameters(m_opt.velMGNumSmooth, m_opt.velMGNumSmooth, m_opt.velMGNumSmooth,
                                              m_opt.velMGNumMG, m_opt.VelMGMaxIter, m_opt.velMGTolerance, m_opt.velMGHang, m_opt.velMGNormThresh);
    }
  }
}

void AMRLevelMushyLayer::defineUstarSolver(     Vector<RefCountedPtr<LevelBackwardEuler> >& UstarBE,
                                                Vector<RefCountedPtr<LevelTGA> >& UstarTGA)
{
  CH_TIME("AMRLevelMushyLayer::defineUstarSolver");

  Vector<AMRLevelMushyLayer*> hierarchy;
  Vector<DisjointBoxLayout> allGrids;
  Vector<int> refRat;
  ProblemDomain lev0Dom;
  Real lev0Dx;
  getHierarchyAndGrids(hierarchy, allGrids, refRat, lev0Dom, lev0Dx);

  defineUstarMultigrid();

  UstarBE.resize(SpaceDim);
  UstarTGA.resize(SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    // Now define Backward Euler for timestepping
    UstarBE[idir] = RefCountedPtr<LevelBackwardEuler> (new LevelBackwardEuler(allGrids, refRat, lev0Dom, m_uStarOpFact[idir], m_uStarAMRMG[idir]));
    UstarTGA[idir] = RefCountedPtr<LevelTGA> (new LevelTGA(allGrids, refRat, lev0Dom, m_uStarOpFact[idir], m_uStarAMRMG[idir]));
  }

}


void AMRLevelMushyLayer::defineSolvers(Real a_time)
{

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::defineSolvers" << endl;
  }
  CH_TIME("AMRLevelMushyLayer::defineSolvers");

  //  bool a_homogeneous = false;
  Real alpha = 1.0;
  Real beta = 1.0;

  IntVect ivGhost = m_numGhost * IntVect::Unit;

  s_botSolverHC.m_imax = m_opt.HCMultigridBottomSolveIterations;

  Vector<AMRLevelMushyLayer*> hierarchy;
  Vector<DisjointBoxLayout> grids;
  Vector<int> refRat;
  ProblemDomain lev0Dom;
  Real lev0Dx;
  getHierarchyAndGrids(hierarchy, grids, refRat, lev0Dom, lev0Dx);

  int numLevels = grids.size();

  const DisjointBoxLayout* crseGridsPtr = NULL;
  int nRefCrse = -1;

  if (m_level > 0)
  {
    AMRLevelMushyLayer* coarserAMRLevel = getCoarserLevel();
    crseGridsPtr = &(coarserAMRLevel->m_grids);
    nRefCrse = coarserAMRLevel->m_ref_ratio;
  }

  s_botSolverUStar.m_verbosity = max(m_opt.HCMultigridVerbosity - 2, 0);
  s_botSolverHC.m_verbosity = max(m_opt.HCMultigridVerbosity - 2, 0);

  {

    CH_TIME("AMRLevelMushyLayer::defineSolvers::defineViscousOp");
    for (int dir = 0; dir < SpaceDim; dir++)
    {
      BCHolder viscousBCcomp = m_physBCPtr->velFuncBC(dir, m_opt.viscousBCs);

      m_viscousOp[dir] = RefCountedPtr<AMRPoissonOp>(new AMRPoissonOp());
      m_viscousOp[dir]->define(m_grids, crseGridsPtr, m_dx, nRefCrse,
                               m_problem_domain, viscousBCcomp);

    }
  }


  {
    CH_TIME("AMRLevelMushyLayer::defineSolvers::defineUStarAMRMG");
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      // I think we only have to define this once
      if (m_uStarAMRMG[idir] == NULL)
      {
        m_uStarAMRMG[idir] =
            RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >(
                new AMRMultiGrid<LevelData<FArrayBox> >());
        m_uStarAMRMG[idir]->setSolverParameters(m_opt.HCMultigridNumSmoothUp, m_opt.HCMultigridNumSmoothUp, m_opt.HCMultigridNumSmoothUp,
                                                m_opt.HCMultigridNumMG, m_opt.HCMultigridMaxIter, m_opt.HCMultigridTolerance, m_opt.HCMultigridHang, m_opt.HCMultigridNormThresh);
        m_uStarAMRMG[idir]->m_verbosity = m_opt.HCMultigridVerbosity;
      }
    }
  }

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::defineSolvers - finished Ustar" << endl;
  }

  Vector<RefCountedPtr<LevelData<FluxBox> > > porosityFace(numLevels);

  Vector<RefCountedPtr<LevelData<FArrayBox> > > porosity(numLevels);

  Vector<RefCountedPtr<LevelData<FArrayBox> > > enthalpy(numLevels);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > bulkConcentration(numLevels);

  Vector<RefCountedPtr<LevelData<FArrayBox> > > enthalpySolidus(numLevels);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > enthalpyLiquidus(numLevels);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > enthalpyEutectic(numLevels);

  //Get coarsest level
  AMRLevelMushyLayer* amrML = this;
  while(amrML->getCoarserLevel())
  {
    amrML = amrML->getCoarserLevel();
  }

  // Fill all levels
  int lev = 0;

  while(amrML != NULL && lev < numLevels)
  {

    porosityFace[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(grids[lev], 1, IntVect::Zero));
    //    porosity[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], 1, ivGhost));

//    enthalpy[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], 1, ivGhost));
//    bulkConcentration[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], 1, ivGhost));

    enthalpySolidus[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], 1, ivGhost));
    enthalpyLiquidus[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], 1, ivGhost));
    enthalpyEutectic[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], 1, ivGhost));

    // Can only fill these at the current level or coarser
    bool secondOrder = true; // Should probably use quadratic CF interp here

    // Only fill at levels which have reached m_time
    // in general, than means this level and coarser
    // however in postTimeStep(), all finer levels will have reached m_time too
    //    if (lev <= m_level)
    if (amrML->m_time <= (m_time + TIME_EPS))
    {
      // Should be new time, as we evaluate operators at new time
      // Actually, make this an argument as it should be different in different cases
      //      Real a_time = m_time;
      // NO! setup solvers at the start of a timestep so have to use fields at a_time-dt
      // else we can't set up the

      amrML->fillScalarFace(*porosityFace[lev], a_time, ScalarVars::m_porosity, true, secondOrder);

//      amrML->fillScalars(*enthalpy[lev], a_time, ScalarVars::m_enthalpy, true , secondOrder);
//      amrML->fillScalars(*bulkConcentration[lev], a_time, ScalarVars::m_bulkConcentration, true, secondOrder);

      amrML->fillScalars(*enthalpySolidus[lev],a_time, ScalarVars::m_enthalpySolidus, true, secondOrder);
      amrML->fillScalars(*enthalpyLiquidus[lev], a_time, ScalarVars::m_enthalpyLiquidus, true, secondOrder);
      amrML->fillScalars(*enthalpyEutectic[lev],a_time, ScalarVars::m_enthalpyEutectic, true ,secondOrder);

    }

    amrML = amrML->getFinerLevel();
    lev = lev + 1;
  }

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::defineSolvers - finished filling scalars" << endl;
  }


  EdgeVelBCHolder porosityEdgeBC(m_physBCPtr->porosityFaceBC());

  // two components: enthalpy and salinity
  int numComps = 2;
  int Hcomp = 0;
  int Ccomp = 1;

  Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(numLevels);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(numLevels);

  for (int lev=0; lev<numLevels; lev++)
  {

    bCoef[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(grids[lev], numComps, ivGhost));
    aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], numComps, ivGhost));

    for (DataIterator dit = bCoef[lev]->dataIterator(); dit.ok(); ++dit)
    {
      (*aCoef[lev])[dit].setVal(1.0);
      (*bCoef[lev])[dit].setVal(1.0);

      // bCoef for salt solve
      (*bCoef[lev])[dit].mult((*porosityFace[lev])[dit], (*porosityFace[lev])[dit].box(), 0, Ccomp);

      // for heat solve
      (*bCoef[lev])[dit].minus((*porosityFace[lev])[dit], (*porosityFace[lev])[dit].box(), 0, Hcomp);
      for (int dir=0; dir<SpaceDim; dir++)
      {
        (*bCoef[lev])[dit][dir].mult(m_parameters.heatConductivityRatio, Hcomp);
      }
      (*bCoef[lev])[dit].plus((*porosityFace[lev])[dit], (*porosityFace[lev])[dit].box(), 0, Hcomp);


      for (int dir=0; dir<SpaceDim; dir++)
      {
        (*bCoef[lev])[dit][dir].mult(-m_scalarDiffusionCoeffs[ScalarVars::m_enthalpy], Hcomp);
        (*bCoef[lev])[dit][dir].mult(-m_scalarDiffusionCoeffs[ScalarVars::m_bulkConcentration], Ccomp);
      }
    }

  } // end loop over levels

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::defineSolvers - finished setting coefficients" << endl;
  }

  MushyLayerParams* mlParamsPtr = &m_parameters;

  //        EnthalpyVariable calcTemperature = computeTemperatureFunc;
  BCHolder temperature_Sl_BC; // = m_physBCPtr->BasicthetaFuncBC();
  temperature_Sl_BC = m_physBCPtr->temperatureLiquidSalinityBC();
  BCHolder HC_BC  = m_physBCPtr->enthalpySalinityBC();



  AMRNonLinearMultiCompOpFactory* HCop = new AMRNonLinearMultiCompOpFactory();
  HCop->define(lev0Dom, grids, refRat, lev0Dx, HC_BC,
               alpha, aCoef, beta, bCoef,
               enthalpySolidus, enthalpyLiquidus, enthalpyEutectic,
               mlParamsPtr, temperature_Sl_BC,
               m_opt.HCMultigridRelaxMode, porosityEdgeBC, m_opt.apply_diagnostic_bcs);

  if (m_opt.nonlinearHCOpSuperOptimised)
  {
    HCop->setSuperOptimised(true);
  }

  m_HCOpFact = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(HCop);

  int maxAMRlevels = hierarchy.size();

  if (m_opt.MGtype == MGmethod::MGTypeStandard)
  {
    MayDay::Error("Standard multigrid no longer supported");
  }
  else if (m_opt.MGtype == MGmethod::MGTypeFAS)
  {
    //FAS multigrid

    m_multiCompFASMG = RefCountedPtr<AMRFASMultiGrid<LevelData<FArrayBox> > >(
        new AMRFASMultiGrid<LevelData<FArrayBox> >());

    if (m_opt.HCMultigridUseRelaxBottomSolverForHC)
    {
      m_multiCompFASMG->define(lev0Dom, *m_HCOpFact, &s_botSolverHC,
                               maxAMRlevels);
    }
    else
    {
      // Borrow the BiCGStab bottom solver from U star
      m_multiCompFASMG->define(lev0Dom, *m_HCOpFact, &s_botSolverUStar,
                               maxAMRlevels);
    }

    m_multiCompFASMG->setSolverParameters(m_opt.HCMultigridNumSmoothDown, m_opt.HCMultigridNumSmoothUp, m_opt.HCMultigridNumSmoothUp, m_opt.HCMultigridNumMG,
                                          m_opt.HCMultigridMaxIter, m_opt.HCMultigridTolerance, m_opt.HCMultigridHang, m_opt.HCMultigridNormThresh);
    m_multiCompFASMG->m_verbosity = m_opt.HCMultigridVerbosity;



    // Not doing TGA yet
    m_enthalpySalinityTGA = RefCountedPtr<LevelTGA>(
        new LevelTGA(grids, refRat, lev0Dom, m_HCOpFact,
                     m_multiCompFASMG));


    m_enthalpySalinityBE = RefCountedPtr<LevelBackwardEuler>(
        new LevelBackwardEuler(grids, refRat, lev0Dom, m_HCOpFact,
                               m_multiCompFASMG));

  }
  else
  {
    MayDay::Error("Unknown multigrid type specified");
  }

}

// Refactored this so we can call it in each timestep if necessary
void AMRLevelMushyLayer::setAdvectionBCs()
{

  CH_TIME("AMRLevelMushyLayer::setAdvectionBCs");
  PhysIBC* velIBC = m_physBCPtr->advectionVelIBC();
  PhysIBC* hcIBC = m_physBCPtr->scalarTraceHC_IBC();
  PhysIBC* TSlIBC = m_physBCPtr->scalarTraceTSl_IBC();
  PhysIBC* H_IBC = m_physBCPtr->scalarTraceH_IBC();
  PhysIBC* C_IBC = m_physBCPtr->scalarTraceC_IBC();
  PhysIBC* T_IBC = m_physBCPtr->scalarTrace_IBC(m_parameters.bcTypeTemperatureHi, m_parameters.bcTypeTemperatureLo,
                                                m_parameters.bcValTemperatureHi, m_parameters.bcValTemperatureLo,
                                                m_parameters.thetaPlumeInflow);
  PhysIBC* Sl_IBC = m_physBCPtr->scalarTrace_IBC(m_parameters.bcTypeLiquidConcentrationHi, m_parameters.bcTypeLiquidConcentrationLo,
                                                 m_parameters.bcValLiquidConcentrationHi, m_parameters.bcValLiquidConcentrationLo,
                                                 m_parameters.ThetaLPlumeInflow);

  m_advectionPhysicsVelocity.setPhysIBC(velIBC);

  m_advPhysHC.setPhysIBC(hcIBC);
  m_advPhysTSl.setPhysIBC(TSlIBC);

  m_advPhysH.setPhysIBC(H_IBC);
  m_advPhysC.setPhysIBC(C_IBC);

  m_advPhysT.setPhysIBC(T_IBC);
  m_advPhysSl.setPhysIBC(Sl_IBC);

  // setPhysIBC takes the PhysIBC pointers and uses them to create new IBCS,
  // so we can now safely delete pointers to prevent a memory leak
  delete velIBC;
  delete hcIBC;
  delete TSlIBC;
  delete H_IBC;
  delete C_IBC;
  delete T_IBC;
  delete Sl_IBC;
}

void
AMRLevelMushyLayer::defineIBCs ()
{
  for (int var = 0; var < m_numScalarVars; var++) {
    if (m_makeFluxRegForScalarVar[var]) {
      m_scalarIBC[var] = getScalarIBCs (var);
    }
  }
}

void AMRLevelMushyLayer::define(MushyLayerOptions a_opt, MushyLayerParams a_params)
{
  m_isDefined = true;

  m_opt = a_opt;

  // AMRLevel owns these
  s_verbosity = a_opt.verbosity;
  m_initial_dt_multiplier = a_opt.initial_dt_multiplier;

  // We should own everything else
  // All the following member variables are things which can change during simulations
  // anything that's static should be wrapped into either
  //    m_opt
  // or m_parameters if it's to do with the model physics
  m_timestepReduced = false;
  m_timestepFailed = false;

  m_adv_vel_centering = 0.5;
  m_dtReduction = -1;

//  m_parameters.getParameters();
  m_parameters = a_params;

  m_physBCPtr = new PhysBCUtil(m_parameters, m_dx);
  m_physBCPtr->setAdvVel(&m_advVel);



  Real diag_timescale = 0;
  if (m_parameters.nonDimVel != 0)
  {
    // Set timescale as mushy layer height (0.1) / V
    diag_timescale = 0.1/abs(m_parameters.nonDimVel);
  }
  else if (m_parameters.rayleighTemp > 0)
  {
    diag_timescale = 0.1/sqrt(m_parameters.rayleighTemp);
  }
  else
  {
    diag_timescale = 1.0;
  }

  m_diagnostics.define(diag_timescale, s_verbosity, m_opt.steadyStateCondition/10);

  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::define - made diagnostics object" << endl;
  }

  /// Porosity for Darcy-Brinkman, permeability for Darcy
  m_pressureScaleVar = solvingFullDarcyBrinkman() ? m_porosity : m_permeability ;

  // This shouldn't be an option.
//  ppMain.query("pressureScaleVar", m_pressureScaleVar);

  s_implicit_reflux = (m_parameters.prandtl > 0);

  // Use inviscid BCs if darcy term likely to dominate
  // This may need more care, particularly re:specific boundaries
//  m_parameters.isViscous() = (m_parameters.m_viscosityCoeff > 0);

  // For mushy layer calculations, we want to try and kick off the instability if we've converged
//  m_doAutomaticRestart = m_opt.initiallyDoAutomaticRestart;


  // Regridding stuff
  m_regrid_smoothing_done = false; // standard initial value


  // There isn't a nondimensional constant we can change to stop doing the U del (U/chi) bit, so have a switch here instead
  //  m_doEulerPart = false;
  //  m_doProjection = true;
  //  m_opt.doSyncOperations = true;

  // This can be changed during initialisation
  m_usePrevPressureForUStar = m_opt.usePrevPressureForUStar;
//  ppMain.query("addSubtractGradP", m_addSubtractGradP);

  //I don't think it makes sense to do projection *and* calculate grad(P)
  m_enforceGradP = (!m_opt.doProjection);


  // Define which scalar fields we want flux registers for
  for (int var = 0; var < m_numScalarVars; var++)
  {
    m_makeFluxRegForScalarVar[var] = false;
    m_scalarDiffusionCoeffs[var] = 0.0;
  }

  // Anything I want to be able to advect
  m_makeFluxRegForScalarVar[ScalarVars::m_lambda] = true;
  m_makeFluxRegForScalarVar[ScalarVars::m_lambda_porosity] = true;

  m_scalarDiffusionCoeffs[ScalarVars::m_temperature] = 0.0;
  m_scalarDiffusionCoeffs[ScalarVars::m_liquidConcentration] = 0.0;
  m_scalarDiffusionCoeffs[ScalarVars::m_bulkConcentration] = 0.0;

  if (!(m_parameters.physicalProblem == PhysicalProblems::m_poiseuilleFlow))
  {
    m_scalarDiffusionCoeffs[ScalarVars::m_enthalpy] = m_parameters.m_heatDiffusionCoeff;
    m_scalarDiffusionCoeffs[ScalarVars::m_bulkConcentration] = m_parameters.m_saltDiffusionCoeff;
  }

  for (int var = 0; var < m_numVectorVars; var++)
  {
    m_makeFluxRegForVectorVar[var] = false;

  }

  m_makeFluxRegForVectorVar[VectorVars::m_fluidVel] = true;
  m_makeFluxRegForVectorVar[VectorVars::m_U_porosity] = true;

  defineIBCs ();



}

void AMRLevelMushyLayer::initialDataDefault()
{

  Real H =  m_parameters.bcValEnthalpyLo[SpaceDim-1];
  Real C =  m_parameters.bcValBulkConcentrationLo[SpaceDim-1];

  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    (*m_scalarNew[ScalarVars::m_enthalpy])[dit].setVal(H);
    (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit].setVal(C);

  }

}

void AMRLevelMushyLayer::initialDataSidewallHeating()
{
  // heating in x direction, i.e. left/right walls;
  int dir = 0;
  Real alpha = 0.1;

  Real deltaH =   m_parameters.bcValEnthalpyHi[dir] - m_parameters.bcValEnthalpyLo[dir];

  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {

    BoxIterator bit((*m_scalarNew[0])[dit].box());

    for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      //Real x_dir = loc[dir];
      //Real y = loc[1];


      Real perturbation = alpha*cos(loc[1]*M_PI)*sin(loc[0]*M_PI);


      (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) =  m_parameters.bcValEnthalpyLo[dir] + deltaH * loc[dir] + perturbation;
    }
  }

}

void AMRLevelMushyLayer::initialDataHRL()
{

  Real dir = 1;

  Real deltaH = m_parameters.bcValEnthalpyHi[dir] - m_parameters.bcValEnthalpyLo[dir];
  Real deltaC = m_parameters.bcValBulkConcentrationHi[dir] - m_parameters.bcValBulkConcentrationLo[dir];



  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {

    BoxIterator bit((*m_scalarNew[0])[dit].box());

    for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      Real x = loc[0];
      Real y = loc[1];




      Real alpha = 0.3; //0.02
      Real perturbation = alpha*cos(x*M_PI)*sin((y/m_domainHeight)*M_PI);

      int numRands = 2*m_opt.domainWidth/m_dx;
      Vector<Real> rands(numRands);
      for (int i = 1; i < numRands; i++)
      {
        rands[i] = ((double) rand()/ (RAND_MAX));
      }


      // if we have periodic bcs, the perturbation must be continuous at the boundary
      if (m_problem_domain.isPeriodic(0))
      {
        // White noise initialization

        int N = 0.5*m_opt.domainWidth/m_dx;
        Real weight = alpha/N;

        perturbation = 0;
        for (int n = 1; n <= N; n++)
        {

          Real additionalWeight = 1*rands[n];
          Real phaseShift = rands[n];
          //                if (n == 2)
          //                {
          //                  additionalWeight = 5;
          //                }
          perturbation += additionalWeight*weight*cos(2*n*M_PI*(x+phaseShift)/(m_opt.domainWidth))*sin((y/m_domainHeight)*M_PI);
        }

        // Random noise
        //              Real r = ((double) rand()/ (RAND_MAX));
        //              perturbation = alpha * r * sin((y/domainHeight)*M_PI);

        // No perturbation
        //              perturbation = 0;


        // Preferred wavelength
        alpha = 0.1;
        int wavenumber = 2; // 2 for vertical profiles (matching 3d profile)
        wavenumber = 3; // 3 for Nu(Ra)
        wavenumber = 1;
        perturbation = alpha*cos(2*wavenumber*M_PI*x/(m_opt.domainWidth))*sin((y/m_domainHeight)*M_PI);

        //perturbation = alpha*cos(8*(x/domainWidth)*M_PI)*sin((y/domainHeight)*M_PI);
      }



      (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) = m_parameters.bcValEnthalpyLo[dir] + deltaH * loc[dir]/m_domainHeight + perturbation;
      (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit](iv) = m_parameters.bcValBulkConcentrationLo[dir] + deltaC * loc[dir]/m_domainHeight;



    }
  }

}


void AMRLevelMushyLayer::initialDataConvectionMixedPorous()
{

  // Now do data

  int dir=0;
  if (m_parameters.bcValEnthalpyHi[0] == m_parameters.bcValEnthalpyLo[0])
  {
    dir = 1;
  }

  Real Htop = m_parameters.bcValEnthalpyHi[dir];
  Real Hbottom =  m_parameters.bcValEnthalpyLo[dir];

  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {

    BoxIterator bit((*m_scalarNew[0])[dit].box());

    for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      Real x = loc[0];
      Real y = loc[1];

      Real perturbation = 0;

      Real arg = 2*m_opt.perturbationWavenumber*x*M_PI/m_opt.domainWidth;

      if (m_problem_domain.isPeriodic(0))
      {
        //          perturbation = alpha*cos(2*x*M_PI/m_domainWidth)*sin(y*M_PI);
        perturbation = m_opt.initialPerturbation*sin(y*M_PI/m_domainHeight);

        if (m_opt.perturbationSin)
        {
          perturbation = perturbation*sin(arg);
        }
        else
        {
          perturbation = perturbation*cos(arg);
        }
      }
      else
      {
        //          Real arg = 2*pertWavenumber*x*M_PI/m_domainWidth;
        perturbation = -m_opt.initialPerturbation*sin(y*M_PI);

        if (m_opt.perturbationSin)
        {
          perturbation = perturbation*sin(arg);
        }
        else
        {
          perturbation = perturbation*cos(arg);
        }

      }
      // Make sure we don't take H above H max
      perturbation -= m_opt.initialPerturbation;


      (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) =Hbottom  +  (Htop - Hbottom) * loc[dir]/m_opt.domainWidth + perturbation;
      (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit](iv) = m_parameters.bcValBulkConcentrationLo[dir];

      if (!m_opt.doScalarAdvectionDiffusion)
      {
        //        (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) = (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) - 0.1*sin(M_PI*x/m_domainWidth)*cos(M_PI*y/m_domainHeight);
        //        (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) = Htop;
      }

      //      for (int dir=0;dir<SpaceDim;dir++)
      //      {
      //        (*m_vectorNew[VectorVars::m_fluidVel])[dit](iv, dir) = 0.0;
      //      }

      (*m_vectorNew[VectorVars::m_fluidVel])[dit](iv, 0) = m_opt.initVel*cos(M_PI*y/m_domainHeight)*sin(M_PI*x/m_opt.domainWidth);
      (*m_vectorNew[VectorVars::m_fluidVel])[dit](iv, 1) = -m_opt.initVel*cos(M_PI*x/m_opt.domainWidth)*sin(M_PI*y/m_domainHeight);

    }
  }
}

void AMRLevelMushyLayer::initialDataRayleighBenard()
{


  int numRands = 2*m_opt.domainWidth/m_dx;
  Vector<Real> rands(numRands);
  for (int i = 1; i < numRands; i++)
  {
    rands[i] = ((double) rand()/ (RAND_MAX));
  }

  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {

    BoxIterator bit((*m_scalarNew[0])[dit].box());

    for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      Real x = loc[0];
      Real y = loc[1];

      Real perturbation = m_opt.initialPerturbation*cos(x*M_PI)*sin((y/m_domainHeight)*M_PI);

      // if we have periodic bcs, the perturbation must be continuous at the boundary
      if (m_problem_domain.isPeriodic(0))
      {
        // White noise initialization

        int N = 0.5*m_opt.domainWidth/m_dx;
        Real weight = m_opt.initialPerturbation/N;

        perturbation = 0;
        for (int n = 1; n <= N; n++)
        {

          Real additionalWeight = 1*rands[n];
          Real phaseShift = rands[n];
          //                if (n == 2)
          //                {
          //                  additionalWeight = 5;
          //                }
          perturbation += additionalWeight*weight*cos(2*n*M_PI*(x+phaseShift)/(m_opt.domainWidth))*sin((y/m_domainHeight)*M_PI);
        }

        // Random noise
        //              Real r = ((double) rand()/ (RAND_MAX));
        //              perturbation = alpha * r * sin((y/domainHeight)*M_PI);

        // No perturbation
        //              perturbation = 0;


        // Preferred wavelength
        //        alpha = 0.1;
        //        int wavenumber = 2; // 2 for vertical profiles (matching 3d profile)
        //        wavenumber = 3; // 3 for Nu(Ra)

        perturbation = m_opt.initialPerturbation*cos(2*m_opt.perturbationWavenumber*M_PI*x/(m_opt.domainWidth))*sin((y/m_domainHeight)*M_PI);

        //perturbation = alpha*cos(8*(x/domainWidth)*M_PI)*sin((y/domainHeight)*M_PI);
      }

      (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) = m_parameters.bcValEnthalpyLo[1] + (m_parameters.bcValEnthalpyHi[1] - m_parameters.bcValEnthalpyLo[1]) * y/m_domainHeight + perturbation;
      //            (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit](iv) = ;


    }
  }

}

void AMRLevelMushyLayer::initialDataSoluteFlux()
{

  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {

    BoxIterator bit((*m_scalarNew[0])[dit].box());

    for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      //      Real x = loc[0];
      Real y = loc[1];


      //      (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv)  = (m_parameters.Hinitial - 0.2*m_parameters.stefan*exp((y-m_domainHeight)/0.2))*(1-0.05*sin(5*M_PI * x / m_domainWidth));
      (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv)  = m_parameters.Hinitial;

      //        valNew = m_parameters.ThetaInitial + sin(M_PI * y / m_domainHeight)*(1-0.05*cos(5*M_PI * x / m_domainWidth));
      (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit](iv) = -1 + y;

      // Shift y values because flux will be face centred

      Real yFace = y - 0.5*m_dx;
      Real magU = 10;
      Real freq = M_PI/m_opt.domainWidth;
      Real w = -magU*((1/freq)*(1-cos(freq)-0.5))*sin(freq*yFace); //horizontally averaged w
      //w = 1.0;
      w = -magU*(1/freq)*(cos(freq) - 1.0 + 0.5*freq)*exp(-freq*yFace);
      //    w = 0;
      (*m_scalarNew[ScalarVars::m_soluteFluxAnalytic])[dit](iv) = m_parameters.m_saltDiffusionCoeff + (1-yFace)*(m_parameters.nonDimVel + w);


    }
  }
}

void AMRLevelMushyLayer::initialDataBurgers()
{
  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {

    BoxIterator bit((*m_scalarNew[0])[dit].box());

    for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      Real x = loc[0];
      Real y = loc[1];

      for (int dir = 0; dir<SpaceDim; dir++)
      {

        if (m_problem_domain.isPeriodic(1))
        {
          burgersPeriodicInit( x,  y,  dir,  m_parameters);
        }
        else
        {
          burgersSinInit( x,  y,  dir,  m_parameters);
        }


      }

    }
  }

}

void AMRLevelMushyLayer::initialDataRefluxTest()
{
  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {

    BoxIterator bit((*m_scalarNew[0])[dit].box());

    for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      //      Real x = loc[0];
      Real y = loc[1];


      (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit](iv)  = -1.0;


      (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) = m_parameters.bcValEnthalpyLo[1] + (m_parameters.bcValEnthalpyHi[1] - m_parameters.bcValEnthalpyLo[1])*y;

    }

  }
}

void AMRLevelMushyLayer::addVortex(RealVect center, Real strength, Real radius)
{
  DataIterator dit = m_grids.dataIterator();

  Real r = 0;
  Real theta = 0;
  Real u_theta;

  for (dit.reset(); dit.ok(); ++dit)
  {
    FArrayBox& vel = (*m_vectorNew[VectorVars::m_fluidVel])[dit];
    Box b = vel.box();
    for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      Real x = loc[0];
      Real y = loc[1];

      Real x_rel = x-center[0];

      r = sqrt(pow(x-center[0], 2) + pow(y-center[1], 2));
      Real y_x = (y-center[1])/(x-center[0]);
      //        if ((y-center[1]) > 0)
      //        {
      //          theta = atan(y_x);
      //        }
      //        else
      //        {
      //          theta = -atan(y_x);
      //        }
      theta = atan(y_x);

      // Compute velocity in polar coords
      if (r < radius)
      {
        u_theta = strength*((8/3)*pow(r, 4)/pow(radius, 5) - 5*pow(r, 3)/pow(radius, 4) + (10/3)*r/pow(radius, 2));
      }
      else
      {
        u_theta = strength/r;
      }

      // Transform to cartesian coords (dr/dt = 0)
      Real sign = x_rel/abs(x_rel);

      vel(iv, 0) += -sign*r*u_theta*sin(theta);
      vel(iv, 1) += sign*r*u_theta*cos(theta);


    }
  }
}

void AMRLevelMushyLayer::initialDataVortexPair()
{

  setValLevel(*m_vectorNew[VectorVars::m_fluidVel], 0.0);

  // Vortex 1:
  RealVect vortexCenter;
  Real radius = 0.15;
  Real vortexStrength = 0.35;
  vortexCenter[0] = 0.3;
  vortexCenter[1] = 0.65;


  addVortex(vortexCenter, vortexStrength, radius);

  // Second vortex
  vortexStrength = -0.35;
  vortexCenter[0] = 0.3;
  vortexCenter[1] = 0.35;
  //

  //  vortexStrength = 1.0;
  //  vortexCenter[0] = 0.5;
  //  vortexCenter[1] = 0.5;

  addVortex(vortexCenter, vortexStrength, radius);


}

void AMRLevelMushyLayer::initialDataDiffusion()
{
  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    (*m_scalarNew[ScalarVars::m_enthalpy])[dit].setVal(m_parameters.Hinitial);

    (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit].setVal(m_parameters.bcValBulkConcentrationLo[1]);
  }
}


void AMRLevelMushyLayer::initialDataZeroPorosityTest()
{

  //  Real deltaC = -(m_parameters.compositionRatio+1);

  Real Cav = -m_parameters.compositionRatio/2;


  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {

    BoxIterator bit((*m_scalarNew[0])[dit].box());

    for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      Real x = loc[0];
      Real y = loc[1];

      (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) = m_parameters.stefan*(1+1.5*sin(M_PI*(x/m_opt.domainWidth*2 + 0.5))*sin(M_PI*y/m_domainHeight)); //(m_parameters.Hinitial - (m_parameters.stefan+10)*exp((y-m_domainHeight)/0.2));;
      (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit](iv) = Cav + Cav*sin(M_PI*x/m_opt.domainWidth*2)*sin(M_PI*y/m_domainHeight);


    }
  }
}

void AMRLevelMushyLayer::initialDataPorousHole()
{

  DataIterator dit = m_grids.dataIterator();

  Real HTop = m_parameters.bcValEnthalpyHi[1];
  Real HBottom = m_parameters.bcValEnthalpyLo[1];
  Real ThetaTop = m_parameters.bcValBulkConcentrationHi[1];
  Real ThetaBottom = m_parameters.bcValBulkConcentrationLo[1];

  Real thetaTop = m_parameters.bcValTemperatureHi[1];
  Real ThetaLTop = m_parameters.bcValLiquidConcentrationHi[1];

  Real cx = m_opt.domainWidth/2;
  Real cy = m_domainHeight/2;
  Real cz = m_domainHeight/2;



  for (dit.reset(); dit.ok(); ++dit)
  {
    BoxIterator bit((*m_scalarNew[0])[dit].box());
    for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);

      Real distanceFactor = 0;
      if (SpaceDim == 2)
      {
        distanceFactor = exp(- ( pow(loc[0]-cx,2) + pow(loc[1]-cy, 2))/(pow(m_opt.porousHoleRadius,2)));
      }
      else if(SpaceDim == 3)
      {
        distanceFactor = exp(- ( pow(loc[0]-cx,2) + pow(loc[1]-cy, 2) + pow(loc[2]-cz, 2))/(pow(m_opt.porousHoleRadius,2)));
      }
      (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) = HTop + (HBottom-HTop)*distanceFactor;
      (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit](iv) = ThetaTop + (ThetaBottom-ThetaTop)*distanceFactor;

      (*m_vectorNew[VectorVars::m_fluidVel])[dit](iv, SpaceDim-1) = m_opt.initVelScale*sin(M_PI*loc[0]/m_opt.domainWidth);
      (*m_vectorNew[VectorVars::m_fluidVel])[dit](iv, 0) = 0;

      (*m_vectorNew[VectorVars::m_bodyForce])[dit](iv, SpaceDim-1) =  m_parameters.m_buoyancySCoeff*ThetaLTop \
          - m_parameters.m_buoyancyTCoeff*thetaTop ;

    }
  }


}

void AMRLevelMushyLayer::initialDataMushyLayer()
{


  Real perturbation;

  DataIterator dit = m_grids.dataIterator();

  Real HTop = m_parameters.bcValEnthalpyHi[1];
  Real HBottom = m_parameters.bcValEnthalpyLo[1];
  Real ThetaBottom = m_parameters.bcValBulkConcentrationLo[1];
  Real thetaBottom = m_parameters.bcValTemperatureLo[1];

  if (m_opt.summerProfile > 0)
  {
    Real blScale = 0.00001;

    if (m_opt.summerProfile == 1)
    {
      blScale = 0.01*m_domainHeight;
    }

    for (dit.reset(); dit.ok(); ++dit)
    {
      BoxIterator bit((*m_scalarNew[0])[dit].box());
      for (bit.reset(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, loc, m_dx);
        Real y = loc[1];

        Real chiMin = 0.2;
        Real chiAv = (1+chiMin)/2;
        Real chi = chiAv - (1-chiAv)*tanh((y-m_opt.mushHeight)/(blScale));

        Real thetaMush = 1.01;
        Real thetaAv = (thetaMush + thetaBottom)/2;
        Real theta = thetaAv - (thetaBottom-thetaAv)*tanh((y-m_opt.mushHeight)/(blScale));

        Real H = m_parameters.stefan*chi + theta;

        Real S = 0;

        Real sice = -1.12; Real Sav = (sice + ThetaBottom)/2;
        S = Sav - (Sav-sice)*tanh((y-(m_opt.mushHeight*1.0))/(blScale));

        (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) = H;
        (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit](iv) = S;

      }
    }


  }
  else
  {
    // Not summer profile

    for (dit.reset(); dit.ok(); ++dit)
    {

      BoxIterator bit((*m_scalarNew[0])[dit].box());

      for (bit.reset(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, loc, m_dx);
        Real x = loc[0];
        Real y = loc[1];
        Real Hmin = 1.5;
        Real Hav = (Hmin+HBottom )/2;

        Real H = 0;

        if ((m_domainHeight - y) < m_opt.mushHeight/5)
        {
          H = Hmin + (y - (m_domainHeight -  m_opt.mushHeight/5))*((HTop-Hmin)/(m_opt.mushHeight/5));
        }
        else
        {
          H = Hav - (HBottom  - Hav)*tanh((y-m_opt.mushHeight)/(0.02));
        }

        Real arg = 2*m_opt.perturbationWavenumber*x*M_PI/m_opt.domainWidth;

        if (m_opt.linearGradient)
        {

          H = m_parameters.bcValEnthalpyLo[1] + (m_parameters.bcValEnthalpyHi[1] - m_parameters.bcValEnthalpyLo[1]) * (y/m_domainHeight);
        }
        else
        {
          H = m_parameters.Hinitial;
        }

        // Now compute perturbation
        // Only add perturbation to liquid region
        if (H > m_parameters.stefan)
        {

          if (m_problem_domain.isPeriodic(0))
          {
            //          perturbation = alpha*cos(2*x*M_PI/m_domainWidth)*sin(y*M_PI);
            perturbation = m_opt.initialPerturbation*sin(y*M_PI/m_domainHeight);

            if (m_opt.perturbationSin)
            {
              perturbation = perturbation*sin(arg);
            }
            else
            {
              perturbation = perturbation*cos(arg);
            }
          }
          else
          {
            perturbation = -m_opt.initialPerturbation*sin(y*M_PI);

            if (m_opt.perturbationSin)
            {
              perturbation = perturbation*sin(arg);
            }
            else
            {
              perturbation = perturbation*cos(arg);
            }

          }
          // Make sure we don't take H above H max
          perturbation -= m_opt.initialPerturbation;
        }
        else
        {
          perturbation = 0;
        }

        H = H + perturbation;


        (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) = H;

        if (m_opt.initVel)
        {
          (*m_vectorNew[VectorVars::m_fluidVel])[dit](iv, 0) = sin(2*M_PI*y);
          (*m_vectorNew[VectorVars::m_fluidVel])[dit](iv, 1) = sin(2*M_PI*x);
        }

      }
    }

  }


}

void AMRLevelMushyLayer::initialDataIceBlock()
{

  MayDay::Error("AMRLevelMushyLayer:: initialDataIceBlock Not implemented");

}


void AMRLevelMushyLayer::initialDataPoiseuille()
{
  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {

    BoxIterator bit((*m_scalarNew[0])[dit].box());

    for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      Real x = loc[0];
      Real y = loc[1];

      for (int dir=0; dir <SpaceDim; dir++)
      {
        stokesDarcyInit(x, y, dir, m_parameters);
      }


      if (m_opt.porosityFunction == PorosityFunctions::constantLiquid)
      {
        (*m_scalarNew[ScalarVars::m_porosity])[dit](iv) = 1.0;
      }
      else if (m_opt.porosityFunction == PorosityFunctions::linear)
      {
        (*m_scalarNew[ScalarVars::m_porosity])[dit](iv) = 0.1 + 0.9*x;
      }
      else if (m_opt.porosityFunction == PorosityFunctions::gaussian)
      {
        (*m_scalarNew[ScalarVars::m_porosity])[dit](iv) =  0.1+0.9*exp(-(x-0.5)*(x-0.5)*100);
      }
      else if (m_opt.porosityFunction == PorosityFunctions::gaussianSinusoidal)
      {
        (*m_scalarNew[ScalarVars::m_porosity])[dit](iv) =(exp(-(x-0.5)*(x-0.5)))*(0.5 + 0.1*sin(2*M_PI*y));
      }
      else if (m_opt.porosityFunction == PorosityFunctions::constantSmall)
      {
        (*m_scalarNew[ScalarVars::m_porosity])[dit](iv) = 0.1;
      }
      else
      {
        MayDay::Error("AMRLevelMushyLayer::inititalData - No porosity function specified");
      }

    }
  } // end if porosity

}

void AMRLevelMushyLayer::initialDataCornerFlow()
{

  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {

    BoxIterator bit((*m_scalarNew[0])[dit].box());

    for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      Real x = loc[0];
      Real y = loc[1];

      (*m_scalarNew[ScalarVars::m_porosity])[dit](iv) = (m_parameters.Hinitial - (m_parameters.stefan+10)*exp((y-m_domainHeight)/0.2));;
      (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit](iv) = m_parameters.ThetaInitial + (y-1)*exp(-(y-1)*(y-1)/0.01);

      //      Real fixed_porosity = -1.0;


      if (m_opt.fixedPorosity < 0)
      {

        if (m_opt.porosityFunction ==  PorosityFunctions::constantLiquid)
        {
          (*m_scalarNew[ScalarVars::m_porosity])[dit](iv)  = 1.0;
        }
        else if (m_opt.porosityFunction ==  PorosityFunctions::linear)
        {
          (*m_scalarNew[ScalarVars::m_porosity])[dit](iv)  = 0.1 + 0.9*y;
        }
        else if (m_opt.porosityFunction ==  PorosityFunctions::cubic)
        {
          (*m_scalarNew[ScalarVars::m_porosity])[dit](iv)  = 0.01 + pow(y-0.01, 3);
        }
        else if (m_opt.porosityFunction == PorosityFunctions::hyperbolicTan)
        {
          (*m_scalarNew[ScalarVars::m_porosity])[dit](iv)  = 0.5*( (1-tanh(20*y)) + (1-x));
        }
        else
        {
          MayDay::Error("No valid porosity function specified");
        }
      }
      else
      {
        (*m_scalarNew[ScalarVars::m_porosity])[dit](iv)  = m_opt.fixedPorosity;
      }


      Real f = - 1;
      Real fprime = 0;


      (*m_vectorNew[VectorVars::m_fluidVel])[dit](iv, 0)  = -x*fprime;
      (*m_vectorNew[VectorVars::m_fluidVel])[dit](iv, 1)  = f;


    }
  }

}

/*******/
void AMRLevelMushyLayer::initialData()
{
  pout() << "AMRLevelMushyLayer::initialData - setting initial data on level " << m_level << endl;

  // For some reason AMR calls initialData() on levels
  // which don't have any grids yet. Can't fill empty grids
  // so just do nothing in this case.
  if (m_level_grids.size() == 0)
  {
    return;
  }

  // Defaults:
  DataIterator dit = m_grids.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {

    (*m_scalarNew[ScalarVars::m_lambda])[dit].setVal(1.0);
    (*m_scalarOld[ScalarVars::m_lambda])[dit].setVal(1.0);
    //    (*m_dScalar[ScalarVars::m_lambda])[dit].setVal(1.0);

    (*m_scalarNew[ScalarVars::m_enthalpy])[dit].setVal(m_parameters.Hinitial);
    (*m_scalarOld[ScalarVars::m_enthalpy])[dit].setVal(m_parameters.Hinitial);
    //    (*m_dScalar[ScalarVars::m_enthalpy])[dit].setVal(m_parameters.Hinitial);

    (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit].setVal(m_parameters.ThetaInitial);
    (*m_scalarOld[ScalarVars::m_bulkConcentration])[dit].setVal(m_parameters.ThetaInitial);
    //    (*m_dScalar[ScalarVars::m_bulkConcentration])[dit].setVal(m_parameters.ThetaInitial);

    (*m_vectorNew[VectorVars::m_fluidVel])[dit].setVal(0.0);
    (*m_vectorOld[VectorVars::m_fluidVel])[dit].setVal(0.0);
    (*m_dVector[VectorVars::m_fluidVel])[dit].setVal(0.0);

    (*m_vectorNew[VectorVars::m_bodyForce])[dit].setVal(m_parameters.body_force);
    (*m_vectorOld[VectorVars::m_bodyForce])[dit].setVal(m_parameters.body_force);
    (*m_dVector[VectorVars::m_bodyForce])[dit].setVal(m_parameters.body_force);

    m_advVel[dit].setVal(0.0);


  }


  // If we've specified some custom initial data, deal with it here
  if (m_opt.customInitData >= 0)
  {
    switch(m_opt.customInitData)
    {
      case 1:
        initialDataPorousHole();
        break;
      default:
        initialDataDefault();
    }

  }
  else
  {

    // May want to overwrite default options for some problems
    switch(m_parameters.physicalProblem)
    {
      case PhysicalProblems::m_sidewallHeating:
        initialDataSidewallHeating();
        break;
      case PhysicalProblems::m_HRL:
        initialDataHRL();
        break;
      case PhysicalProblems::m_convectionMixedPorous:
        initialDataConvectionMixedPorous();
        break;
      case PhysicalProblems::m_rayleighBenard:
        initialDataRayleighBenard();
        break;
      case PhysicalProblems::m_cornerFlow:
        initialDataCornerFlow();
        break;
      case  PhysicalProblems::m_mushyLayer:
        initialDataMushyLayer();
        break;
      case PhysicalProblems::m_poiseuilleFlow:
        initialDataPoiseuille();
        break;
      case PhysicalProblems::m_soluteFluxTest:
        initialDataSoluteFlux();
        break;
      case PhysicalProblems::m_refluxTest:
        initialDataRefluxTest();
        break;
      case PhysicalProblems::m_zeroPorosityTest:
        initialDataZeroPorosityTest();
        break;
      case PhysicalProblems::m_meltingIceBlock:
        initialDataIceBlock();
        break;
      case PhysicalProblems::m_diffusion:
        initialDataDiffusion();
        break;
      case PhysicalProblems::m_vortexPair:
        initialDataVortexPair();
        break;


      default:
        initialDataDefault();
    }

  }
  addMeltPond();

//  if(m_parameters.physicalProblem == PhysicalProblems::m_poiseuilleFlow)
//  {
//    stokesDarcyForcing(*m_scalarNew[ScalarVars::m_temperature], 0);
//  }
//  else if (m_parameters.physicalProblem == PhysicalProblems::m_cornerFlow)
//  {
//    // Don't do anything
//  }
//  else
//  {
    updateEnthalpyVariables();
    copyNewToOldStates();
    updateEnthalpyVariables();
//  }

  calculatePermeability(); //make sure this is up to date

  // need this to get the initial timestep right
  //Initialise advVel
  CellToEdge(*m_vectorNew[VectorVars::m_fluidVel], m_advVel);
  CellToEdge(*m_vectorNew[VectorVars::m_fluidVel], m_advVelOld);
  CellToEdge(*m_vectorNew[VectorVars::m_fluidVel], m_advVelNew);

  m_advVel.exchange();
  m_frameAdvVel.exchange();

  calculateAnalyticSolns(m_opt.enforceAnalyticSoln);

  if (m_opt.initAnalyticVel)
  {
    fillAnalyticVel(m_advVel);

    EdgeToCell(m_advVel, *m_vectorNew[VectorVars::m_fluidVel]);
  }



}

void AMRLevelMushyLayer::addPerturbation(int a_var, Real alpha, int waveNumber, Real phaseShift)
{
  // ignore small perturbations
  if (alpha < 1e-20)
  {
    return;
  }

  pout() << "Adding perturbation " << alpha << " with wavenumber " << waveNumber << endl;

  Real domainWidth = m_opt.domainWidth;


  Vector<Real> rands1(m_numCells[0]);
  Vector<Real> rands2(m_numCells[1]);
  Vector<Real> randsk(m_opt.maxRestartWavenumbers);
  for (int i = 0; i < m_numCells[0]; i++)
  {
    rands1[i] = ((double) rand()/ (RAND_MAX));
  }
  for (int i = 0; i < m_numCells[1]; i++)
  {
    rands2[i] = ((double) rand()/ (RAND_MAX));
  }
  for (int i = 0; i < m_opt.maxRestartWavenumbers; i++)
  {
    randsk[i] = ((double) rand()/ (RAND_MAX));
  }

  for (DataIterator dit = m_scalarNew[a_var]->dataIterator(); dit.ok(); ++dit)
  {
    Box b = (*m_scalarNew[a_var])[dit].box();
    b &= m_problem_domain.domainBox();

    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc;
      getLocation(iv, loc, m_dx);
      Real x = loc[0];
      Real y = loc[1];

      Real perturbation = 0.0;
      //Real Nwaves = round(domainWidth);
      //Nwaves = Nwaves*10;
      if (waveNumber > 0)
      {
        // sine wave perturbation
        Real wavelength = m_opt.domainWidth/waveNumber;
        if (m_problem_domain.isPeriodic(0))
        {
          //                         perturbation = alpha*cos(2*x*M_PI/domainWidth)*sin(y*M_PI);
          perturbation = alpha*cos(M_PI*(phaseShift + 2*x/wavelength))*cos(2*y*M_PI/wavelength);
        }
        else
        {
          perturbation = -alpha*cos(M_PI*(phaseShift + 2*x/domainWidth))*sin(y*M_PI/wavelength);
        }

      }
      else if (waveNumber == 0)
      {
        // Sum of wavenumbers 1:50 in x direction
        // random in z direction, going to zero at top and bottom
        for (int k=1; k<m_opt.maxRestartWavenumbers; k++)
        {
          Real random = ((double) rand()/ (RAND_MAX));

          Real thisPert = (alpha/m_opt.maxRestartWavenumbers)*cos(M_PI*(k*x/domainWidth + randsk[k])); // x component

          thisPert = thisPert * random * sin(y*M_PI/m_domainHeight); // z component
          perturbation += thisPert;
        }

      }
      else
      {
        // Random perturbation if no wave number specified
        Real random = ((double) rand()/ (RAND_MAX));
        perturbation = alpha*pow(random, 2);
      }

      (*m_scalarNew[a_var])[dit](iv) += perturbation;
    }

  }
}

void AMRLevelMushyLayer::fillFrameVelocity()
{
  for (DataIterator dit = m_frameAdvVel.dataIterator(); dit.ok(); ++dit)
  {
    m_frameAdvVel[dit][SpaceDim-1].setVal(m_parameters.nonDimVel);
    //    m_frameVel[dit].setVal(m_parameters.nonDimVel, SpaceDim-1);
    for (int dir = 0; dir < SpaceDim-1; dir++)
    {
      m_frameAdvVel[dit][dir].setVal(0);
      //      m_frameVel[dit].setVal(0, dir);
    }

  }
}

void AMRLevelMushyLayer::calculateAnalyticSolns(bool enforceSolutions)
{
  CH_TIME("AMRLevelMushyLayer::calculateAnalyticSolns");

  // Calculate some analytic solutions if appropriate
  IntVect ivGhost = IntVect::Unit*m_numGhost;

  LevelData<FArrayBox> enthalpyAnalytic(m_grids, 1, ivGhost);
  LevelData<FArrayBox> bulkCAnalytic(m_grids, 1, ivGhost);


  if (m_opt.analyticSolution == PhysicalProblems::m_solidificationNoFlow)
  {

    LevelData<FArrayBox> CLanalytic(m_grids, 1, ivGhost);
    LevelData<FArrayBox> CSanalytic(m_grids, 1, ivGhost);

    ::analyticSolnSolidificationNoFlow(enthalpyAnalytic, bulkCAnalytic,
                                       m_domainHeight, m_dx,
                                       m_parameters);

    ::updateEnthalpyVariables(enthalpyAnalytic, bulkCAnalytic,
                              *m_scalarNew[ScalarVars::m_temperatureAnalytic], CLanalytic,
                              CSanalytic, *m_scalarNew[ScalarVars::m_porosityAnalytic],
                              m_parameters);

    ::updateEnthalpyVariables(enthalpyAnalytic, bulkCAnalytic,
                              *m_scalarOld[ScalarVars::m_temperatureAnalytic], CLanalytic,
                              CSanalytic, *m_scalarOld[ScalarVars::m_porosityAnalytic],
                              m_parameters);

    if (enforceSolutions)
    {


      for (DataIterator dit = m_scalarOld[ScalarVars::m_porosityAnalytic]->dataIterator(); dit.ok(); ++dit)
      {
        (*m_scalarOld[ScalarVars::m_enthalpy])[dit].setVal(0);
        (*m_scalarOld[ScalarVars::m_enthalpy])[dit] += enthalpyAnalytic[dit];
        (*m_scalarNew[ScalarVars::m_enthalpy])[dit].copy((*m_scalarOld[ScalarVars::m_enthalpy])[dit]);

        (*m_scalarOld[ScalarVars::m_bulkConcentration])[dit].setVal(0);
        (*m_scalarOld[ScalarVars::m_bulkConcentration])[dit] += bulkCAnalytic[dit];
        (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit].copy((*m_scalarOld[ScalarVars::m_bulkConcentration])[dit]);

        if (m_opt.initialPerturbation != 0.0)
        {
          Box b = (*m_scalarNew[ScalarVars::m_enthalpy])[dit].box();
          if (m_opt.perturbationWavenumber > -1)
          {
            for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc;
              ::getLocation(iv, loc, m_dx);

              Real pert;
              if (m_opt.perturbationSin)
              {
                pert = m_opt.initialPerturbation*sin(M_PI*2*(loc[1]/m_domainHeight)*m_opt.perturbationWavenumber)*sin(M_PI*(loc[0]/m_opt.domainWidth)*m_opt.perturbationWavenumber*2);
              }
              else
              {
                pert = m_opt.initialPerturbation*sin(M_PI*2*(loc[1]/m_domainHeight)*m_opt.perturbationWavenumber)*cos(M_PI*(loc[0]/m_opt.domainWidth)*m_opt.perturbationWavenumber*2);
              }

              pert *= (*m_scalarNew[ScalarVars::m_porosity])[dit](iv);
              (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) += pert;
            }
          }
          else
          {

            // Add a small random perturbation everywhere

            for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              Real r = ((double) rand()/ (RAND_MAX));
              (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) += r*m_opt.initialPerturbation;
            }

          }
        }
      }
      updateEnthalpyVariablesOld();
      updateEnthalpyVariables();
    }
  }
  else if (m_opt.analyticSolution == PhysicalProblems::m_poiseuilleFlow)
  {

    for (DataIterator dit = m_vectorNew[VectorVars::m_fluidVel]->dataIterator(); dit.ok(); ++dit)
    {
      Box b = (*m_vectorNew[VectorVars::m_fluidVel])[dit].box();
      for (BoxIterator bit(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, loc, m_dx);
        Real x = loc[0];
        //				Real y = loc[1];


        //        pp.query("fixed_porosity", fixedPorosity);
        (*m_vectorNew[VectorVars::m_fluidVelAnalytic])[dit](iv, 0) = 0;

        if(m_opt.fixedPorosity >= 0)
        {
          Real perm  = m_parameters.calculatePermeability(m_opt.fixedPorosity);
          Real L = sqrt(m_opt.fixedPorosity/(m_parameters.darcy*perm));
          (*m_vectorNew[VectorVars::m_fluidVelAnalytic])[dit](iv, 1) = m_parameters.rayleighTemp*perm*(1-(cosh(L*(0.5-x)))/(cosh(0.5*L)));
        }
        else
        {
          (*m_vectorNew[VectorVars::m_fluidVelAnalytic])[dit](iv, 1) = m_parameters.rayleighTemp*(1-4*(x-0.5)*(x-0.5))/8;
        }

        (*m_vectorOld[VectorVars::m_fluidVelAnalytic])[dit](iv, 0) = (*m_vectorNew[VectorVars::m_fluidVelAnalytic])[dit](iv, 0);
        (*m_vectorOld[VectorVars::m_fluidVelAnalytic])[dit](iv, 1) = (*m_vectorNew[VectorVars::m_fluidVelAnalytic])[dit](iv, 1);


        if (enforceSolutions)
        {
          (*m_vectorOld[VectorVars::m_fluidVel])[dit](iv, 0) = (*m_vectorNew[VectorVars::m_fluidVelAnalytic])[dit](iv, 0);
          (*m_vectorOld[VectorVars::m_fluidVel])[dit](iv, 1) = (*m_vectorNew[VectorVars::m_fluidVelAnalytic])[dit](iv, 1);

          (*m_vectorNew[VectorVars::m_fluidVel])[dit](iv, 0) = (*m_vectorNew[VectorVars::m_fluidVelAnalytic])[dit](iv, 0);
          (*m_vectorNew[VectorVars::m_fluidVel])[dit](iv, 1) = (*m_vectorNew[VectorVars::m_fluidVelAnalytic])[dit](iv, 1);
        }


        setPorosity((*m_scalarNew[ScalarVars::m_porosity])[dit]);

      }

    }
  }
  else if (m_opt.analyticSolution == PhysicalProblems::m_diffusion)
  {
    for (DataIterator dit = m_vectorNew[VectorVars::m_fluidVel]->dataIterator(); dit.ok(); ++dit)
    {
      Box b = (*m_vectorNew[VectorVars::m_fluidVel])[dit].box();
      for (BoxIterator bit(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, loc, m_dx);
        //        Real x = loc[0];
        Real y = loc[1];

        (*m_scalarNew[ScalarVars::m_temperatureAnalytic])[dit](iv) = m_parameters.bcValTemperatureLo[0] - 0.5*y*(y-1);
      }
    }
  }
  else if (m_opt.analyticSolution == PhysicalProblems::m_mushyLayer)
  {
    // when we ask for an analytic solution here, just generate a rough chimney shape

    ::channelFieldMushyLayer(enthalpyAnalytic, bulkCAnalytic,
                             m_domainHeight, m_opt.domainWidth, m_dx,
                             m_parameters);

    if (enforceSolutions)
    {

      for (DataIterator dit = m_scalarOld[ScalarVars::m_porosityAnalytic]->dataIterator(); dit.ok(); ++dit)
      {
        (*m_scalarOld[ScalarVars::m_enthalpy])[dit].setVal(0);
        (*m_scalarOld[ScalarVars::m_enthalpy])[dit] += enthalpyAnalytic[dit];
        (*m_scalarNew[ScalarVars::m_enthalpy])[dit].copy((*m_scalarOld[ScalarVars::m_enthalpy])[dit]);

        (*m_scalarOld[ScalarVars::m_bulkConcentration])[dit].setVal(0);
        (*m_scalarOld[ScalarVars::m_bulkConcentration])[dit] += bulkCAnalytic[dit];
        (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit].copy((*m_scalarOld[ScalarVars::m_bulkConcentration])[dit]);

      }
      updateEnthalpyVariablesOld();
      updateEnthalpyVariables();
    }

  }


}

void AMRLevelMushyLayer::setPorosity(FArrayBox& a_porosity)
{
  // Enforce correct porosity

  Real porosity = 1.0;

  for (BoxIterator bit(a_porosity.box()); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    //          porosity = fixedPorosity;

    RealVect loc;
    getLocation(iv, loc, m_dx);
    Real x = loc[0];
    Real y = loc[1];

    if (m_opt.fixedPorosity >= 0)
    {
      porosity = m_opt.fixedPorosity;
    }
    else
    {

      if (m_opt.porosityFunction == PorosityFunctions::constantLiquid)
      {
        porosity = 1.0;
      }
      else if (m_opt.porosityFunction == PorosityFunctions::linear)
      {
        porosity = 0.1 + 0.9*x;
      }
      else if (m_opt.porosityFunction == PorosityFunctions::gaussian)
      {
        porosity =  0.1+0.9*exp(-(x-0.5)*(x-0.5)*100);
      }
      else if (m_opt.porosityFunction == PorosityFunctions::gaussianSinusoidal)
      {
        porosity =(exp(-(x-0.5)*(x-0.5)))*(0.5 + 0.1*sin(2*M_PI*y));
      }
      else if (m_opt.porosityFunction == PorosityFunctions::constantSmall)
      {
        porosity = 0.1;
      }
      else
      {
        MayDay::Error("AMRLevelMushyLayer::analyticsolns - No porosity function specified");
      }

    }

    a_porosity(iv) = porosity;
    a_porosity(iv) = porosity;
  } // end loop over IVs

}



void AMRLevelMushyLayer::fillAMRVelPorosity(Vector<LevelData<FArrayBox>*> & amrVel,
                                            Vector<RefCountedPtr< LevelData<FluxBox> > > &  amrPorosityFace,
                                            Vector<RefCountedPtr< LevelData<FArrayBox> > >& amrPorosity)
{
  AMRLevelMushyLayer* thisLevelData = getCoarsestLevel();
  int base_level = 0;
  CH_assert(thisLevelData->m_level == base_level);

  while (thisLevelData->hasFinerLevel())
  {
    thisLevelData = thisLevelData->getFinerLevel();
  }
  //          CH_assert(thisLevelData->finestLevel());
  int numLevels = thisLevelData->m_level + 1;

  //    Vector<LevelData<FArrayBox>*> amrVel(numLevels);
  //    Vector<RefCountedPtr< LevelData<FluxBox> > > amrPorosityFace(numLevels);
  //    Vector<RefCountedPtr< LevelData<FArrayBox> > > amrPorosity(numLevels);


  thisLevelData = getCoarsestLevel();

  for (int lev = base_level; lev < numLevels; lev++)
  {
    amrVel[lev] = &(*thisLevelData->m_vectorNew[VectorVars::m_fluidVel]);
    // Fill CF BCs
    thisLevelData->fillVectorField(*amrVel[lev], m_time, m_fluidVel, false, true);

    amrPorosity[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(thisLevelData->m_grids, 1));
    thisLevelData->fillScalars(*amrPorosity[lev], m_time, m_porosity, true, true);

    amrPorosityFace[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(thisLevelData->m_grids, 1));
    thisLevelData->fillScalarFace(*amrPorosityFace[lev], m_time, m_porosity, true, true);

    if (!thisLevelData->finestLevel())
      thisLevelData =
          dynamic_cast<AMRLevelMushyLayer*>(thisLevelData->m_finer_level_ptr);
  }

}

/*******/
void AMRLevelMushyLayer::postInitialize()
{

  pout() << "AMRLevelMushyLayer::postInitialize on level " << m_level << endl;


  // initialize data structures which haven't yet been initialized
  //	levelSetup();
  defineSolvers(m_time-m_dt); // Used to be in postinitialgrid

  // set up for initial velocity projection -- do this from level 0,
  // since post_initialize is called from fine->coarse
  if (m_level == 0 && m_opt.doProjection)
  {
    AMRLevelMushyLayer* thisLevelData = this;
    while (thisLevelData->hasFinerLevel())
    {
      thisLevelData = thisLevelData->getFinerLevel();
    }
    //		CH_assert(thisLevelData->finestLevel());
    int numLevels = thisLevelData->m_level + 1;

    Vector<LevelData<FArrayBox>*> amrVel(numLevels);
    Vector<RefCountedPtr< LevelData<FluxBox> > > amrPorosityFace(numLevels);
    Vector<RefCountedPtr< LevelData<FArrayBox> > > amrPorosity(numLevels);

    fillAMRVelPorosity(amrVel, amrPorosityFace, amrPorosity);

    Projector& level0Proj = m_projection;

    // set physical boundary conditions on velocity

    thisLevelData = this;

    VelBCHolder velBC(m_physBCPtr->uStarFuncBC(m_opt.viscousBCs));

    for (int lev = 0; lev < numLevels; lev++)
    {
      const ProblemDomain& levelDomain = thisLevelData->problemDomain();
      Real levelDx = thisLevelData->m_dx;
      LevelData<FArrayBox>& thisAmrVel = *amrVel[lev];
      const DisjointBoxLayout& thisLevelGrids = thisAmrVel.getBoxes();
      velBC.applyBCs(thisAmrVel, thisLevelGrids, levelDomain, levelDx,
                     false); // inhomogeneous

      // This deals with periodic BCs
      thisAmrVel.exchange();

      if (thisLevelData->m_finer_level_ptr != NULL)
      {
        thisLevelData =
            dynamic_cast<AMRLevelMushyLayer*>(thisLevelData->m_finer_level_ptr);
      }
    }

    bool homoBC = false;
    if (m_opt.project_initial_vel)
    {
      level0Proj.initialVelocityProject(amrVel, amrPorosityFace, amrPorosity, homoBC);
    }

    // Compute initial volume discrepancy correction
    // For this we need to advect lambda on each level
    // to get a measure of initial freestream violation
    if (m_opt.compute_initial_VD_corr)
    {
      if (s_verbosity >= 3)
      {
        pout() << "AMRlevelMushyLayer::postInitialize - compute initial VD corr" << endl;
      }

      CH_TIME("compute_initial_VD_corr");

      Vector<LevelData<FArrayBox>*> amrLambda(numLevels);

      //      Vector<RefCountedPtr<LevelData<FluxBox> > > amrPorosity(numLevels);

      AMRLevelMushyLayer* thisML = this;

      // Do a calculation on the coarsest level to get an estimate of the velocity
      // we need a very rough estimate of dt in order to do this
      Real initDt = -1;
      Real initTime = -1;


      if (m_opt.fixedDt > 0)
      {
        initDt = m_opt.fixedDt;
      }
      else
      {
        initDt= thisML->computeInitialDt();
        //        initDt = 1e-4*initDt;
        thisML->dt(initDt);
        thisML->computeInitAdvectionVel();

        // Now get a better estimate of the initial dt
        initDt = thisML->computeInitialDt();
        initDt = 1e-4*initDt; //make it a lot smaller to be safe
      }

      initTime = m_time + initDt;

      pout() << "Computing initial VD correction with dt = " << initDt << ", time = " << m_time << endl;

      for (int lev = 0; lev < numLevels; lev++)
      {
        thisML->dt(initDt);

        thisML->time(initTime);

        // initialize pressures to first
        // guess (most likely 0)
        setValLevel(thisML->m_projection.Pi(), 0.0);

        thisML->computeInitAdvectionVel();

        // Advect lambda
        thisML->advectLambda(true);

        thisML = thisML->getFinerLevel();
      }

      AMRRefluxLambda();


      thisML = this;
      for (int lev = 0; lev < numLevels; lev++)
      {
        // Finally, stick lambda into amrLambda
        amrLambda[lev] = thisML->m_scalarNew[ScalarVars::m_lambda];
        //        amrPorosity[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(thisML->m_grids, 1));
        //        thisML->fillScalarFace(*amrPorosity[lev],m_time, m_porosity, true, true);
        //thisML->m_scalarNew[ScalarVars::m_porosity];
        thisML = thisML->getFinerLevel();
      }

      level0Proj.doInitialSyncOperations(amrVel, amrLambda,
                                         amrPorosityFace, amrPorosity,
                                         initTime, initDt); // false - don't add to grad_elambda, simply reset it

      // Reset original fields incl. lambda
      thisML = this;
      for (int lev = 0; lev < numLevels; lev++)
      {
        thisML->time(initTime - initDt);

        thisML->copyNewToOldStates();

        // Make sure lambda is definitely 1
        thisML->resetLambda();

        thisML = thisML->getFinerLevel();
      }
    }


    thisLevelData = this;

    // need to reset boundary conditions here
    for (int lev = 0; lev < numLevels; lev++)
    {
      const ProblemDomain& levelDomain = thisLevelData->problemDomain();
      Real levelDx = thisLevelData->m_dx;
      LevelData<FArrayBox>& thisAmrVel = *amrVel[lev];
      const DisjointBoxLayout& thisLevelGrids = thisAmrVel.getBoxes();
      velBC.applyBCs(thisAmrVel, thisLevelGrids, levelDomain, levelDx,
                     false); // inhomogeneous
      if (thisLevelData->m_finer_level_ptr != NULL)
      {
        thisLevelData =
            dynamic_cast<AMRLevelMushyLayer*>(thisLevelData->m_finer_level_ptr);
      }
    }

    if (m_opt.initialize_pressures)
    {
      if (s_verbosity >= 3)
      {
        pout() << "AMRlevelMushyLayer::postInitialize - initialize pressures" << endl;
      }

//      if (solvingFullDarcyBrinkman())
      if (m_opt.doEulerPart)
      {
        Real dtInit = computeDtInit(numLevels-1);

        dtInit *= m_opt.init_dt_scale;

        thisLevelData = this->getCoarsestLevel();
        for (int lev = 0; lev < numLevels; lev++)
        {
          thisLevelData->m_projection.setNonSubcycledMACBCs();
          thisLevelData = thisLevelData->getFinerLevel();
        }

        initializeGlobalPressure(dtInit, true);

        // Reset lambda
        thisLevelData = this->getCoarsestLevel();
        for (int lev = 0; lev < numLevels; lev++)
        {
          thisLevelData->resetLambda();

          thisLevelData->m_projection.setSubcycledMACBCs();

          thisLevelData = thisLevelData->getFinerLevel();
        }

      }
      else
      {
        // Compute advection vel on every level, and keep doing so until div u is small enough
        // Don't actually need to do this, as we don't need the previous pressure in order to compute
        // the velocity at a new timestep. It makes the solve a bit quicker though, so we'll do it for level 0.
        // Other levels then use this as a BC so really no need to init them at all

        AMRLevelMushyLayer* lev = this->getCoarsestLevel();
        initTimeIndependentPressure(lev);

      }

    }

  }

  // Calculate initial sum of salt
  if (m_level == 0)
  {
    Vector<LevelData<FArrayBox>* > S_new, H_new;

    Vector<AMRLevelMushyLayer*> mlVect;
    Vector<DisjointBoxLayout> grids;
    Vector<int> refRat;
    ProblemDomain lev0Dom;
    Real lev0Dx;
    getHierarchyAndGrids(mlVect, grids, refRat, lev0Dom, lev0Dx);

    int nLevels = mlVect.size();

    S_new.resize(nLevels);
    H_new.resize(nLevels);

    AMRLevelMushyLayer* ml = this;

    for (int lev = 0; lev < nLevels;  ++lev)
    {
      S_new[lev] = ml->m_scalarNew[ScalarVars::m_bulkConcentration];
      H_new[lev] = ml->m_scalarNew[ScalarVars::m_enthalpy];

      ml = ml->getFinerLevel();
    }

    // Calculate sums over valid regions
    AMRSaltSum_old = ::computeSum(S_new, refRat, lev0Dx, Interval(0,0), 0);
    AMREnthalpySum_old = ::computeSum(H_new, refRat, lev0Dx, Interval(0,0), 0);
  }

  // Store our initialized values in case the first timestep crashes (quite likely)
  backupTimestep();

  // Get all init plot fields (including pressure fields from the projection)
  if (m_level==0)
  {
    AMRLevelMushyLayer* thisML = this;
    while(thisML)
    {
      thisML->getExtraPlotFields();
      thisML = thisML->getFinerLevel();
    }
  }

}

void AMRLevelMushyLayer::initTimeIndependentPressure(AMRLevelMushyLayer* lev, int a_max_num_iter)
{
  Real maxDivU = 1e10;
  int i = 0;
  int maxNumIter = m_opt.num_init_passes;

  // Possible for the argument to override the default number of passes
  if (a_max_num_iter >= 0)
  {
    maxNumIter = a_max_num_iter;
  }

  while(maxDivU > 1e-10 && i < maxNumIter)
  {
    lev->calculateTimeIndAdvectionVel(lev->m_time, lev->m_advVel);

    Divergence::levelDivergenceMAC(*lev->m_scalarNew[ScalarVars::m_divUadv], lev->m_advVel, m_dx);
    maxDivU = ::computeNorm(*lev->m_scalarNew[ScalarVars::m_divUadv], NULL, 1, lev->m_dx, Interval(0,0), 0);
    pout() << "  Pressure init " << i << ", max(div U) = " << maxDivU << endl;

    i = i + 1;
  }
}

void AMRLevelMushyLayer::resetLambda()
{
  setValLevel(*m_scalarOld[ScalarVars::m_lambda], 1.0);
  setValLevel(*m_scalarNew[ScalarVars::m_lambda], 1.0);
}

void AMRLevelMushyLayer::AMRRefluxLambda()
{
  // Reflux lambda to capture freestream (non)-preservation
  // Should do this from finest to coarsest level
  AMRLevelMushyLayer* thisML = this;
  int finest_level = m_level;
  while(thisML->hasFinerLevel())
  {
    thisML = thisML->getFinerLevel();
    finest_level++;
  }

  // thisML is now the finest level
  // want the second finest though, as no reflux to be done on finest level
  thisML = thisML->getCoarserLevel();

  for (int lev = finest_level-1; lev >= 0; lev--)
  {
    thisML->doExplicitReflux(ScalarVars::m_lambda);
    thisML = thisML->getCoarserLevel();
  }
}

void AMRLevelMushyLayer::computeInitAdvectionVel()
{
  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::computeInitAdvectionVel (level " << m_level << ")" << endl;
  }

//  if (solvingFullDarcyBrinkman())
  if (m_opt.doEulerPart)
  {
    // Just initialising so don't advect/diffuse any scalars here

    IntVect ivGhost = m_numGhostAdvection*IntVect::Unit;
    LevelData<FArrayBox> advectionSourceTerm(m_grids, SpaceDim, ivGhost);

    computeAdvectionVelSourceTerm(advectionSourceTerm);

    if (!m_opt.doEulerPart)
    {
      // Compute unprojected cell centred velocity (don't need to project it though)

      computeCCvelocity(advectionSourceTerm, m_time, m_dt,
                        false, // don't do flux register updates
                        false, // don't do projection
                        false, // don't compute u.del(u)
                        false); // not a mac projection

      m_vectorNew[m_UpreProjection]->copyTo(*m_vectorOld[m_UpreProjection]);
    }

    computeAdvectionVelocities(advectionSourceTerm);

  }
  else
  {
    calculateTimeIndAdvectionVel(m_time, m_advVel);
  }
}

// -------------------------------------------------------------
// this function manages the pressure initialization after
// initialization and regridding
void AMRLevelMushyLayer::initializeGlobalPressure(Real dtInit, bool init)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::initializeGlobalPressure (level " << m_level << ")" << endl;
  }

  // This function is only for darcy-brinkman - throw an error if that's not the case
  CH_assert(solvingFullDarcyBrinkman());

  setAdvVelCentering(m_opt.initAdvVelCentering);

  // index through levels to find out what finest level is
  AMRLevelMushyLayer* thisMLPtr = this;
  while (!thisMLPtr->finestLevel())
  {
    thisMLPtr = thisMLPtr->getFinerLevel();
  }
  CH_assert(thisMLPtr->finestLevel());
  int finest_level = thisMLPtr->m_level;

  // first compute dtInit
  Real cur_time = m_time;

  // If we're starting the simulation, we may want to change the time
  if (init && m_opt.restart_new_time >= 0.0)
  {
    cur_time = m_opt.restart_new_time;
  }

  //  Real dtLevel;
  //  Real dtInit = 10000000.0;
  Real max_dtInit = 0.5*computeDtInit(finest_level);
  if (dtInit < 0)
  {
    dtInit = computeDtInit(finest_level);
    dtInit *= 0.5;
  }

  if (m_time == 0)
  {
    dtInit = dtInit/1000;
  }

  // Save dt for each level so it can be reset later
  Vector<Real> dtSave(finest_level + 1, 1.0e8);
  thisMLPtr = this;
  for (int lev = m_level; lev <= finest_level; lev++)
  {
    dtSave[lev] = thisMLPtr->m_dt;

    thisMLPtr = thisMLPtr->getFinerLevel();
  }

  thisMLPtr = this;
  bool orig_add_subtract_grad_p = thisMLPtr->m_usePrevPressureForUStar;

  for (int lev = m_level; lev <= finest_level; lev++)
  {

    // Option to turn off adding/subtracting grad(p) for initialisation
    if (!m_opt.init_use_prev_pressure_for_Ustar)
    {
      thisMLPtr->m_usePrevPressureForUStar=false;
    }

    thisMLPtr = thisMLPtr->getFinerLevel();

  }

  // now loop through levels and compute estimate of pressure
  // gradually increase dt during initialization steps

  thisMLPtr = this;
  int lbase = thisMLPtr->m_level;
  for (int iter = 0; iter < m_opt.num_init_passes; iter++)
  {
    // Set dt for this init step
//    if (m_opt.increaseDt)
//    {
      dtInit = min(dtInit*2, max_dtInit);
//    }
    thisMLPtr = this;
    for (int lev = m_level; lev <= finest_level; lev++)
    {
      thisMLPtr->dt(dtInit);
    }

    if (s_verbosity >= 3)
    {
      pout() << "Initial pressure calculation ("<< iter<< "/" <<  m_opt.num_init_passes <<") with dt = " << dtInit << endl;
    }


    thisMLPtr = this;
    for (int lev = lbase; lev <= finest_level; lev++)
    {

      thisMLPtr->backupTimestep(); // save current 'new' states
      thisMLPtr->initializeLevelPressure(cur_time, dtInit); // advance current new states by a small dt

      bool doLambdaReflux = true;


      // Also advect lambda to see how this behaves
      thisMLPtr->advectLambda(doLambdaReflux);

      // Get plot fields
      thisMLPtr->getExtraPlotFields();

      thisMLPtr = thisMLPtr->getFinerLevel();
    }

    thisMLPtr = this;

    // Now let's recompute freestream preservation
    AMRRefluxLambda();
    int numLevels = finest_level+1;
    Vector<LevelData<FArrayBox>*> amrLambda(numLevels);
    Vector<RefCountedPtr<LevelData<FluxBox> > > amrPorosityFace(numLevels);
    Vector<RefCountedPtr<LevelData<FArrayBox> > > temp(numLevels);
    Vector<LevelData<FArrayBox> * > temp2(numLevels);

    fillAMRVelPorosity(temp2, amrPorosityFace, temp);
    fillAMRLambda(amrLambda);

    if (m_opt.writePressureInitFields)
    {
      char filename[300];
      sprintf(filename, "pressureInit.%d.hdf5", iter);
      writeAMRHierarchy(filename);
    }

    thisMLPtr = this;
    for (int lev = lbase; lev <= finest_level; lev++)
    {

      bool dontReplacePressure = true;

      if (m_opt.initResetStates)
      {
        thisMLPtr->restartTimestepFromBackup(dontReplacePressure); // refill 'new' states
      }
      thisMLPtr->time(cur_time);
      thisMLPtr->dt(dtSave[lev]);

      thisMLPtr->setFluxRegistersZero();

      thisMLPtr = thisMLPtr->getFinerLevel();
    }


    if (s_verbosity >= 5)
    {
      pout() << "AMRLevelMushyLayer::initializeGlobalPressure - finished init pass " << iter  << endl;
    }

  } // end loop over init passes

  // Get this for writing out
  thisMLPtr = this;

  for (int lev = lbase; lev <= finest_level; lev++)
  {
    // turns this back on
    thisMLPtr->m_usePrevPressureForUStar=orig_add_subtract_grad_p;

    thisMLPtr->getExtraPlotFields();
    if (m_opt.doEulerPart)
    {
    thisMLPtr->m_projection.unscaledPi(*thisMLPtr->m_scalarNew[ScalarVars::m_pressure], m_dt);
    }
    thisMLPtr = thisMLPtr->getFinerLevel();
  }

}

Real AMRLevelMushyLayer::computeDtInit(int finest_level)
{
  Real dtLevel;
  Real dtInit = 10000000.0;
  AMRLevelMushyLayer* thisMLPtr = this;
  for (int lev = m_level; lev <= finest_level; lev++)
  {
    if (m_time - m_dt <= 0)
    {
      dtLevel = thisMLPtr->computeInitialDt();
    }
    else
    {
      dtLevel = thisMLPtr->computeDt();
    }

    if (dtLevel < dtInit && dtLevel != 0)
    {
      dtInit = dtLevel;
    }

    thisMLPtr = thisMLPtr->getFinerLevel();
  }

  return dtInit;

}

// -------------------------------------------------------------
// this function does all the level-based operations for
// initializing the pressure
void AMRLevelMushyLayer::initializeLevelPressure(Real a_currentTime,
                                                 Real a_dtInit)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::initializeLevelPressure " << m_level  << endl;
  }

  Real new_time = a_currentTime + a_dtInit;
  m_time = new_time;
  m_dt = a_dtInit;

  // We're not going to bother doing refluxing after initialisation,
  // do doesn't really matter if we store FR updates or not
  bool doFRupdates = false;

  // Computes advection velocities, and CC velocities if applicable
  if (solvingFullDarcyBrinkman())
  {
    IntVect ivGhost = m_numGhostAdvection*IntVect::Unit;
    LevelData<FArrayBox> advectionSourceTerm(m_grids, SpaceDim, ivGhost);
    computeAdvectionVelSourceTerm(advectionSourceTerm);

    computeAdvectionVelocities(advectionSourceTerm);

    // Just initialising so don't advect/diffuse any scalars here
    computeCCvelocity(advectionSourceTerm, m_time-m_dt, m_dt, doFRupdates, m_opt.init_compute_uDelu);

  }
  else
  {
    calculateTimeIndAdvectionVel(m_time, m_advVel);
  }

  // Get the advection velocity as a CC variable so we can write it out
  EdgeToCell(m_advVel, *m_vectorNew[VectorVars::m_advectionVel]);

}

void AMRLevelMushyLayer::computeAllVelocities(bool doFRupdates)
{

  if (solvingFullDarcyBrinkman())
  {
    IntVect ivGhost = m_numGhostAdvection*IntVect::Unit;
    LevelData<FArrayBox> advectionSourceTerm(m_grids, SpaceDim, ivGhost);
    computeAdvectionVelSourceTerm(advectionSourceTerm);

    computeAdvectionVelocities(advectionSourceTerm);

    // Just initialising so don't advect/diffuse any scalars here

    computeCCvelocity(advectionSourceTerm, m_time-m_dt, m_dt, doFRupdates);

  }
  else
  {
    calculateTimeIndAdvectionVel(m_time, m_advVel);
  }
}



/*******/
void AMRLevelMushyLayer::initialGrid(const Vector<Box>& a_newGrids)
{

  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::initialGrid (level " << m_level << ")" << endl;
  }
  // Save original grids and load balance
  m_level_grids = a_newGrids;
  Vector<int> procs;
  LoadBalance(procs, a_newGrids);
  m_grids = DisjointBoxLayout(a_newGrids, procs, m_problem_domain);

  // We need to perform a check on the time here.
  // Basically, the time on this level when initialising can't be before the time on the coarser
  // level, however when restarting if this level isn't defined yet it will have time = 0 which causes issues
  if (m_level > 0)
  {
    Real crseTime = getCoarserLevel()->m_time;
    if (crseTime > m_time)
    {
      m_time = crseTime;
    }
  }

  // Define old and new state data structures
  initDataStructures();

  // Create data structures
  createDataStructures();

  // Set up operators and stuff
  levelSetup();
}

// Refactored this so we can use when restarting
void AMRLevelMushyLayer::initDataStructures()
{
  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::initDataStructures (level " << m_level << ")" << endl;
  }
  m_scalarNew.resize(m_numScalarVars);
  m_scalarOld.resize(m_numScalarVars);
  //  m_dScalar.resize(m_numScalarVars);
  m_scalarRestart.resize(m_numScalarVars);
  //  m_scalarFluxNew.resize(m_numScalarVars);
  //  m_scalarFluxOld.resize(m_numScalarVars);

  m_vectorNew.resize(m_numVectorVars);
  m_vectorOld.resize(m_numVectorVars);
  m_dVector.resize(m_numVectorVars);
  m_vectorRestart.resize(m_numVectorVars);
}

void AMRLevelMushyLayer::createDataStructures()
{
  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::createDataStructures()" << endl;
  }

  IntVect ivGhost = m_numGhost * IntVect::Unit;
  IntVect advectionGhost = m_numGhostAdvection *IntVect::Unit;

  m_advVel.define(m_grids, 1, advectionGhost);
  m_advVelOld.define(m_grids, 1, advectionGhost);
  m_advVelNew.define(m_grids, 1, advectionGhost);

  m_frameAdvVel.define(m_grids, 1, advectionGhost);
  m_totalAdvVel.define(m_grids, 1, advectionGhost);

  m_saltFluxTop.define(m_grids, SpaceDim);
  m_saltFluxBottom.define(m_grids, SpaceDim);
  m_dPorosity_dt.define(m_grids, 1, advectionGhost);

  setValLevel(m_saltFluxTop, 0.0);
  setValLevel(m_saltFluxBottom, 0.0);
  setValLevel(m_dPorosity_dt, 0.0);

  for (int scalarVar = 0; scalarVar < m_numScalarVars; scalarVar++)
  {
    m_scalarNew[scalarVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(m_grids, 1, ivGhost));
    m_scalarOld[scalarVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(m_grids, 1, ivGhost));
    m_scalarRestart[scalarVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(m_grids, 1, ivGhost));
  }

  for (int vectorVar = 0; vectorVar < m_numVectorVars; vectorVar++)
  {
    IntVect ghost;
    if (vectorVar == m_advSrcLapU)
    {
      ghost = m_numGhostAdvection*ivGhost;
    }
    else
    {
      ghost = ivGhost;
    }

    m_vectorNew[vectorVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(m_grids, SpaceDim, ghost));
    m_vectorOld[vectorVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(m_grids, SpaceDim, ghost));
    m_dVector[vectorVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(m_grids, SpaceDim, ghost));
    m_vectorRestart[vectorVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(m_grids, SpaceDim, ghost));
  }

  // When created, all our variables will have some bogus value 10^300
  // Leave most of our variables like this, so we that if we don't initialise memory in the algorithm we find out
  // There are a few exception though, which we want to set to 0 initially, do this here:
  for (DataIterator dit = m_advVel.dataIterator(); dit.ok(); ++dit)
  {
    m_advVel[dit].setVal(0.0);
    m_advVelOld[dit].setVal(0.0);
    m_advVelNew[dit].setVal(0.0);

    (*m_scalarNew[ScalarVars::m_streamfunction])[dit].setVal(0.0);

    (*m_vectorNew[VectorVars::m_bodyForce])[dit].setVal(m_parameters.body_force);
  }

  fillFrameVelocity();



}
void AMRLevelMushyLayer::addMeltPond(int depth, Real salinity, Real enthalpy, bool rescaleExistingSolution)
{
  // This only works on level 0 for now, as we decide where to place the pond based on number
  // of cells from the top of the domain. Should change this to do it based on distance at some point.

  CH_assert(m_level==0);

  // Basic option just replaces top cells with water.
  // More advanced option will rescale existing solution to the remaining domain

  if (rescaleExistingSolution)
  {
    // Rescale existing solute to the reduced number of cells
    /**
     * We just do a very simple linear rescaling for now, so the old domain
     * 0 < z < h becomes 0 < z < h-dh, and new variables are denoted by a `, i.e. we need a mapping
     * for H(0 < z < h) -> H`(0 < z < h-dh), and S->S`
     * We assume that the boundary conditions on the
     * new domain are the same as the old domain at the new extent, i.e.
     *  H`(h-dh) = H(h), H`(h-dh) = S(h).
     * then, H`(z) = H(z*(h-dh)/h) = H(z*(1-dh/h))
     *
     */

    // First make a duplicate of the existing solution
    // Add enough ghost vectors (num ghost cells = depth, where depth is the depth of the pond being added measured in number of grid cells)
    // to ensure that, if the domain has been decomposed into multiple boxes,
    // each box has access to enough adjacent cells to compute the new solution
    IntVect ghost_vector = depth*IntVect::Unit;
    LevelData<FArrayBox> old_enthalpy(m_grids, 1, ghost_vector);
    LevelData<FArrayBox> old_bulk_conc(m_grids, 1, ghost_vector);

    m_scalarNew[ScalarVars::m_enthalpy]->copyTo(old_enthalpy);
    m_scalarNew[ScalarVars::m_bulkConcentration]->copyTo(old_bulk_conc);

    // Fill ghost cells
    old_enthalpy.exchange();
    old_bulk_conc.exchange();


    // Work out how much we're squishing the old domain by
    /**
     * Stretching maps z->z`, via z` = z*stretching
     */
    Box domBox = m_problem_domain.domainBox();
    int top_j = domBox.bigEnd(1);
    int boxHeight = domBox.bigEnd(1)-domBox.smallEnd(1);
    Real stretching = 1 - float(depth)/float(boxHeight);

    // Don't allow use of ghost cells
    int lowest_j = domBox.smallEnd(1);

    Box shrunkBox(domBox);
    shrunkBox.growHi(1, -depth);


    // Now iterate over the new grid, and fill with values interpolated from old mesh
    for (DataIterator dit = m_scalarNew[ScalarVars::m_enthalpy]->dataIterator(); dit.ok(); ++dit)
    {
      Box b = m_grids[dit];
      b &= shrunkBox;

      FArrayBox& newEnthalpy = (*m_scalarNew[ScalarVars::m_enthalpy])[dit];
      FArrayBox& oldEnthalpy = old_enthalpy[dit];
      FArrayBox& newBulkConc = (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit];
      FArrayBox& oldBulkConc = old_bulk_conc[dit];

      int vertical_index = SpaceDim-1;

      for (BoxIterator bit(b); bit.ok(); ++bit)
      {


        // Need to find the two cells to interpolate between, along with which fraction of each to take
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, loc, m_dx);

        Real z = loc[vertical_index];
        Real z_new = z*stretching;

        Real dz = z-z_new;
        Real fractional_shift = dz/m_dx;

        int num_cells_shift = floor(fractional_shift);
        Real extra_shift = fractional_shift - num_cells_shift;

        pout() << fractional_shift << ", " << extra_shift << endl;

        int z_j = iv[vertical_index] + num_cells_shift+1;
        int z_j_lower = iv[vertical_index] + num_cells_shift;

        z_j = max(lowest_j, z_j);
        z_j_lower = max(lowest_j, z_j_lower);

        z_j = min(top_j, z_j);
        z_j_lower = min(top_j, z_j_lower);


        IntVect upper = iv;
        upper[vertical_index] = z_j;

        IntVect lower = iv;
        lower[vertical_index] = z_j_lower;

        newEnthalpy(iv) = oldEnthalpy(lower)*(1-extra_shift) + extra_shift*oldEnthalpy(upper);
        newBulkConc(iv) = oldBulkConc(lower)*(1-extra_shift) + extra_shift*oldBulkConc(upper);

      }
    }

  }

  for (DataIterator dit = m_scalarNew[ScalarVars::m_enthalpy]->dataIterator(); dit.ok(); ++dit)
  {
    Box b = m_grids[dit];
    Box domBox = m_problem_domain.domainBox();
    int top_j = domBox.bigEnd(1);
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (iv[1] > top_j-depth)
      {
        (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) = enthalpy;
        (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit](iv) = salinity;
      }
    }
  }
}

void AMRLevelMushyLayer::addMeltPond()
{
  // This is an old method which I haven't used for a while

  if (m_opt.meltPondDepth > 0)
  {
    addMeltPond(m_opt.meltPondDepth, m_parameters.bcValBulkConcentrationHi[1], m_parameters.bcValEnthalpyHi[1]);
  }

}

void AMRLevelMushyLayer::postInitialGrid(const bool a_restart)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevel::postInitialGrid (level " << m_level << ")" << endl;
  }

  // Set up operators and stuff
  levelSetup();
  defineSolvers(m_time);

  if (a_restart)
  {
    // Turn this off so we can do continuation easier
//    m_doAutomaticRestart = false;

    // This will add a meltpond if m_opt.meltPondDepth is > 0
    addMeltPond();

    // if we've asked to horizontally average the solution before restarting, do it here
    // only need to do enthalpy and bulk salinity (assuming we're solving Darcy's equation).
    // this won't work with Darcy-Brinkman at the moment as we don't average velocity yet.
    // For Darcy-Brinkman we should actually probably just set U=0
    if (m_opt.horizAverageRestart)
    {
      horizontallyAverage(*m_scalarNew[ScalarVars::m_enthalpy], *m_scalarNew[ScalarVars::m_enthalpy]);
      horizontallyAverage(*m_scalarNew[ScalarVars::m_bulkConcentration], *m_scalarNew[ScalarVars::m_bulkConcentration]);
    }


    if (m_opt.restartPerturbation != 0)
    {

      // For now just apply to enthalpy - eventually let this be a choice
      addPerturbation(m_opt.restart_perturbation_var, m_opt.restartPerturbation, m_opt.perturbationWavenumber);
      copyNewToOldStates();
      updateEnthalpyVariables();
    }

    // Calculate initial sum of salt
    if (m_level == 0)
    {
      Vector<LevelData<FArrayBox>* > S_new;

      Vector<AMRLevelMushyLayer*> mlVect;
      Vector<DisjointBoxLayout> grids;
      Vector<int> refRat;
      ProblemDomain lev0Dom;
      Real lev0Dx;
      getHierarchyAndGrids(mlVect, grids, refRat, lev0Dom, lev0Dx);

      int nLevels = mlVect.size();

      S_new.resize(nLevels);

      AMRLevelMushyLayer* ml = this;
      int lev = 0;
      while(ml && ml->m_scalarNew.size() > 0)
      {
        S_new[lev] = ml->m_scalarNew[ScalarVars::m_bulkConcentration];

        ml = ml->getFinerLevel();
        lev++;
      }

      // Calculate sums over valid regions
      AMRSaltSum_old = ::computeSum(S_new, refRat, lev0Dx, Interval(0,0), 0);
    }
  }

  // Should also initialise pressure if we're restarting

  if (a_restart)
  {

    // Get the pressure from the checkpoint file

    // Note that the pressure we store in the checkpoint file is not scaled by dt,
    // whilst that stored in m_projection is.

    LevelData<FArrayBox>& pi = m_projection.Pi();
    LevelData<FArrayBox>& phi = m_projection.phi();
    m_scalarNew[ScalarVars::m_pressure]->copyTo(pi);
    m_scalarNew[ScalarVars::m_pressure]->copyTo(phi);

    // Fill pressure BCs
    this->fillScalars(pi, m_time, m_pressure, false, false);
    this->fillScalars(phi, m_time, m_pressure, false, false);

    if (solvingFullDarcyBrinkman())
    {
      for (DataIterator dit = pi.dataIterator(); dit.ok(); ++dit)
      {
        pi[dit].divide(m_dt);
      }
    }

    if(m_opt.initialize_pressures)
    {
      if (solvingFullDarcyBrinkman())
      {
        initializeGlobalPressure(m_dt,
                                 false); // not initialising, but restarting
      }
      else
      {
        AMRLevelMushyLayer* lev = this->getCoarsestLevel();
        initTimeIndependentPressure(lev);
      }
    }

    // Only do this on level 0 to ensure all other levels are setup
    if (m_level == 0)
    {

      // Need to do some CF BCs here
      AMRLevelMushyLayer* ml = getFinerLevel();
      int lev = 1;

      // This will loop over all levels, even those which aren't currently active
      while(ml)
      {
        AMRLevelMushyLayer* amrMLcrse = ml->getCoarserLevel();
        const DisjointBoxLayout& crseGrids = amrMLcrse->m_grids;
        if (crseGrids.size() > 0)
        {
          RefCountedPtr<LevelData<FArrayBox> > crsePressurePtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(crseGrids, 1));

          amrMLcrse->fillScalars(*crsePressurePtr, m_time, ScalarVars::m_pressure, true);

          // Make sure this is defined before trying to use it
          if (ml->m_quadCFInterpScalar.isDefined())
          {
            ml->m_quadCFInterpScalar.coarseFineInterp(ml->m_projection.Pi(), *crsePressurePtr);
          }
        }

        ml = ml->getFinerLevel();
        lev++;
      }
    }
  }
  // Setup diagnostics
  Vector<DiagnosticNames> diagsToPrint;
  diagsToPrint.push_back(DiagnosticNames::diag_time);
  diagsToPrint.push_back(DiagnosticNames::diag_dt);

  diagsToPrint.push_back(DiagnosticNames::diag_dSdt);
  diagsToPrint.push_back(DiagnosticNames::diag_dUdt);
  diagsToPrint.push_back(DiagnosticNames::diag_dTdt);

  switch(m_parameters.physicalProblem)
  {
    case PhysicalProblems::m_convectionMixedPorous:


      diagsToPrint.push_back(DiagnosticNames::diag_Nu);
      diagsToPrint.push_back(DiagnosticNames::diag_NuLeft);
      diagsToPrint.push_back(DiagnosticNames::diag_NuMiddle);
      diagsToPrint.push_back(DiagnosticNames::diag_NuRight);

      diagsToPrint.push_back(DiagnosticNames::diag_maxVhalf);
      diagsToPrint.push_back(DiagnosticNames::diag_maxUhalf);

      break;

    case PhysicalProblems::m_mushyLayer:

      diagsToPrint.push_back(DiagnosticNames::diag_averageVerticalSaltFlux);
      diagsToPrint.push_back(DiagnosticNames::diag_L2FsVertDiffusion);
      diagsToPrint.push_back(DiagnosticNames::diag_L2FsVertFluid);
      diagsToPrint.push_back(DiagnosticNames::diag_L2FsVertFrame);
      diagsToPrint.push_back(DiagnosticNames::diag_L1FsVertDiffusion);
      diagsToPrint.push_back(DiagnosticNames::diag_L1FsVertFluid);
      diagsToPrint.push_back(DiagnosticNames::diag_L1FsVertFrame);
      diagsToPrint.push_back(DiagnosticNames::diag_L0FsVertDiffusion);
      diagsToPrint.push_back(DiagnosticNames::diag_L0FsVertFluid);
      diagsToPrint.push_back(DiagnosticNames::diag_L0FsVertFrame);
      diagsToPrint.push_back(DiagnosticNames::diag_soluteFluxTop);
      diagsToPrint.push_back(DiagnosticNames::diag_soluteFluxBottom);
      diagsToPrint.push_back(DiagnosticNames::diag_saltFluxAbsMismatch);
      diagsToPrint.push_back(DiagnosticNames::diag_heatFluxAbsMismatch);
//      diagsToPrint.push_back(DiagnosticNames::diag_soluteFluxSponge);
      diagsToPrint.push_back(DiagnosticNames::diag_heatFluxBottom);
      diagsToPrint.push_back(DiagnosticNames::diag_heatFluxTop);
      diagsToPrint.push_back(DiagnosticNames::diag_dUdt);
      diagsToPrint.push_back(DiagnosticNames::diag_dSdt);
      diagsToPrint.push_back(DiagnosticNames::diag_dTdt);
      diagsToPrint.push_back(DiagnosticNames::diag_chimneySpacing);
      diagsToPrint.push_back(DiagnosticNames::diag_chimneyWidth);
      diagsToPrint.push_back(DiagnosticNames::diag_HorizAvSalinity0);
      diagsToPrint.push_back(DiagnosticNames::diag_HorizAvSalinity20);
      diagsToPrint.push_back(DiagnosticNames::diag_HorizAvSalinity40);
      diagsToPrint.push_back(DiagnosticNames::diag_HorizAvSalinity60);
      diagsToPrint.push_back(DiagnosticNames::diag_avSalinity);
      diagsToPrint.push_back(DiagnosticNames::diag_mushDepth);
      diagsToPrint.push_back(DiagnosticNames::diag_mushyAverageBulkConc);
      diagsToPrint.push_back(DiagnosticNames::diag_mushyAveragePorosity);
      diagsToPrint.push_back(DiagnosticNames::diag_mushyVol);
      diagsToPrint.push_back(DiagnosticNames::diag_Fs10);
      diagsToPrint.push_back(DiagnosticNames::diag_Fs20);
      diagsToPrint.push_back(DiagnosticNames::diag_Fs30); // diag_Fs40
      diagsToPrint.push_back(DiagnosticNames::diag_Fs40);
      diagsToPrint.push_back(DiagnosticNames::diag_Fs50);

      break;

    default:
      // Don't do anything else by default
      break;

  }

  m_diagnostics.setPrintDiags(diagsToPrint);

  if (m_level == 0 && procID() == 0)
  {

    m_diagnostics.printHeader();
  }

}


DisjointBoxLayout AMRLevelMushyLayer::grids()
{
  return m_grids;
}

void AMRLevelMushyLayer::shiftData(int dir, int distance)
{
  // Backup all data
  Vector<RefCountedPtr<LevelData<FArrayBox> > > previousScal, previousVect;
  LevelData<FluxBox> prevAdvVel;

  previousScal.resize(m_numScalarVars);
  previousVect.resize(m_numVectorVars);

  IntVect advectionGhost = m_numGhostAdvection*IntVect::Unit;
  IntVect ivGhost = IntVect::Unit;

  Interval scalInterval(0,0);
  Interval vectInterval(0, SpaceDim-1);

  prevAdvVel.define(m_grids, 1, advectionGhost);

  m_advVel.copyTo(scalInterval, prevAdvVel, scalInterval);

  for (int scalarVar = 0; scalarVar < m_numScalarVars; scalarVar++)
  {
    previousScal[scalarVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(m_grids, 1, ivGhost));

    m_scalarNew[scalarVar]->copyTo(scalInterval, *previousScal[scalarVar], scalInterval);

  }

  for (int vectorVar = 0; vectorVar < m_numVectorVars; vectorVar++)
  {
    previousVect[vectorVar] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(m_grids, SpaceDim, ivGhost));

    m_vectorNew[vectorVar]->copyTo(vectInterval, *previousVect[vectorVar], vectInterval);
  }

  // only works when we have one grid
  CH_assert(m_grids.size() == 1);

  Box domBox = m_problem_domain.domainBox();

  // Loop over all grid points and replace data
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    Box b = m_grids[dit];
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect ivTo = bit();
      IntVect ivFrom = ivTo - distance*BASISV(dir);

      // Need to ensure ivFrom wraps around periodic domain correctly
      if (!domBox.contains(ivFrom))
      {
        if (ivFrom[dir] >  domBox.bigEnd(dir))
        {
          int wrapDistance = ivFrom[dir] - domBox.bigEnd(dir);
          ivFrom[dir] = wrapDistance - 1;
        }
        else
        {
          int wrapDistance =  domBox.smallEnd(dir) - ivFrom[dir];
          ivFrom[dir] = domBox.bigEnd(dir) - (wrapDistance-1);
        }
      }

      for (int scalarVar = 0; scalarVar < m_numScalarVars; scalarVar++)
      {
        (*m_scalarNew[scalarVar])[dit](ivTo) = (*previousScal[scalarVar])[dit](ivFrom);
      }
      for (int var = 0; var < m_numVectorVars; var++)
      {
        (*m_vectorNew[var])[dit](ivTo) = (*previousVect[var])[dit](ivFrom);
      }
    }
  }
}

/// Changed to allow changes in either x or y directions

void AMRLevelMushyLayer::reshapeData(DisjointBoxLayout newGrids, ProblemDomain newDomain)
{
  // Need data structures to backup old data
  RefCountedPtr<LevelData<FArrayBox> > previousScal, previousVect;
  LevelData<FluxBox> prevAdvVel;

  IntVect advectionGhost = m_numGhostAdvection*IntVect::Unit;
  IntVect ivGhost = IntVect::Unit;

  Interval scalInterval(0,0);
  Interval vectInterval(0, SpaceDim-1);

  prevAdvVel.define(m_grids, 1, advectionGhost);
  previousScal = RefCountedPtr<LevelData<FArrayBox> >(
      new LevelData<FArrayBox>(m_grids, 1, ivGhost));
  previousVect = RefCountedPtr<LevelData<FArrayBox> >(
      new LevelData<FArrayBox>(m_grids, SpaceDim, ivGhost));

  m_advVel.copyTo(scalInterval, prevAdvVel, scalInterval);
  m_advVel.define(newGrids, 1, advectionGhost); // reshape
  prevAdvVel.copyTo(scalInterval,m_advVel, scalInterval); // copy back

  m_totalAdvVel.define(newGrids, 1, advectionGhost);
  m_frameAdvVel.define(newGrids, 1, advectionGhost);
  fillFrameVelocity();

  pout() << "Reshaping scalars" << endl;
  for (int scalarVar = 0; scalarVar < m_numScalarVars; scalarVar++)
  {
    m_scalarNew[scalarVar]->copyTo(scalInterval, *previousScal, scalInterval);
    m_scalarNew[scalarVar]->define(newGrids, 1, ivGhost); //reshape
    previousScal->copyTo(scalInterval, *m_scalarNew[scalarVar], scalInterval); // copy back
  }

  pout() << "Reshaping vectors" << endl;
  for (int vectorVar = 0; vectorVar < m_numVectorVars; vectorVar++)
  {
    m_vectorNew[vectorVar]->copyTo(vectInterval, *previousVect, vectInterval);
    m_vectorNew[vectorVar]->define(newGrids, SpaceDim, ivGhost); //reshape
    previousVect->copyTo(vectInterval, *m_vectorNew[vectorVar], vectInterval); // copy back
  }

  ProblemDomain oldDomain = m_problem_domain;
  m_problem_domain = newDomain;

  m_grids = newGrids;

  // Extend in both directions as necessary
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    Box newDomBox = newDomain.domainBox();

    if (newDomBox.contains(oldDomain.domainBox()))
    {
      pout() << "New domain contains old domain box" << endl;
      pout() << "New domain is " << newDomBox << endl;
      pout() << "Old domain is " << oldDomain.domainBox() << endl;
      Box filledBox = oldDomain.domainBox();

      // Keep looping until the filledBox has the same extent in this direction as the new domain box
      //      while(!filledBox.eq(newDomBox))
      while(filledBox.size(dir) != newDomBox.size(dir))
      {
        // Grow in lo and hi directions separately, as we may have
        // changed the domain assymetrically

        for (SideIterator sit = SideIterator(); sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();
          Box borderBox;

          if (side == Side::Lo)
          {
            borderBox = adjCellLo(filledBox, dir, 1);
          }
          else
          {
            borderBox = adjCellHi(filledBox, dir, 1);
          }

          // Only need to fill these cells if they fall in the new domain
          bool needToGrowDir = newDomBox.contains(borderBox);


          for (DataIterator dit = newGrids.dataIterator(); dit.ok(); ++dit)
          {
            Box b = newGrids[dit];
            Box thisBorderBox = borderBox;
            thisBorderBox &= b;

            Box copyBox;
            if (side == Side::Lo)
            {
              copyBox = adjCellHi(thisBorderBox, dir, 1);
            }
            else
            {
              copyBox = adjCellLo(thisBorderBox, dir, 1);
            }


            // Scalars
            for (int scalarVar = 0; scalarVar < m_numScalarVars; scalarVar++)
            {
              FArrayBox& phi = (*m_scalarNew[scalarVar])[dit];

              phi.copy(phi, copyBox, 0, thisBorderBox, 0, 1);

            }

            // Also need to do vectors in case U is time dependent
            for (int var = 0; var < m_numVectorVars; var++)
            {
              FArrayBox& phi = (*m_vectorNew[var])[dit];

              phi.copy(phi, copyBox, 0, thisBorderBox, 0, SpaceDim);

            }

          } // loops over grids

          //      filledBox.grow(dir, 1);
          if (needToGrowDir)
          {
            filledBox.growDir(dir, sit(), 1);

          }
        } // loop over sides

      } // while filled box != new domain box

    }
  } // end loop over directions


}


void AMRLevelMushyLayer::dx(Real newDx)
{
  m_dx = newDx;
}

Real AMRLevelMushyLayer::dx()
{
  return m_dx;
}


