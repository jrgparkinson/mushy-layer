#include "computeSum.H"
#include "AMRLevelMushyLayer.H"
#include "analyticSolns.H"
#include "SetValLevel.H"

// For custom math expressions
//#include <iostream>
//#include "muParser.h"


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
      define(   amrMLcrse->m_cfl, amrMLcrse->m_domainWidth,
                amrMLcrse->m_refineThresh, amrMLcrse->m_tagBufferSize, amrMLcrse->m_useLimiting);
    } else {
      MayDay::Error(
          "AMRLevelMushyLayer::define: a_coarserLevelPtr is not castable to AMRLevelMushyLayer*");
    }
  }

//  if (m_level == 0)
//  {
//    m_parameters.printParameters();
//  }

  // Why wasn't this line here before?
  m_problem_domain = a_problemDomain;

  m_numCells = m_problem_domain.domainBox().size();

  m_domainHeight = m_domainWidth*m_numCells[SpaceDim-1]/m_numCells[0];

  // Compute the grid spacing
  m_dx = m_domainWidth / m_numCells[0];

  m_numGhost = 1;
  m_numGhostAdvection = 4;

  m_scalarVarNames = Vector<string>(m_numScalarVars, string("scalar"));
  m_scalarVarNames[m_enthalpy] = string("Enthalpy");
  m_scalarVarNames[m_bulkConcentration] = string("Bulk concentration");
  m_scalarVarNames[m_temperature] = string("Temperature");
  m_scalarVarNames[m_porosity] = string("Porosity");
  m_scalarVarNames[m_liquidConcentration] = string("Liquid concentration");
  m_scalarVarNames[m_solidConcentration] = string("Solid concentration");
  m_scalarVarNames[m_pressure] = string("Pressure");
  m_scalarVarNames[m_permeability] = string("Permeability");
  m_scalarVarNames[m_lambda] = string("lambda");
  m_scalarVarNames[m_lambda_porosity] = string("lambda_porosity");
  m_scalarVarNames[m_enthalpySolidus] = string("Enthalpy solidus");
  m_scalarVarNames[m_enthalpyLiquidus] = string("Enthalpy liquidus");
  m_scalarVarNames[m_enthalpyEutectic] = string("Enthalpy eutectic");
  m_scalarVarNames[m_temperatureAnalytic] = string("T analytic");
  m_scalarVarNames[m_porosityAnalytic] = string("Porosity analytic");
  m_scalarVarNames[m_saltEqnSrcGodunov] = string("Salt src term godunov");
  m_scalarVarNames[m_saltEqnSrcFiniteDiff] = string("Salt src term finite");
  m_scalarVarNames[m_saltEqnOp] = string("Salt op");
  m_scalarVarNames[m_Terr] = string("T err");
  m_scalarVarNames[m_enthalpyOp] = string("H opp");
  m_scalarVarNames[m_enthalpySrc] = string("H src");
  m_scalarVarNames[m_divUadv] = string("div U face");
  m_scalarVarNames[m_dHdt] = string("dHdt");
  m_scalarVarNames[m_dSdt] = string("dSdt");
  m_scalarVarNames[m_averageVerticalFlux] = string("average vertical solute flux");
  m_scalarVarNames[m_soluteFluxAnalytic] = string("analytic vertical solute flux");
  m_scalarVarNames[m_verticalFlux] = string("vertical solute flux");
  m_scalarVarNames[m_saltResidual] = string("salt residual");
  m_scalarVarNames[m_divU] = string("div U cell");
  m_scalarVarNames[m_averageHeatFlux] = string("average vertical heat flux");
  m_scalarVarNames[m_streamfunction] = string("streamfunction");
  m_scalarVarNames[m_vorticity] = string("vorticity");
  m_scalarVarNames[m_FsVertFluid] = string("FsVertFluid");
  m_scalarVarNames[m_FsVertFrame] = string("FsVertFrame");
  m_scalarVarNames[m_FsVertDiffusion] = string("FsVertDiffusion");
  m_scalarVarNames[m_divUcorr] = string("Div U correction");
  m_scalarVarNames[m_pi] = string("pi");
  m_scalarVarNames[m_phi] = string("phi");
  m_scalarVarNames[m_MACBC] = string("MAC projection BC");
  m_scalarVarNames[m_MACrhs] = string("MAC projection rhs");
  m_scalarVarNames[m_CCrhs] = string("CC projection rhs");

  m_vectorVarNames = Vector<string>(m_numVectorVars, string("vector"));
  m_vectorVarNames[m_fluidVel] = string("Darcy velocity");
  m_vectorVarNames[m_U_porosity] = string("U divided by porosity");
  m_vectorVarNames[m_Ustar] = string("U star");
  m_vectorVarNames[m_advUstar] = string("adv U star");
  m_vectorVarNames[m_advectionVel] = string("Advection velocity");
  m_vectorVarNames[m_viscousSolveSrc] = string("Viscous solve src");
  m_vectorVarNames[m_UdelU] = string("u del U");
  m_vectorVarNames[m_advectionSrc] = string("Advection Source");
  m_vectorVarNames[m_fluidVelAnalytic] = string("Darcy Vel Analytic");
  m_vectorVarNames[m_fluidVelErr] = string("Darcy Vel Error");
  //  m_vectorVarNames[m_fluidRefluxCorr] = string("Vel reflux correction");
  m_vectorVarNames[m_dUdt] = string("dUdt");
  m_vectorVarNames[m_FsDiffusion] = string("FsDiffusion");
  m_vectorVarNames[m_FsFluid] = string("FsFluid");
  m_vectorVarNames[m_Fs] = string("Fs");
  m_vectorVarNames[m_freestreamCorrection] = string("Freestream correction");
  m_vectorVarNames[m_advUpreProjection] = string("Unprojected advection vel");
  m_vectorVarNames[m_UpreProjection] = string("Unprojected CC vel");
  m_vectorVarNames[m_advectionImplicitSrc] = string("Implicit Advection Src");
  m_vectorVarNames[m_MACcorrection] = string("MAC correction");
  m_vectorVarNames[m_CCcorrection] = string("CC correction");
  //  m_vectorVarNames[m_lapUinit] = string("Lap U init");
  //  m_vectorVarNames[m_lapUfinal] = string("Lap U final");
  m_vectorVarNames[m_advSrcLapU] = string("AdvSrc Lap U");

  m_vectorVarNames[m_bodyForce] = string("Body force");
  m_vectorVarNames[m_advVelCorr] = string("Explicit advVel corr");


  bool debug = false;
  bool minimalOutput = false;
  ParmParse pp("main");
  pp.query("debug", debug);
  if (!debug)
  {
    pp.query("minimalOutput", minimalOutput);
  }

  if (minimalOutput)
  {
    m_outputScalarVars.push_back(m_enthalpy);
    m_outputScalarVars.push_back(m_bulkConcentration);
  }
  else
  {
    // These are the default output vars
    m_outputScalarVars.push_back(m_enthalpy);
    m_outputScalarVars.push_back(m_bulkConcentration);
    m_outputScalarVars.push_back(m_temperature);
    m_outputScalarVars.push_back(m_porosity);
    m_outputScalarVars.push_back(m_liquidConcentration);

    m_outputScalarVars.push_back(m_streamfunction);
    m_outputScalarVars.push_back(m_permeability);
    m_outputScalarVars.push_back(m_lambda);
//    m_outputScalarVars.push_back(m_lambda_porosity);


    if (debug)
    {


      m_outputScalarVars.push_back(m_solidConcentration);

      m_outputScalarVars.push_back(m_saltEqnSrcGodunov);
      if (m_parameters.physicalProblem == m_parameters.m_solidificationNoFlow)
      {
        m_outputScalarVars.push_back(m_temperatureAnalytic);
        m_outputScalarVars.push_back(m_porosityAnalytic);
        m_outputScalarVars.push_back(m_Terr);
      }

      m_outputScalarVars.push_back(m_enthalpyOp);
      m_outputScalarVars.push_back(m_enthalpySrc);
      m_outputScalarVars.push_back(m_divUadv);
      m_outputScalarVars.push_back(m_dHdt);
      m_outputScalarVars.push_back(m_averageVerticalFlux);
      m_outputScalarVars.push_back(m_dSdt);
      m_outputScalarVars.push_back(m_soluteFluxAnalytic);
      m_outputScalarVars.push_back(m_verticalFlux);
      m_outputScalarVars.push_back(m_saltResidual);
      m_outputScalarVars.push_back(m_divU);
      m_outputScalarVars.push_back(m_averageHeatFlux);

      m_outputScalarVars.push_back(m_vorticity);
      m_outputScalarVars.push_back(m_FsVertFrame);
      m_outputScalarVars.push_back(m_FsVertFluid);
      m_outputScalarVars.push_back(m_FsVertDiffusion);
      //  m_outputScalarVars.push_back(m_divUcorr);
      m_outputScalarVars.push_back(m_pi);
      m_outputScalarVars.push_back(m_phi);
      m_outputScalarVars.push_back(m_MACBC);

      m_outputScalarVars.push_back(m_MACrhs);
      m_outputScalarVars.push_back(m_CCrhs);

    }

  }

  if (minimalOutput)
  {
    m_outputVectorVars.push_back(m_advectionVel);
  }
  else
  {
    // These are the default output vars

    m_outputVectorVars.push_back(m_fluidVel);
    m_outputVectorVars.push_back(m_advectionVel);

//    m_outputVectorVars.push_back(m_FsDiffusion);
//    m_outputVectorVars.push_back(m_FsFluid);
    m_outputVectorVars.push_back(m_Fs);


    if (debug)
    {
      m_outputVectorVars.push_back(m_FsDiffusion);
      m_outputVectorVars.push_back(m_FsFluid);


      //  m_outputVectorVars.push_back(m_U_porosity);
      m_outputVectorVars.push_back(m_Ustar);
      m_outputVectorVars.push_back(m_advUstar);

      //  m_outputVectorVars.push_back(m_fluidVelAnalytic);
      //  m_outputVectorVars.push_back(m_fluidVelErr);
      //  m_outputVectorVars.push_back(m_fluidRefluxCorr);
      m_outputVectorVars.push_back(m_advectionSrc);
      m_outputVectorVars.push_back(m_viscousSolveSrc);
      //  m_outputVectorVars.push_back(m_dUdt);
      m_outputVectorVars.push_back(m_freestreamCorrection);
      m_outputVectorVars.push_back(m_advUpreProjection);
      m_outputVectorVars.push_back(m_UpreProjection);
      m_outputVectorVars.push_back(m_advectionImplicitSrc);
      m_outputVectorVars.push_back(m_CCcorrection);
      m_outputVectorVars.push_back(m_MACcorrection);

      m_outputVectorVars.push_back(m_UdelU);

      m_outputVectorVars.push_back(m_advVelCorr);

      //  m_outputVectorVars.push_back(m_advSrcLapU);
      //  m_outputVectorVars.push_back(m_lapUfinal);


    }

  }

  m_numOutputComps = m_outputScalarVars.size() + m_outputVectorVars.size()*SpaceDim;

  // Now sort out what we need for checkpoint files
  // This should be a fairly small list - trying to make these files as small as possible
  // to save storage space on the disk. Should only output the bare minimum needed to
  // restart simulations
  m_chkVectorVars.push_back(m_fluidVel);

  m_chkScalarVars.push_back(m_enthalpy);
  m_chkScalarVars.push_back(m_bulkConcentration);
  m_chkScalarVars.push_back(m_pressure);
  m_chkScalarVars.push_back(m_lambda);

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::define - finished (level " << m_level << ")" << endl;
  }
}


void AMRLevelMushyLayer::levelSetup()
{
  CH_TIME("AMRLevelMushyLayer::levelSetup");

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::levelSetup (level " << m_level << ")" << endl;
  }

  m_physBCPtr->Dx(m_dx);

  AMRLevelMushyLayer* amrMLCoarserPtr = getCoarserLevel();
  AMRLevelMushyLayer* amrMLFinerPtr = getFinerLevel();

  m_hasCoarser = (amrMLCoarserPtr != NULL);
  m_hasFiner = (amrMLFinerPtr != NULL);

  int nRefCrse = 2;
  DisjointBoxLayout* crseGridsPtr = NULL;

  CCProjector *crsProj = NULL;
  CCProjector *fineProj = NULL;
  CCProjector *crsProjBackup = NULL;
  CCProjector *fineProjBackup = NULL;
  //    DisjointBoxLayout *crseGrids = NULL;

  bool scaleFineFluxes = true;

  int projectionVerbosity = 0;
  ParmParse ppProjection("projection");
  //  ParmParse ppProj("projector");
  ppProjection.query("verbosity", projectionVerbosity);


  bool usePiAdvectionBCs = true;
  ppProjection.query("usePiAdvectionBCs", usePiAdvectionBCs);

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
      m_projection.verbosity(projectionVerbosity);
      m_projection.define(m_grids, crseGridsPtr, m_problem_domain, m_dx,
                          fineProj, crsProj, nRefCrse, m_level, *m_physBCPtr, usePiAdvectionBCs);

      m_projectionBackup.verbosity(projectionVerbosity);
      m_projectionBackup.define(m_grids, crseGridsPtr, m_problem_domain, m_dx,
                                fineProjBackup, crsProjBackup, nRefCrse, m_level, *m_physBCPtr, usePiAdvectionBCs);

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
                        fineProj, crsProj, m_ref_ratio, m_level, *m_physBCPtr, usePiAdvectionBCs);
    m_projection.verbosity(projectionVerbosity);

    m_projectionBackup.define(m_grids, crseGridsPtr, m_problem_domain, m_dx,
                              fineProj, crsProj, m_ref_ratio, m_level, *m_physBCPtr, usePiAdvectionBCs);
    m_projectionBackup.verbosity(projectionVerbosity);


  }

  if (!hasFinerLevel())
  {
    m_projection.isFinestLevel(true);
    m_projectionBackup.isFinestLevel(true);

    //    if (hasCoarserLevel())
    //    {
    //      AMRLevelMushyLayer* mlCrse = getCoarserLevel();
    //      mlCrse->m_projection.isFinestLevel(false);
    //      mlCrse->m_projectionBackup.isFinestLevel(false);
    //    }
  }

  m_quadCFInterpScalar.define(m_grids, crseGridsPtr, m_dx, nRefCrse, 1,
                              m_problem_domain);
  m_quadCFInterpVector.define(m_grids, crseGridsPtr, m_dx, nRefCrse, SpaceDim,
                              m_problem_domain);

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

  ParmParse ppPatchGodunov("patchGodunov");

  // 1 -> PLM, 2 -> PPM
  int normalPredOrder = 1;

  // Use 4th order slope computations
  bool useFourthOrderSlopes = true;

  // Don't do slope limiting. Not possible to do more than one of these at the same time

  // Primitive limiting stops 2nd order convergence?
  // However definitely need it to do advection properly. Without we get porosities > 1 !
  bool usePrimLimiting = true && m_useLimiting;
  bool useCharLimiting = false && m_useLimiting;
  bool useFlattening = false; // Can't do this

  // No artificial viscosity
  bool useArtVisc = false;
  Real artVisc = -0.0;

  ppPatchGodunov.query("velOrder", normalPredOrder);
  ppPatchGodunov.query("velFourthOrderSlopes", useFourthOrderSlopes);
  ppPatchGodunov.query("velUseArtVisc", useArtVisc);
  ppPatchGodunov.query("velArtVisc", artVisc);

  m_patchGodVelocity.define(m_problem_domain, m_dx, &m_advectionPhysicsVelocity,
                            normalPredOrder, useFourthOrderSlopes, usePrimLimiting,
                            useCharLimiting, useFlattening, useArtVisc, artVisc);

  // For scalars
  int normalPredOrderHC = normalPredOrder;
  bool useFourthOrderSlopesHC = useFourthOrderSlopes;
  bool usePrimLimitingHC = usePrimLimiting;

  bool useCharLimitingHC = useCharLimiting;
  bool useFlatteningHC = useFlattening;
  bool useArtViscHC = useArtVisc;
  Real artViscHC = artVisc;

  // This speeds up convergence a little.
  // The eutectic boundary causes issues with higher order methods.
  normalPredOrderHC = 1;
  useFourthOrderSlopesHC = false;

  // this is crucial for getting convergence.
  usePrimLimitingHC = m_useLimiting;

  // This doesn't make a difference
  //  useArtViscHC = false;
  //  artViscHC = 1.0;

  ppPatchGodunov.query("HCOrder", normalPredOrderHC);
  ppPatchGodunov.query("HCFourthOrderSlopes", useFourthOrderSlopesHC);
  ppPatchGodunov.query("HCUseArtVisc", useArtVisc);
  ppPatchGodunov.query("HCArtVisc", artVisc);

  m_patchGodHC.define(m_problem_domain, m_dx, &m_advPhysHC,
                      normalPredOrderHC, useFourthOrderSlopesHC, usePrimLimitingHC,
                      useCharLimitingHC, useFlatteningHC, useArtViscHC, artViscHC);

  m_patchGodTSl.define(m_problem_domain, m_dx, &m_advPhysTSl,
                       normalPredOrderHC, useFourthOrderSlopesHC, usePrimLimitingHC,
                       useCharLimitingHC, useFlatteningHC, useArtViscHC, artViscHC);
  //                       normalPredOrder, useFourthOrderSlopes, usePrimLimiting,
  //                       useCharLimiting, useFlattening, useArtVisc, artVisc);

  // Single component solve
  m_patchGodH.define(m_problem_domain, m_dx, &m_advPhysH,
                     normalPredOrderHC, useFourthOrderSlopesHC, usePrimLimitingHC,
                     useCharLimitingHC, useFlatteningHC, useArtViscHC, artViscHC);
  m_patchGodC.define(m_problem_domain, m_dx, &m_advPhysC,
                     normalPredOrderHC, useFourthOrderSlopesHC, usePrimLimitingHC,
                     useCharLimitingHC, useFlatteningHC, useArtViscHC, artViscHC);

  m_patchGodT.define(m_problem_domain, m_dx, &m_advPhysT,
                     normalPredOrderHC, useFourthOrderSlopesHC, usePrimLimitingHC,
                     useCharLimitingHC, useFlatteningHC, useArtViscHC, artViscHC);
  m_patchGodSl.define(m_problem_domain, m_dx, &m_advPhysSl,
                      normalPredOrderHC, useFourthOrderSlopesHC, usePrimLimitingHC,
                      useCharLimitingHC, useFlatteningHC, useArtViscHC, artViscHC);


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
                                     normalPredOrder,
                                     useFourthOrderSlopes,
                                     usePrimLimiting,
                                     useCharLimiting,
                                     useFlattening,
                                     useArtVisc,
                                     artVisc);
    }
  }


}


// Refactored this so we can call it in each timestep if necessary
void AMRLevelMushyLayer::setAdvectionBCs()
{
  // todo make scalarTraceHC_IBC, scalarTraceTSl_IBC member variables,
  // and delete the old pointer when redefining them. (at the moment I think we have a memory leak)
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


void AMRLevelMushyLayer::define(const Real& a_cfl,
                                const Real& a_domainWidth, const Real& a_refineThresh,
                                const int& a_tagBufferSize,
                                const bool& a_useLimiting)
{

  ParmParse ppMain("main");
  ParmParse ppParams("parameters");
  ParmParse ppMG("amrmultigrid");

  ppMain.query("verbosity", s_verbosity);

  //  m_nusselt.resize(3);

  //  m_useSubcycling = true;
  ppMain.query("use_subcycling", m_useSubcycling);

  m_isDefined = true;
  m_cfl = a_cfl;
  m_domainWidth = a_domainWidth;
  m_refineThresh = a_refineThresh;
  m_tagBufferSize = a_tagBufferSize;
  m_useLimiting = a_useLimiting;

  // 1st/2nd order interpolation for advection
  CFinterpOrder_advection = 1;
  ppMain.query("advectionInterpOrder", CFinterpOrder_advection);

  // 1 for volume averaged, 0 for max
  m_steadyStateNormType = 1;
  ppMain.query("steadyStateNormType", m_steadyStateNormType);

  //  m_fixedDt = -1;
  ppMain.query("fixed_dt", m_fixedDt);

  //Initialise to something large
  //  m_max_dt_growth = 1.1;
  ppMain.query("max_dt_growth", m_max_dt_growth);

  m_timestepReduced = false;
  m_timestepFailed = false;
  //  m_ignoreSolverFails = true; // already defaults to true
  ppMain.query("ignoreSolverFails", m_ignoreSolverFails);
  m_solverFailRestartMethod = m_restartHalveDt;
  ppMain.query("solverFailRestartMethod", m_solverFailRestartMethod);

  m_adv_vel_centering = 0.5;
  m_adv_vel_centering_growth = 1.01;
  ppMain.query("adv_vel_centering_growth", m_adv_vel_centering_growth);

  m_dtReduction = -1;

  ppMain.query("initial_cfl", m_initial_dt_multiplier);

  m_parameters.getParameters();

  // Delete old m_physBCPtr to prevent memory issues
  // this doesn't seem to work though
  //    if (m_physBCPtr != NULL)
  //    {
  //      delete m_physBCPtr;
  //      m_physBCPtr = NULL;
  //    }
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
  Real convCrit = 1e-4;
  ppMain.query("steady_state", convCrit);
  m_diagnostics.define(diag_timescale, s_verbosity, convCrit/10);


  m_computeDiagnostics = true;
  ppMain.query("computeDiagnostics", m_computeDiagnostics);

  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::define - made diagnostics object" << endl;
  }

  m_scalePwithPorosity = true;
  ppMain.query("scalePwithChi", m_scalePwithPorosity);

  // Makes more sense to store this parameter under projection
  ParmParse ppProjection("projection");
  ppProjection.query("scalePressureWithPorosity", m_scalePwithPorosity);

  explicitDarcyTerm = false;
  ppMain.query("explicitDarcyTerm", explicitDarcyTerm);


  m_pressureScaleVar = (m_parameters.prandtl > 0) ? m_porosity : m_permeability;

  s_implicit_reflux = (m_parameters.prandtl > 0);

  m_parameters.m_nondimensionalisation = 0;
  ppMain.query("nondimensionalisation", m_parameters.m_nondimensionalisation);

  m_implicitAdvectionSolve = false;
  ppMain.query("implicitAdvectionSolve", m_implicitAdvectionSolve);

  m_solidPorosity = 0.05;


  // Default option
  m_advectionMethod =  m_porosityOutsideAdvection; // change to porosity outside advection
  ppMain.query("advectionMethod", m_advectionMethod);

  m_skipTrickySourceTerm = -1;
  ppMain.query("skipTrickySourceTime", m_skipTrickySourceTerm);

  if (m_parameters.m_nondimensionalisation == m_parameters.m_darcyTime_advectiveVel)
  {
    if (m_level == 0)
    {
    pout() << "Darcy timescale, advective velocity scale" << endl;
    }
    // To avoid dividing by 0 when Da = Pr = 0
    if (m_parameters.darcy == m_parameters.prandtl)
    {
      m_parameters.m_heatDiffusionCoeff = 1;
    }
    else
    {
      m_parameters.m_heatDiffusionCoeff = m_parameters.darcy/m_parameters.prandtl;
    }
    //    m_parameters.m_saltDiffusionCoeff = m_parameters.darcy/(m_parameters.prandtl*m_parameters.lewis);
    m_parameters.m_saltDiffusionCoeff = m_parameters.m_heatDiffusionCoeff/m_parameters.lewis;
    m_parameters.m_viscosityCoeff = m_parameters.darcy;

    m_parameters.m_buoyancyTCoeff = m_parameters.rayleighTemp*m_parameters.darcy*m_parameters.darcy*m_parameters.prandtl;
    m_parameters.m_buoyancySCoeff = m_parameters.rayleighComposition*m_parameters.darcy*m_parameters.darcy*m_parameters.prandtl;
    m_parameters.m_darcyCoeff = 1.0;
    m_parameters.m_advectionCoeff = 1.0;
  }
  else if (m_parameters.m_nondimensionalisation == m_parameters.m_diffusiveTime_advectiveVel)
  {
    if (m_level == 0)
        {
    pout() << "Diffusive timescale, advective velocity scale" << endl;
        }
    m_parameters.m_heatDiffusionCoeff = 1.0;
    m_parameters.m_saltDiffusionCoeff = 1/m_parameters.lewis;
    m_parameters.m_viscosityCoeff = m_parameters.prandtl;
    m_parameters.m_buoyancyTCoeff = m_parameters.prandtl*m_parameters.rayleighTemp;
    m_parameters.m_buoyancySCoeff = m_parameters.prandtl*m_parameters.rayleighComposition;
    m_parameters.m_darcyCoeff = m_parameters.prandtl/m_parameters.darcy;
    m_parameters.m_advectionCoeff = 1.0;
  }
  else if (m_parameters.m_nondimensionalisation == m_parameters.m_darcyTime_darcyVel)
  {
    if (m_level == 0)
        {
    pout() << "Darcy timescale, darcy velocity scale" << endl;
        }

    m_parameters.m_heatDiffusionCoeff = m_parameters.darcy/m_parameters.prandtl;
    m_parameters.m_saltDiffusionCoeff = m_parameters.m_heatDiffusionCoeff/m_parameters.lewis;
    m_parameters.m_viscosityCoeff = m_parameters.darcy;
    m_parameters.m_buoyancyTCoeff = 1.0; //m_parameters.rayleighTemp*m_parameters.darcy*m_parameters.darcy*m_parameters.prandtl;
    m_parameters.m_buoyancySCoeff = m_parameters.rayleighComposition/m_parameters.rayleighTemp; //m_parameters.rayleighComposition*m_parameters.darcy*m_parameters.darcy*m_parameters.prandtl;
    m_parameters.m_darcyCoeff = 1.0;
    m_parameters.m_advectionCoeff = m_parameters.rayleighTemp*m_parameters.darcy*m_parameters.darcy/m_parameters.prandtl;
  }
  else if (m_parameters.m_nondimensionalisation == m_parameters.m_advectiveTime_darcyVel)
  {
    if (m_level == 0)
        {
    pout() << "Advective timescale, darcy velocity scale" << endl;
        }

    m_parameters.m_heatDiffusionCoeff = 1/(m_parameters.darcy*m_parameters.rayleighTemp);///m_parameters.prandtl;
    m_parameters.m_saltDiffusionCoeff = m_parameters.m_heatDiffusionCoeff/m_parameters.lewis;
    m_parameters.m_viscosityCoeff = m_parameters.prandtl/(m_parameters.darcy*m_parameters.rayleighTemp);
    m_parameters.m_buoyancyTCoeff = m_parameters.prandtl/(m_parameters.darcy*m_parameters.rayleighTemp); //m_parameters.rayleighTemp*m_parameters.darcy*m_parameters.darcy*m_parameters.prandtl;
    m_parameters.m_buoyancySCoeff = m_parameters.m_buoyancyTCoeff*(m_parameters.rayleighComposition/m_parameters.rayleighTemp); //m_parameters.rayleighComposition*m_parameters.darcy*m_parameters.darcy*m_parameters.prandtl;
    m_parameters.m_darcyCoeff = m_parameters.prandtl/(m_parameters.darcy*m_parameters.darcy*m_parameters.rayleighTemp);
    m_parameters.m_advectionCoeff = 1.0;
  }
  else if (m_parameters.m_nondimensionalisation == m_parameters.m_buoyancyTime_advectiveVel)
  {
    if (m_level == 0)
        {
    pout() << "Buoyancy timescale, advective velocity scale" << endl;
        }

    Real R = m_parameters.rayleighComposition;
    if (R == 0)
    {
      if (m_parameters.rayleighTemp == 0)
      {
        MayDay::Error("AMRLevelMushyLayer::define - Can't nondimensionalise with buoyancy timescale if RaT=RaC=0!!");
      }
      R = m_parameters.rayleighTemp;
    }

    m_parameters.m_heatDiffusionCoeff = 1/sqrt(R*m_parameters.prandtl);
    m_parameters.m_saltDiffusionCoeff = m_parameters.m_heatDiffusionCoeff/m_parameters.lewis;

    m_parameters.m_viscosityCoeff = sqrt(m_parameters.prandtl/R);

    m_parameters.m_darcyCoeff = (1/m_parameters.darcy)*sqrt(m_parameters.prandtl/R);
    m_parameters.m_advectionCoeff = 1.0;

    // Avoid dividing by zero
    if (m_parameters.rayleighComposition != 0)
    {
      m_parameters.m_buoyancySCoeff = 1;
      m_parameters.m_buoyancyTCoeff = m_parameters.rayleighTemp/m_parameters.rayleighComposition;
    }
    else
    {
      m_parameters.m_buoyancySCoeff = 0;
      m_parameters.m_buoyancyTCoeff = 1;
    }

  }
  else
  {
    MayDay::Error("Unknown non dimensionalisation");
  }

  // Finally, option to manually set certain terms if we want (for testing purposes)
  ppMain.query("heatDiffusionCoeff", m_parameters.m_heatDiffusionCoeff );
  ppMain.query("saltDiffusionCoeff",  m_parameters.m_saltDiffusionCoeff);
  ppMain.query("viscosityCoeff", m_parameters.m_viscosityCoeff );
  ppMain.query("buoyancyTCoeff", m_parameters.m_buoyancyTCoeff );
  ppMain.query("buoyancySCoeff", m_parameters.m_buoyancySCoeff);
  ppMain.query("darcyCoeff",  m_parameters.m_darcyCoeff );
  ppMain.query("advectionCoeff",  m_parameters.m_advectionCoeff );

  // Use inviscid BCs if darcy term likely to dominate
  // This may need more care, particularly re:specific boundaries
  m_isViscous = (m_parameters.m_viscosityCoeff > 0); //(m_parameters.prandtl*m_parameters.darcy > 1e-4);

  m_viscousBCs = m_isViscous; //(m_parameters.prandtl > 0);
  ppMain.query("viscousBCs", m_viscousBCs);

  // For mushy layer calculations, we want to try and kick off the instability if we've converged
  m_doAutomaticRestart = false; //(m_parameters.physicalProblem == MushyLayerParams::m_mushyLayer);
  ppMain.query("doAutomaticRestart",m_doAutomaticRestart);

  // Regridding stuff
  // default: no smoothing
  s_do_postRegrid_smoothing = false;
  ppMain.query("do_postRegrid_smoothing", s_do_postRegrid_smoothing);
  m_regrid_smoothing_done = false; // standard initial value
  s_regrid_smoothing_coeff = 0.05;
  ppMain.query("regrid_smoothing_coeff", s_regrid_smoothing_coeff);

  m_variable_eta_factor = 1.0;
  ppMain.query("variable_eta_factor", m_variable_eta_factor);
  CH_assert(m_variable_eta_factor >= 1); // This must be >= 1, else eta will increase when it should be decreasing (and vice-versa)

  s_reflux_momentum = true;
  ppMain.query("reflux_momentum", s_reflux_momentum);

  s_reflux_normal_momentum = true;
  ppMain.query("reflux_normal_momentum", s_reflux_normal_momentum);

  //  s_reflux_scal = true;
  s_reflux_enthalpy = true;
  s_reflux_concentration = true;
  s_reflux_lambda = true;
  //  ppMain.query("reflux_scalar", s_reflux_scal);
  ppMain.query("reflux_enthalpy", s_reflux_enthalpy);
  ppMain.query("reflux_concentration", s_reflux_concentration);
  ppMain.query("reflux_lambda", s_reflux_lambda);

  s_viscous_solver_tol = 1e-10;
  ppMain.query("viscous_solver_tol", s_viscous_solver_tol);

  s_viscous_num_smooth_down = 8;
  ppMain.query("viscous_num_smooth_down", s_viscous_num_smooth_down);

  s_viscous_num_smooth_up = 8;
  ppMain.query("viscous_num_smooth_up", s_viscous_num_smooth_up);

  //	m_physicalProblem = MushyLayerParams::m_mushyLayer;
  //	ppParams.query("problem_type", m_physicalProblem);

  m_iter_plot_interval = -1;
  ppMain.query("iter_plot_interval", m_iter_plot_interval);

  //Initialization
  if (solvingFullDarcyBrinkman())
  {
    s_initialize_pressures = true;
    s_project_initial_vel = true;
    ppMain.query("initialize_pressure", s_initialize_pressures);
  }
  else
  {
    s_initialize_pressures = false;
    s_project_initial_vel = false;
  }

  s_compute_initial_VD_corr = true;

  ppMain.query("project_initial_vel", s_project_initial_vel);
  ppMain.query("initialize_VD_corr", s_compute_initial_VD_corr);

  s_regrid_init_pressure = s_initialize_pressures;
  ppMain.query("regrid_init_pressure", s_regrid_init_pressure);

  m_timeIntegrationOrder = 1;
  ppMain.query("time_integration_order", m_timeIntegrationOrder);


  // There isn't a nondimensional constant we can change to stop doing the U del (U/chi) bit, so have a switch here instead
  //  m_doEulerPart = false;
  m_doProjection = true;
  m_doSyncOperations = true;
  m_addSubtractGradP = true;
  m_enforceAnalyticSoln = false;
  m_maxDivUFace = 1e-10;

  ppMain.query("doEuler", m_doEulerPart);
  ppMain.query("doProjection", m_doProjection);
  ppMain.query("doSyncOperations", m_doSyncOperations);
  ppMain.query("addSubtractGradP", m_addSubtractGradP);
  ppMain.query("enforceAnalyticSoln", m_enforceAnalyticSoln);
  ppMain.query("maxDivUFace", m_maxDivUFace);

  Real analyticSoln=-1;
  ppMain.query("analyticSoln", analyticSoln);
  if (analyticSoln > -1)
  {
    m_enforceAnalyticSoln = true;
  }

  m_enforceGradP = (!m_doProjection); //I don't think it makes sense to do projection *and* calculate grad(P)

  // Define which scalar fields we want flux registers for
  for (int var = 0; var < m_numScalarVars; var++)
  {
    m_makeFluxRegForScalarVar[var] = false;
    m_scalarDiffusionCoeffs[var] = 0.0;
  }

  // Anything I want to be able to advect
  //  m_makeFluxRegForScalarVar[m_temperature] = true;
  m_makeFluxRegForScalarVar[m_lambda] = true;
  m_makeFluxRegForScalarVar[m_lambda_porosity] = true;
  //  m_makeFluxRegForScalarVar[m_liquidConcentration] = true;
  //  m_makeFluxRegForScalarVar[m_enthalpy] = true;
  //  m_makeFluxRegForScalarVar[m_bulkConcentration] = true;

  m_scalarDiffusionCoeffs[m_temperature] = 0.0;
  m_scalarDiffusionCoeffs[m_liquidConcentration] = 0.0;
  m_scalarDiffusionCoeffs[m_bulkConcentration] = 0.0;

  if (!(m_parameters.physicalProblem == MushyLayerParams::m_poiseuilleFlow))
  {
    //		m_scalarDiffusionCoeffs[m_temperature] = 1.0;
    m_scalarDiffusionCoeffs[m_enthalpy] = m_parameters.m_heatDiffusionCoeff;
    //		m_scalarDiffusionCoeffs[m_liquidConcentration] = 1/m_parameters.lewis;
    m_scalarDiffusionCoeffs[m_bulkConcentration] = m_parameters.m_saltDiffusionCoeff;
  }

  for (int var = 0; var < m_numVectorVars; var++)
  {
    m_makeFluxRegForVectorVar[var] = false;

  }
  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::define - about to define advectPhysics" << endl;
  }

  m_makeFluxRegForVectorVar[m_fluidVel] = true;
  m_makeFluxRegForVectorVar[m_U_porosity] = true;

  for (int var=0; var<m_numScalarVars; var++)
  {
    if (m_makeFluxRegForScalarVar[var])
    {

      m_scalarIBC[var] = getScalarIBCs(var);


    }
  }

  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::define - defined advectPhysics (level" << m_level << ")" << endl;
  }



  m_MGtype = m_FAS;
  ppMG.query("MGtype", m_MGtype);

  m_verbosity_multigrid = 0;
  ppMG.query("multigrid", m_verbosity_multigrid);

  m_lowerPorosityLimit = 1e-15;
  ppMain.query("lowPorosityLimit", m_lowerPorosityLimit);

  ppMain.query("initialPerturbation", m_initialPerturbation);

  ppMain.query("delayedPerturbaation", m_delayedPerturbation);
  ppMain.query("perturbationTime", m_perturbationTime);
  //  pp.query("initialPerturbation", m_restartPerturbation);
  ppMain.query("perturbationWavenumber", m_perturbationWavenumber);
  ppMain.query("perturbationSin", m_perturbationSin);
  ppMain.query("fixed_porosity", m_fixedPorosity);
  ppMain.query("porosity_function", m_porosityFunction);

  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::define - finished (level " << m_level << ")" << endl;
  }


}

void AMRLevelMushyLayer::initialDataDefault()
{

  Real H =  m_parameters.bcValEnthalpyLo[SpaceDim-1];
  Real C =  m_parameters.bcValBulkConcentrationLo[SpaceDim-1];

  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    (*m_scalarNew[m_enthalpy])[dit].setVal(H);
    (*m_scalarNew[m_bulkConcentration])[dit].setVal(C);

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


      (*m_scalarNew[m_enthalpy])[dit](iv) =  m_parameters.bcValEnthalpyLo[dir] + deltaH * loc[dir] + perturbation;
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

      int numRands = 2*m_domainWidth/m_dx;
      Vector<Real> rands(numRands);
      for (int i = 1; i < numRands; i++)
      {
        rands[i] = ((double) rand()/ (RAND_MAX));
      }


      // if we have periodic bcs, the perturbation must be continuous at the boundary
      if (m_problem_domain.isPeriodic(0))
      {
        // White noise initialization

        int N = 0.5*m_domainWidth/m_dx;
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
          perturbation += additionalWeight*weight*cos(2*n*M_PI*(x+phaseShift)/(m_domainWidth))*sin((y/m_domainHeight)*M_PI);
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
        perturbation = alpha*cos(2*wavenumber*M_PI*x/(m_domainWidth))*sin((y/m_domainHeight)*M_PI);

        //perturbation = alpha*cos(8*(x/domainWidth)*M_PI)*sin((y/domainHeight)*M_PI);
      }



      (*m_scalarNew[m_enthalpy])[dit](iv) = m_parameters.bcValEnthalpyLo[dir] + deltaH * loc[dir]/m_domainHeight + perturbation;
      (*m_scalarNew[m_bulkConcentration])[dit](iv) = m_parameters.bcValBulkConcentrationLo[dir] + deltaC * loc[dir]/m_domainHeight;



    }
  }

}


void AMRLevelMushyLayer::initialDataConvectionMixedPorous()
{
  // First setup diagnostics
  Vector<int> diagsToPrint;

  diagsToPrint.push_back(m_diagnostics.m_time);
  diagsToPrint.push_back(m_diagnostics.m_dt);
  diagsToPrint.push_back(m_diagnostics.m_Nu);
  m_diagnostics.setPrintDiags(diagsToPrint);
  // Now do data
  Real initVel = 0;
  ParmParse pp("init");
  pp.query("initVel", initVel);
  bool doScalarAdvectionDiffusion = true;
  ParmParse ppmain("main");
  ppmain.query("doScalarAdvectionDiffusion", doScalarAdvectionDiffusion);

  bool perturbationSin=false;
  ppmain.query("perturbationSin", perturbationSin);

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

      Real arg = 2*m_perturbationWavenumber*x*M_PI/m_domainWidth;

      if (m_problem_domain.isPeriodic(0))
      {
        //          perturbation = alpha*cos(2*x*M_PI/m_domainWidth)*sin(y*M_PI);
        perturbation = m_initialPerturbation*sin(y*M_PI/m_domainHeight);

        if (perturbationSin)
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
        perturbation = -m_initialPerturbation*sin(y*M_PI);

        if (perturbationSin)
        {
          perturbation = perturbation*sin(arg);
        }
        else
        {
          perturbation = perturbation*cos(arg);
        }

      }
      // Make sure we don't take H above H max
      perturbation -= m_initialPerturbation;


      (*m_scalarNew[m_enthalpy])[dit](iv) =Hbottom  +  (Htop - Hbottom) * loc[dir]/m_domainWidth + perturbation;
      (*m_scalarNew[m_bulkConcentration])[dit](iv) = m_parameters.bcValBulkConcentrationLo[dir];

      if (!doScalarAdvectionDiffusion)
      {
        //        (*m_scalarNew[m_enthalpy])[dit](iv) = (*m_scalarNew[m_enthalpy])[dit](iv) - 0.1*sin(M_PI*x/m_domainWidth)*cos(M_PI*y/m_domainHeight);
        //        (*m_scalarNew[m_enthalpy])[dit](iv) = Htop;
      }

      //      for (int dir=0;dir<SpaceDim;dir++)
      //      {
      //        (*m_vectorNew[m_fluidVel])[dit](iv, dir) = 0.0;
      //      }

      (*m_vectorNew[m_fluidVel])[dit](iv, 0) = initVel*cos(M_PI*y/m_domainHeight)*sin(M_PI*x/m_domainWidth);
      (*m_vectorNew[m_fluidVel])[dit](iv, 1) = -initVel*cos(M_PI*x/m_domainWidth)*sin(M_PI*y/m_domainHeight);

    }
  }
}

void AMRLevelMushyLayer::initialDataRayleighBenard()
{


  int numRands = 2*m_domainWidth/m_dx;
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

      Real perturbation = m_initialPerturbation*cos(x*M_PI)*sin((y/m_domainHeight)*M_PI);

      // if we have periodic bcs, the perturbation must be continuous at the boundary
      if (m_problem_domain.isPeriodic(0))
      {
        // White noise initialization

        int N = 0.5*m_domainWidth/m_dx;
        Real weight = m_initialPerturbation/N;

        perturbation = 0;
        for (int n = 1; n <= N; n++)
        {

          Real additionalWeight = 1*rands[n];
          Real phaseShift = rands[n];
          //                if (n == 2)
          //                {
          //                  additionalWeight = 5;
          //                }
          perturbation += additionalWeight*weight*cos(2*n*M_PI*(x+phaseShift)/(m_domainWidth))*sin((y/m_domainHeight)*M_PI);
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

        perturbation = m_initialPerturbation*cos(2*m_perturbationWavenumber*M_PI*x/(m_domainWidth))*sin((y/m_domainHeight)*M_PI);

        //perturbation = alpha*cos(8*(x/domainWidth)*M_PI)*sin((y/domainHeight)*M_PI);
      }

      (*m_scalarNew[m_enthalpy])[dit](iv) = m_parameters.bcValEnthalpyLo[1] + (m_parameters.bcValEnthalpyHi[1] - m_parameters.bcValEnthalpyLo[1]) * y/m_domainHeight + perturbation;
      //            (*m_scalarNew[m_bulkConcentration])[dit](iv) = ;


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


      //      (*m_scalarNew[m_enthalpy])[dit](iv)  = (m_parameters.Hinitial - 0.2*m_parameters.stefan*exp((y-m_domainHeight)/0.2))*(1-0.05*sin(5*M_PI * x / m_domainWidth));
      (*m_scalarNew[m_enthalpy])[dit](iv)  = m_parameters.Hinitial;

      //        valNew = m_parameters.ThetaInitial + sin(M_PI * y / m_domainHeight)*(1-0.05*cos(5*M_PI * x / m_domainWidth));
      (*m_scalarNew[m_bulkConcentration])[dit](iv) = -1 + y;

      // Shift y values because flux will be face centred

      Real yFace = y - 0.5*m_dx;
      Real magU = 10;
      Real freq = M_PI/m_domainWidth;
      Real w = -magU*((1/freq)*(1-cos(freq)-0.5))*sin(freq*yFace); //horizontally averaged w
      //w = 1.0;
      w = -magU*(1/freq)*(cos(freq) - 1.0 + 0.5*freq)*exp(-freq*yFace);
      //    w = 0;
      (*m_scalarNew[m_soluteFluxAnalytic])[dit](iv) = m_parameters.m_saltDiffusionCoeff + (1-yFace)*(m_parameters.nonDimVel + w);


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


      (*m_scalarNew[m_bulkConcentration])[dit](iv)  = -1.0;


      (*m_scalarNew[m_enthalpy])[dit](iv) = m_parameters.bcValEnthalpyLo[1] + (m_parameters.bcValEnthalpyHi[1] - m_parameters.bcValEnthalpyLo[1])*y;

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
    FArrayBox& vel = (*m_vectorNew[m_fluidVel])[dit];
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

  setValLevel(*m_vectorNew[m_fluidVel], 0.0);

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
    (*m_scalarNew[m_enthalpy])[dit].setVal(m_parameters.Hinitial);

    (*m_scalarNew[m_bulkConcentration])[dit].setVal(m_parameters.bcValBulkConcentrationLo[1]);
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

      (*m_scalarNew[m_enthalpy])[dit](iv) = m_parameters.stefan*(1+1.5*sin(M_PI*(x/m_domainWidth*2 + 0.5))*sin(M_PI*y/m_domainHeight)); //(m_parameters.Hinitial - (m_parameters.stefan+10)*exp((y-m_domainHeight)/0.2));;
      (*m_scalarNew[m_bulkConcentration])[dit](iv) = Cav + Cav*sin(M_PI*x/m_domainWidth*2)*sin(M_PI*y/m_domainHeight);


    }
  }
}

void AMRLevelMushyLayer::initialDataPorousHole()
{
  ParmParse pp("main");

  Real radius = 0.1*m_domainWidth;
    pp.query("radius", radius);

  DataIterator dit = m_grids.dataIterator();

    Real HTop = m_parameters.bcValEnthalpyHi[1];
    Real HBottom = m_parameters.bcValEnthalpyLo[1];
    Real ThetaTop = m_parameters.bcValBulkConcentrationHi[1];
    Real ThetaBottom = m_parameters.bcValBulkConcentrationLo[1];

    Real thetaBottom = m_parameters.bcValTemperatureLo[1];
    Real thetaTop = m_parameters.bcValTemperatureHi[1];
    Real ThetaLTop = m_parameters.bcValLiquidConcentrationHi[1];


     Real cx = m_domainWidth/2;
     Real cy = m_domainHeight/2;

     Real initVelScale = 1e-3;
     pp.query("initVelScale", initVelScale);

     for (dit.reset(); dit.ok(); ++dit)
     {
       BoxIterator bit((*m_scalarNew[0])[dit].box());
       for (bit.reset(); bit.ok(); ++bit)
       {
         IntVect iv = bit();
         RealVect loc;
         getLocation(iv, loc, m_dx);

         (*m_scalarNew[m_enthalpy])[dit](iv) = HTop + (HBottom-HTop)*exp(- ( pow(loc[0]-cx,2) + pow(loc[1]-cy, 2))/(pow(radius,2)));
         (*m_scalarNew[m_bulkConcentration])[dit](iv) = ThetaTop + (ThetaBottom-ThetaTop)*exp(- ( pow(loc[0]-cx,2) + pow(loc[1]-cy, 2))/(pow(radius,2)));

         (*m_vectorNew[m_fluidVel])[dit](iv, 1) = initVelScale*sin(M_PI*loc[0]/m_domainWidth);
         (*m_vectorNew[m_fluidVel])[dit](iv, 0) = 0;

         (*m_vectorNew[m_bodyForce])[dit](iv, 1) =  m_parameters.m_buoyancySCoeff*ThetaLTop \
             - m_parameters.m_buoyancyTCoeff*thetaTop ;

       }
     }


}

void AMRLevelMushyLayer::initialDataMushyLayer()
{

  //  Real alpha = 0.0;
  ParmParse pp("main");

  bool perturbationSin=false;
  pp.query("perturbationSin", perturbationSin);

  int summerProfile = -1;
  pp.query("initSummerProfile", summerProfile);
  Real mushHeight = m_domainHeight/2;
  pp.query("summerProfileMushHeight", mushHeight);

//  bool porousHole = false;
//  pp.query("porousHole", porousHole);

  bool linearGradient = false;
  pp.query("initLinear", linearGradient);

  bool initVel = false;
  pp.query("initVel", initVel);



  Real perturbation;
  //  Real Nwaves = round(m_domainWidth);
  //                     Nwaves = Nwaves*10;
  //  Real wavelength = m_domainWidth/Nwaves;
  DataIterator dit = m_grids.dataIterator();

  Real HTop = m_parameters.bcValEnthalpyHi[1];
  Real HBottom = m_parameters.bcValEnthalpyLo[1];
//  Real ThetaTop = m_parameters.bcValBulkConcentrationHi[1];
  Real ThetaBottom = m_parameters.bcValBulkConcentrationLo[1];

  Real thetaBottom = m_parameters.bcValTemperatureLo[1];
//  Real thetaTop = m_parameters.bcValTemperatureHi[1];
//  Real ThetaLTop = m_parameters.bcValLiquidConcentrationHi[1];

 if (summerProfile > 0)
  {
    Real blScale = 0.00001;

    if (summerProfile == 1)
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
        //        Real x = loc[0];
        Real y = loc[1];


        Real chiMin = 0.2;
        Real chiAv = (1+chiMin)/2;
        Real chi = chiAv - (1-chiAv)*tanh((y-mushHeight)/(blScale));

        //        Real topHeight = mushHeight/5;
        Real thetaMush = 1.01;
        Real thetaAv = (thetaMush + thetaBottom)/2;
        Real theta = thetaAv - (thetaBottom-thetaAv)*tanh((y-mushHeight)/(blScale));

        //      if (m_domainHeight-y < topHeight)
        //      {
        //        Real fractionalPos = (y - (m_domainHeight-topHeight))/topHeight;
        //        chi = chiMin*(1-fractionalPos);
        //        theta = theta*(1-fractionalPos);
        //
        //        Real tanhArg = m_domainHeight*(y - (m_domainHeight-topHeight));
        //        chi = chiMin/2 - (chiMin/2)*tanh(tanhArg);
        //        theta = theta/2 - (theta/2)*tanh(tanhArg);
        //      }


        Real H = m_parameters.stefan*chi + theta;

        //    valNew = m_parameters.ThetaBottom + (m_parameters.ThetaTop - m_parameters.ThetaBottom) * (1-exp(y/m_domainHeight))/(1-exp(1));
        Real S = 0;

        //        Real sdiff = 0.06;
        Real sice = -1.12; Real Sav = (sice + ThetaBottom)/2;
        //Real valMush = -1-2*sdiff;
        S = Sav - (Sav-sice)*tanh((y-(mushHeight*1.0))/(blScale));

        (*m_scalarNew[m_enthalpy])[dit](iv) = H;
        (*m_scalarNew[m_bulkConcentration])[dit](iv) = S;

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

        if ((m_domainHeight - y) < mushHeight/5)
        {
          H = Hmin + (y - (m_domainHeight -  mushHeight/5))*((HTop-Hmin)/(mushHeight/5));
        }
        else
        {
          H = Hav - (HBottom  - Hav)*tanh((y-mushHeight)/(0.02));
        }

        Real arg = 2*m_perturbationWavenumber*x*M_PI/m_domainWidth;

        //        H = m_parameters.Hinitial + alpha*(HTop) * y +  perturbation*y*y;
        //                     valNew = m_parameters.HBottom + (m_parameters.HTop - m_parameters.HBottom) * y * y * y +  perturbation;
        // valNew = m_parameters.HBottom + (m_parameters.HTop - m_parameters.HBottom) * y *y * y* y;
        //    valNew = m_parameters.HBottom +  perturbation*y*y;
        // H = m_parameters.Hinitial;

        if (linearGradient)
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
            perturbation = m_initialPerturbation*sin(y*M_PI/m_domainHeight);

            if (perturbationSin)
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
            perturbation = -m_initialPerturbation*sin(y*M_PI);

            if (perturbationSin)
            {
              perturbation = perturbation*sin(arg);
            }
            else
            {
              perturbation = perturbation*cos(arg);
            }

          }
          // Make sure we don't take H above H max
          perturbation -= m_initialPerturbation;
        }
        else
        {
          perturbation = 0;
        }

        H = H + perturbation;

        //    valNew = m_parameters.HBottom + (m_parameters.HTop - m_parameters.HBottom) * (1-exp(y/m_domainHeight))/(1-exp(1));

        (*m_scalarNew[m_enthalpy])[dit](iv) = H;

        if (initVel)
        {
          (*m_vectorNew[m_fluidVel])[dit](iv, 0) = sin(2*M_PI*y);
          (*m_vectorNew[m_fluidVel])[dit](iv, 1) = sin(2*M_PI*x);
        }

      }
    }

  }


}

void AMRLevelMushyLayer::initialDataIceBlock()
{

  MayDay::Error("AMRLevelMushyLayer:: initialDataIceBlock Not implemented");
  //
  //  Vector<Real> blockWidth(2);
  //  Vector<Real> blockHeight(2);
  //
  //  // For a unit domain!
  //  blockWidth[0] = 0.3; blockWidth[1] = 0.7;
  //  blockHeight[0] = 0.6; blockHeight[1] = 1.0;
  //
  //  if (scalarVar == m_enthalpy)
  //  {
  //
  //    if (x > blockWidth[0] && x < blockWidth[1] &&
  //        y > blockHeight[0] && y < blockHeight[1])
  //    {
  //      valNew = m_parameters.HTop;
  //    }
  //    else
  //    {
  //      valNew = m_parameters.HBottom;
  //    }
  //
  //
  //  } // end if enthalpy
  //  else if (scalarVar == m_bulkConcentration)
  //  {
  //    if (x > blockWidth[0] && x < blockWidth[1] &&
  //        y > blockHeight[0] && y < blockHeight[1])
  //    {
  //      // Think is fresh water?
  //      valNew = -m_parameters.compositionRatio;
  //    }
  //    else
  //    {
  //      valNew = -1.0;
  //    }
  //
  //  }

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

      if (m_porosityFunction == 0)
      {
        (*m_scalarNew[m_porosity])[dit](iv) = 1.0;
      }
      else if (m_porosityFunction == 1)
      {
        (*m_scalarNew[m_porosity])[dit](iv) = 0.1 + 0.9*x;
      }
      else if (m_porosityFunction == 2)
      {
        (*m_scalarNew[m_porosity])[dit](iv) =  0.1+0.9*exp(-(x-0.5)*(x-0.5)*100);
      }
      else if (m_porosityFunction == 3)
      {
        (*m_scalarNew[m_porosity])[dit](iv) =(exp(-(x-0.5)*(x-0.5)))*(0.5 + 0.1*sin(2*M_PI*y));
      }
      else if (m_porosityFunction == 4)
      {
        (*m_scalarNew[m_porosity])[dit](iv) = 0.1;
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

      (*m_scalarNew[m_porosity])[dit](iv) = (m_parameters.Hinitial - (m_parameters.stefan+10)*exp((y-m_domainHeight)/0.2));;
      (*m_scalarNew[m_bulkConcentration])[dit](iv) = m_parameters.ThetaInitial + (y-1)*exp(-(y-1)*(y-1)/0.01);

      //      Real fixed_porosity = -1.0;

      if (m_fixedPorosity < 0)
      {

        if (m_porosityFunction == 0)
        {
          (*m_scalarNew[m_porosity])[dit](iv)  = 1.0;
        }
        else if (m_porosityFunction == 1)
        {
          (*m_scalarNew[m_porosity])[dit](iv)  = 0.1 + 0.9*y;
        }
        else if (m_porosityFunction == 2)
        {
          (*m_scalarNew[m_porosity])[dit](iv)  = 0.01 + pow(y-0.01, 3);
        }
        else if (m_porosityFunction == 3)
        {
          //                                                                                                                      valNew = 0.5 + 0.45* sin(M_PI*x)*sin(2*M_PI*y);
          //                                                              valNew = (1.2-x)*(1-0.2*exp(1.2*y*y));
          (*m_scalarNew[m_porosity])[dit](iv)  = 0.5*( (1-tanh(20*y)) + (1-x));
        }
        else
        {
          MayDay::Error("No porosity function specified");
        }
      }
      else
      {
        (*m_scalarNew[m_porosity])[dit](iv)  = m_fixedPorosity;
      }


      Real f = - 1;
      Real fprime = 0;


      (*m_vectorNew[m_fluidVel])[dit](iv, 0)  = -x*fprime;
      (*m_vectorNew[m_fluidVel])[dit](iv, 1)  = f;


    }
  }

}

/*******/
void AMRLevelMushyLayer::initialData()
{
  pout() << "AMRLevelMushyLayer::initialData - setting initial data on " << m_level << endl;

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

    (*m_scalarNew[m_lambda])[dit].setVal(1.0);
    (*m_scalarOld[m_lambda])[dit].setVal(1.0);
    //    (*m_dScalar[m_lambda])[dit].setVal(1.0);

    (*m_scalarNew[m_enthalpy])[dit].setVal(m_parameters.Hinitial);
    (*m_scalarOld[m_enthalpy])[dit].setVal(m_parameters.Hinitial);
    //    (*m_dScalar[m_enthalpy])[dit].setVal(m_parameters.Hinitial);

    (*m_scalarNew[m_bulkConcentration])[dit].setVal(m_parameters.ThetaInitial);
    (*m_scalarOld[m_bulkConcentration])[dit].setVal(m_parameters.ThetaInitial);
    //    (*m_dScalar[m_bulkConcentration])[dit].setVal(m_parameters.ThetaInitial);

    (*m_vectorNew[m_fluidVel])[dit].setVal(0.0);
    (*m_vectorOld[m_fluidVel])[dit].setVal(0.0);
    (*m_dVector[m_fluidVel])[dit].setVal(0.0);

    (*m_vectorNew[m_bodyForce])[dit].setVal(0.0);
    (*m_vectorOld[m_bodyForce])[dit].setVal(0.0);
    (*m_dVector[m_bodyForce])[dit].setVal(0.0);



    m_advVel[dit].setVal(0.0);


  }

  ParmParse pp("main");

  // If we've specified some custom initial data, deal with it here
  if (pp.contains("initData"))
  {
    int initData = -1;
    pp.query("initData", initData);

    switch(initData)
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
      case MushyLayerParams::m_sidewallHeating:
        initialDataSidewallHeating();
        break;
      case MushyLayerParams::m_HRL:
        initialDataHRL();
        break;
      case MushyLayerParams::m_convectionMixedPorous:
        initialDataConvectionMixedPorous();
        break;
      case MushyLayerParams::m_rayleighBenard:
        initialDataRayleighBenard();
        break;
      case MushyLayerParams::m_cornerFlow:
        initialDataCornerFlow();
        break;
      case  MushyLayerParams::m_mushyLayer:
        initialDataMushyLayer();
        break;
      case MushyLayerParams::m_poiseuilleFlow:
        initialDataPoiseuille();
        break;
      case MushyLayerParams::m_soluteFluxTest:
        initialDataSoluteFlux();
        break;
      case MushyLayerParams::m_refluxTest:
        initialDataRefluxTest();
        break;
      case MushyLayerParams::m_zeroPorosityTest:
        initialDataZeroPorosityTest();
        break;
      case MushyLayerParams::m_meltingIceBlock:
        initialDataIceBlock();
        break;
      case MushyLayerParams::m_diffusion:
        initialDataDiffusion();
        break;
      case MushyLayerParams::m_vortexPair:
        initialDataVortexPair();
        break;


      default:
        initialDataDefault();
        //          case MushyLayerParams::m_:
        //            initialData(x, y, valNew, scalarVar, -1);
        //            break;
    }

  }
  addMeltPond();

  if(m_parameters.physicalProblem == MushyLayerParams::m_poiseuilleFlow)
  {
    stokesDarcyForcing(*m_scalarNew[m_temperature], 0);
  }
  else if (m_parameters.physicalProblem == MushyLayerParams::m_cornerFlow)
  {
    // Don't do anything
  }
  else
  {
    updateEnthalpyVariables();
    copyNewToOldStates();
    updateEnthalpyVariables();
  }

  calculatePermeability(); //make sure this is up to date

  // need this to get the initial timestep right
  //Initialise advVel
  CellToEdge(*m_vectorNew[m_fluidVel], m_advVel);
  CellToEdge(*m_vectorNew[m_fluidVel], m_advVelOld);
  CellToEdge(*m_vectorNew[m_fluidVel], m_advVelNew);
  //  calculateTimeIndependentVelocity(0);
  //  calculateTimeIndAdvectionVel(0, m_advVel);

  m_advVel.exchange();
  m_frameAdvVel.exchange();
  //  m_frameVel.exchange();

  calculateAnalyticSolns(m_enforceAnalyticSoln);


  ParmParse ppMain("main");
  bool initAnalyticVel=false;
  ppMain.query("initAnalyticVel", initAnalyticVel);

  if (initAnalyticVel)
  {
    fillAnalyticVel(m_advVel);
    //    fillAnalyticVel(*m_vectorNew[m_fluidVel]);

    // I think this is the most consistent way to do things

    EdgeToCell(m_advVel, *m_vectorNew[m_fluidVel]);
  }



}

void AMRLevelMushyLayer::addPerturbation(int a_var, Real alpha, int waveNumber, Real phaseShift)
{
  pout() << "Adding perturbation " << alpha << " with wavenumber " << waveNumber << endl;

  Real domainWidth = m_domainWidth;

  int numWavenumbers = 50;
  ParmParse pp("main");
  pp.query("maxRestartWavenumber", numWavenumbers);

  Vector<Real> rands1(m_numCells[0]);
  Vector<Real> rands2(m_numCells[1]);
  Vector<Real> randsk(numWavenumbers);
  for (int i = 0; i < m_numCells[0]; i++)
  {
    rands1[i] = ((double) rand()/ (RAND_MAX));
  }
  for (int i = 0; i < m_numCells[1]; i++)
  {
    rands2[i] = ((double) rand()/ (RAND_MAX));
  }
  for (int i = 0; i < numWavenumbers; i++)
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
        Real wavelength = m_domainWidth/waveNumber;
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
        for (int k=1; k<numWavenumbers; k++)
        {
          Real random = ((double) rand()/ (RAND_MAX));

          Real thisPert = (alpha/numWavenumbers)*cos(M_PI*(k*x/domainWidth + randsk[k])); // x component

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

  int analyticSolution = m_parameters.physicalProblem;
  ParmParse pp("main");
  pp.query("analyticSoln", analyticSolution);

  if (analyticSolution == MushyLayerParams::m_solidificationNoFlow)
  {

    LevelData<FArrayBox> CLanalytic(m_grids, 1, ivGhost);
    LevelData<FArrayBox> CSanalytic(m_grids, 1, ivGhost);

    ::analyticSolnSolidificationNoFlow(enthalpyAnalytic, bulkCAnalytic,
                                       m_domainHeight, m_dx,
                                       m_parameters);

    ::updateEnthalpyVariables(enthalpyAnalytic, bulkCAnalytic,
                              *m_scalarNew[m_temperatureAnalytic], CLanalytic,
                              CSanalytic, *m_scalarNew[m_porosityAnalytic],
                              m_parameters);

    ::updateEnthalpyVariables(enthalpyAnalytic, bulkCAnalytic,
                              *m_scalarOld[m_temperatureAnalytic], CLanalytic,
                              CSanalytic, *m_scalarOld[m_porosityAnalytic],
                              m_parameters);

    if (enforceSolutions)
    {


      for (DataIterator dit = m_scalarOld[m_porosityAnalytic]->dataIterator(); dit.ok(); ++dit)
      {
        (*m_scalarOld[m_enthalpy])[dit].setVal(0);
        (*m_scalarOld[m_enthalpy])[dit] += enthalpyAnalytic[dit];
        (*m_scalarNew[m_enthalpy])[dit].copy((*m_scalarOld[m_enthalpy])[dit]);

        (*m_scalarOld[m_bulkConcentration])[dit].setVal(0);
        (*m_scalarOld[m_bulkConcentration])[dit] += bulkCAnalytic[dit];
        (*m_scalarNew[m_bulkConcentration])[dit].copy((*m_scalarOld[m_bulkConcentration])[dit]);

        if (m_initialPerturbation != 0.0)
        {
          Box b = (*m_scalarNew[m_enthalpy])[dit].box();
          if (m_perturbationWavenumber > -1)
          {
            for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc;
              ::getLocation(iv, loc, m_dx);

              Real pert;
              if (m_perturbationSin)
              {
                pert = m_initialPerturbation*sin(M_PI*2*(loc[1]/m_domainHeight)*m_perturbationWavenumber)*sin(M_PI*(loc[0]/m_domainWidth)*m_perturbationWavenumber*2);
              }
              else
              {
                pert = m_initialPerturbation*sin(M_PI*2*(loc[1]/m_domainHeight)*m_perturbationWavenumber)*cos(M_PI*(loc[0]/m_domainWidth)*m_perturbationWavenumber*2);
              }

              pert *= (*m_scalarNew[m_porosity])[dit](iv);
              (*m_scalarNew[m_enthalpy])[dit](iv) += pert;
            }
          }
          else
          {

            // Add a small random perturbation everywhere

            for (BoxIterator bit = BoxIterator(b); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              Real r = ((double) rand()/ (RAND_MAX));
              (*m_scalarNew[m_enthalpy])[dit](iv) += r*m_initialPerturbation;
            }

          }
        }
      }
      updateEnthalpyVariablesOld();
      updateEnthalpyVariables();
    }
  }
  else if (analyticSolution == MushyLayerParams::m_poiseuilleFlow)
  {
    Real fixedPorosity = -1.0;
    ParmParse pp("parameters");

    for (DataIterator dit = m_vectorNew[m_fluidVel]->dataIterator(); dit.ok(); ++dit)
    {
      Box b = (*m_vectorNew[m_fluidVel])[dit].box();
      for (BoxIterator bit(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, loc, m_dx);
        Real x = loc[0];
        //				Real y = loc[1];


        //        pp.query("fixed_porosity", fixedPorosity);
        (*m_vectorNew[m_fluidVelAnalytic])[dit](iv, 0) = 0;

        if(m_fixedPorosity >= 0)
        {
          Real perm  = m_parameters.calculatePermeability(m_fixedPorosity);
          Real L = sqrt(fixedPorosity/(m_parameters.darcy*perm));
          (*m_vectorNew[m_fluidVelAnalytic])[dit](iv, 1) = m_parameters.rayleighTemp*perm*(1-(cosh(L*(0.5-x)))/(cosh(0.5*L)));
        }
        else
        {
          (*m_vectorNew[m_fluidVelAnalytic])[dit](iv, 1) = m_parameters.rayleighTemp*(1-4*(x-0.5)*(x-0.5))/8;
        }

        (*m_vectorOld[m_fluidVelAnalytic])[dit](iv, 0) = (*m_vectorNew[m_fluidVelAnalytic])[dit](iv, 0);
        (*m_vectorOld[m_fluidVelAnalytic])[dit](iv, 1) = (*m_vectorNew[m_fluidVelAnalytic])[dit](iv, 1);


        if (enforceSolutions)
        {
          (*m_vectorOld[m_fluidVel])[dit](iv, 0) = (*m_vectorNew[m_fluidVelAnalytic])[dit](iv, 0);
          (*m_vectorOld[m_fluidVel])[dit](iv, 1) = (*m_vectorNew[m_fluidVelAnalytic])[dit](iv, 1);

          (*m_vectorNew[m_fluidVel])[dit](iv, 0) = (*m_vectorNew[m_fluidVelAnalytic])[dit](iv, 0);
          (*m_vectorNew[m_fluidVel])[dit](iv, 1) = (*m_vectorNew[m_fluidVelAnalytic])[dit](iv, 1);
        }


        setPorosity((*m_scalarNew[m_porosity])[dit]);

      }

    }
  }
  else if (analyticSolution == MushyLayerParams::m_diffusion)
  {
    for (DataIterator dit = m_vectorNew[m_fluidVel]->dataIterator(); dit.ok(); ++dit)
    {
      Box b = (*m_vectorNew[m_fluidVel])[dit].box();
      for (BoxIterator bit(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, loc, m_dx);
        //        Real x = loc[0];
        Real y = loc[1];

        (*m_scalarNew[m_temperatureAnalytic])[dit](iv) = m_parameters.bcValTemperatureLo[0] - 0.5*y*(y-1);
      }
    }
  }
  else if (analyticSolution == MushyLayerParams::m_mushyLayer)
  {
    // when we ask for an analytic solution here, just generate a rough chimney shape

    ::channelFieldMushyLayer(enthalpyAnalytic, bulkCAnalytic,
                             m_domainHeight, m_domainWidth, m_dx,
                             m_parameters);

    if (enforceSolutions)
    {

      for (DataIterator dit = m_scalarOld[m_porosityAnalytic]->dataIterator(); dit.ok(); ++dit)
      {
        (*m_scalarOld[m_enthalpy])[dit].setVal(0);
        (*m_scalarOld[m_enthalpy])[dit] += enthalpyAnalytic[dit];
        (*m_scalarNew[m_enthalpy])[dit].copy((*m_scalarOld[m_enthalpy])[dit]);

        (*m_scalarOld[m_bulkConcentration])[dit].setVal(0);
        (*m_scalarOld[m_bulkConcentration])[dit] += bulkCAnalytic[dit];
        (*m_scalarNew[m_bulkConcentration])[dit].copy((*m_scalarOld[m_bulkConcentration])[dit]);

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

    if (m_fixedPorosity >= 0)
    {
      porosity = m_fixedPorosity;
    }
    else
    {
      if (m_porosityFunction == 0)
      {
        porosity = 1.0;
      }
      else if (m_porosityFunction == 1)
      {
        porosity = 0.1 + 0.9*x;
      }
      else if (m_porosityFunction == 2)
      {
        porosity =  0.1+0.9*exp(-(x-0.5)*(x-0.5)*100);
      }
      else if (m_porosityFunction == 3)
      {
        porosity =(exp(-(x-0.5)*(x-0.5)))*(0.5 + 0.1*sin(2*M_PI*y));
      }
      else if (m_porosityFunction == 4)
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
    amrVel[lev] = &(*thisLevelData->m_vectorNew[m_fluidVel]);
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
  if (m_level == 0 && m_doProjection)
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

    CCProjector& level0Proj = m_projection;

    // set physical boundary conditions on velocity

    thisLevelData = this;

    VelBCHolder velBC(m_physBCPtr->uStarFuncBC(m_viscousBCs));

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
    if (s_project_initial_vel)
    {
      level0Proj.initialVelocityProject(amrVel, amrPorosityFace, amrPorosity, homoBC);
    }

    // Compute initial volume discrepancy correction
    // For this we need to advect lambda on each level
    // to get a measure of initial freestream violation
    if (s_compute_initial_VD_corr)
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
      Real initDt;
      if (this->m_fixedDt > 0)
      {
        initDt = m_fixedDt;
      }
      else
      {
        initDt= thisML->computeInitialDt();
        initDt = 1e-4*initDt;
        thisML->dt(initDt);
        thisML->computeInitAdvectionVel();

        // Now get a better estimate of the initial dt
        initDt = thisML->computeInitialDt();
        initDt = 1e-4*initDt; //make it a lot smaller to be safe
      }

      pout() << "Computing initial VD correction with dt = " << initDt << endl;

      for (int lev = 0; lev < numLevels; lev++)
      {
        thisML->dt(initDt);

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
        amrLambda[lev] = thisML->m_scalarNew[m_lambda];
        //        amrPorosity[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(thisML->m_grids, 1));
        //        thisML->fillScalarFace(*amrPorosity[lev],m_time, m_porosity, true, true);
        //thisML->m_scalarNew[m_porosity];
        thisML = thisML->getFinerLevel();
      }

      level0Proj.doInitialSyncOperations(amrVel, amrLambda,
                                         amrPorosityFace, amrPorosity,
                                         m_time, m_dt); // false - don't add to grad_elambda, simply reset it

      // Reset original fields incl. lambda
      thisML = this;
      for (int lev = 0; lev < numLevels; lev++)
      {
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

    if (s_initialize_pressures)
    {
      if (s_verbosity >= 3)
            {
              pout() << "AMRlevelMushyLayer::postInitialize - initialize pressures" << endl;
            }

      Real dtInit = computeDtInit(numLevels-1);
      Real scale = 0.1;
      ParmParse ppMain("main");
      ppMain.query("init_dt_scale", scale);
      dtInit *= scale;

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
      S_new[lev] = ml->m_scalarNew[m_bulkConcentration];
      H_new[lev] = ml->m_scalarNew[m_enthalpy];

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

void AMRLevelMushyLayer::resetLambda()
{
  setValLevel(*m_scalarOld[m_lambda], 1.0);
  setValLevel(*m_scalarNew[m_lambda], 1.0);
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
    thisML->doExplicitReflux(m_lambda);
    thisML = thisML->getCoarserLevel();
  }
}

void AMRLevelMushyLayer::computeInitAdvectionVel()
{
  if (solvingFullDarcyBrinkman())
  {

    IntVect ivGhost = m_numGhostAdvection*IntVect::Unit;
    LevelData<FArrayBox> advectionSourceTerm(m_grids, SpaceDim, ivGhost);
    computeAdvectionVelSourceTerm(advectionSourceTerm);

    computeAdvectionVelocities(advectionSourceTerm);

    // Just initialising so don't advect/diffuse any scalars here
  }
  else
  {
    calculateTimeIndAdvectionVel(m_time, m_advVel);
  }
}

//bool AMRLevelMushyLayer::newLevelAdded()
//{
//  AMRLevelMushyLayer* thisMLPtr = this->getCoarsestLevel();
//
//  bool newLevel=false;
//
//  while (thisMLPtr)
//  {
//    newLevel = newLevel || thisMLPtr->m_newLevel;
//    thisMLPtr = thisMLPtr->getFinerLevel();
//  }
//
//  return newLevel;
//}
// -------------------------------------------------------------
// this function manages the pressure initialization after
// initialization and regridding
void AMRLevelMushyLayer::initializeGlobalPressure(Real dtInit, bool init)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::initializeGlobalPressure (level " << m_level << ")" << endl;
  }
  ParmParse ppMain("main");

  // Reset adv vel centering if we're adding new levels?
  // Actually, always reset adv vel if we're initialising pressures
  //  if (newLevelAdded())
  //  {
  Real advVelCentering = 0.1;
  ppMain.query("init_advVel_centering", advVelCentering);
  setAdvVelCentering(advVelCentering);
  //  }

  // this was 5 originally
  int s_num_init_passes = 1;

  ppMain.query("num_init_passes", s_num_init_passes);

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
  if (init)
  {
    ppMain.query("restart_newTime", cur_time);
  }

  bool increaseDt = true;
  ppMain.query("init_increase_dt", increaseDt);

  bool init_add_subtract_grad_p = false;
  ppMain.query("init_add_subtract_grad_p", init_add_subtract_grad_p);

  //  Real dtLevel;
  //  Real dtInit = 10000000.0;
  Real max_dtInit = 0.5*computeDtInit(finest_level);
  if (dtInit < 0)
  {
    dtInit = computeDtInit(finest_level);
    dtInit *= 0.5;
  }

  // Save dt for each level so it can be reset later
  Vector<Real> dtSave(finest_level + 1, 1.0e8);
  thisMLPtr = this;
  for (int lev = m_level; lev <= finest_level; lev++)
  {
    dtSave[lev] = thisMLPtr->m_dt;

    thisMLPtr = thisMLPtr->getFinerLevel();
  }

  // Set add/subtract grad(p) = false ?
  thisMLPtr = this;
  for (int lev = m_level; lev <= finest_level; lev++)
  {
    // This is now done later
    //    thisMLPtr->dt(dtInit);
    //		thisMLPtr->time(dtInit);

    //todo change this properly
    if (!init_add_subtract_grad_p)
    {
      thisMLPtr->m_addSubtractGradP=false;
    }

    thisMLPtr = thisMLPtr->getFinerLevel();

  }

  // now loop through levels and compute estimate of pressure
  // gradually increase dt during initialization steps

  thisMLPtr = this;
  int lbase = thisMLPtr->m_level;
  for (int iter = 0; iter < s_num_init_passes; iter++)
  {
    // Set dt for this init step
    if (increaseDt)
    {
      dtInit = min(dtInit*2, max_dtInit);
    }
    thisMLPtr = this;
    for (int lev = m_level; lev <= finest_level; lev++)
    {
      thisMLPtr->dt(dtInit);
    }

    pout() << "Initial pressure calculation ("<< iter<< "/" <<  s_num_init_passes <<") with dt = " << dtInit << endl;


    thisMLPtr = this;
    for (int lev = lbase; lev <= finest_level; lev++)
    {

      thisMLPtr->backupTimestep(); // save current 'new' states
      thisMLPtr->initializeLevelPressure(cur_time, dtInit); // advance current new states by a small dt

      bool doLambdaReflux = true;
      //      if (iter == s_num_init_passes-1)
      //      {
      //        doLambdaReflux = true;
      //      }

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

    //    thisMLPtr->m_projection.doPostRegridOps(amrLambda, amrPorosityFace, dtInit, m_time);

    bool writePressureInitFields = false;
    ppMain.query("writePressureInitFields", writePressureInitFields);
    if (writePressureInitFields)
    {
      char filename[300];
      sprintf(filename, "pressureInit.%d.hdf5", iter);
      writeAMRHierarchy(filename);
    }

    thisMLPtr = this;
    for (int lev = lbase; lev <= finest_level; lev++)
    {

      bool dontReplacePressure = true;

      bool resetStates = true;
      ppMain.query("init_resetStates", resetStates);
      if (resetStates)
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
    // todo - do this properly
    thisMLPtr->m_addSubtractGradP=true;

    thisMLPtr->getExtraPlotFields();
    thisMLPtr->m_projection.unscaledPi(*thisMLPtr->m_scalarNew[m_pressure], m_dt);
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

  ParmParse ppMain("main");
  Real new_time = a_currentTime + a_dtInit;
  m_time = new_time;
  m_dt = a_dtInit;
  //  swapOldAndNewStates(); // don't do this - advance the old fields

  // We're not going to bother doing refluxing after initialisation,
  // do doesn't really matter if we store FR updates or not
  bool doFRupdates = false;

  // Computes advection velocities, and CC velocities if applicable
  //  computeAllVelocities(doFRupdates);
  if (solvingFullDarcyBrinkman())
  {
    IntVect ivGhost = m_numGhostAdvection*IntVect::Unit;
    LevelData<FArrayBox> advectionSourceTerm(m_grids, SpaceDim, ivGhost);
    computeAdvectionVelSourceTerm(advectionSourceTerm);

    computeAdvectionVelocities(advectionSourceTerm);
    // advVelCentering);

    // todo - this probably shouldn't be turned off
    bool compute_uDelu = true;

    ppMain.query("init_compute_uDelu", compute_uDelu);

    // Just initialising so don't advect/diffuse any scalars here
    computeCCvelocity(advectionSourceTerm, m_time-m_dt, m_dt, doFRupdates, compute_uDelu);

  }
  else
  {
    calculateTimeIndAdvectionVel(m_time, m_advVel);
  }

  // Get the advection velocity as a CC variable so we can write it out
  EdgeToCell(m_advVel, *m_vectorNew[m_advectionVel]);

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
    //    m_dScalar[scalarVar] = RefCountedPtr<LevelData<FArrayBox> >(
    //        new LevelData<FArrayBox>(m_grids, 1, ivGhost));
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

  for (DataIterator dit = m_advVel.dataIterator(); dit.ok(); ++dit)
  {
    m_advVel[dit].setVal(0.0);
    m_advVelOld[dit].setVal(0.0);
    m_advVelNew[dit].setVal(0.0);


    (*m_vectorNew[m_bodyForce])[dit].setVal(0.0);
  }

  fillFrameVelocity();



}


void AMRLevelMushyLayer::addMeltPond()
{
  ParmParse pp("main");

  int meltPondDepth = 0;
  pp.query("meltPondDepth", meltPondDepth);

  if (meltPondDepth > 0)
  {

    for (DataIterator dit = m_scalarNew[m_enthalpy]->dataIterator(); dit.ok(); ++dit)
    {
      Box b = m_grids[dit];
      Box domBox = m_problem_domain.domainBox();
      int top_j = domBox.bigEnd(1);
      for (BoxIterator bit(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        if (iv[1] > top_j-meltPondDepth)
        {
          (*m_scalarNew[m_enthalpy])[dit](iv) = m_parameters.bcValEnthalpyHi[1];
          (*m_scalarNew[m_bulkConcentration])[dit](iv) = m_parameters.bcValBulkConcentrationHi[1];
        }
      }
    }

  }

}

void AMRLevelMushyLayer::postInitialGrid(const bool a_restart)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevel::postInitialGrid (level " << m_level << ")" << endl;
  }

  // Need this here so we call it after restarts

  // Set up operators and stuff
  levelSetup();
  defineSolvers(m_time);

  if (m_level == 0 && procID() == 0)
  {
    m_diagnostics.printHeader();
  }



  if (a_restart)
  {
    // Turn this off so we can do continuation easier
    m_doAutomaticRestart = false;

    Real perturbation = 0.0;
    if (perturbation != 0)
    {
      // For now just apply to enthalpy - eventually let this be a choice
      addPerturbation(m_enthalpy, m_restartPerturbation, m_perturbationWavenumber);
      updateEnthalpyVariables();
    }


    addMeltPond();

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
        S_new[lev] = ml->m_scalarNew[m_bulkConcentration];

        ml = ml->getFinerLevel();
        lev++;
      }

      // Calculate sums over valid regions
      AMRSaltSum_old = ::computeSum(S_new, refRat, lev0Dx, Interval(0,0), 0);
    }
  }

  // Should also initialise pressure if we're restarting

  if (a_restart
      && s_initialize_pressures)
  {

    // Get the pressure from the checkpoint file

    // Note that the pressure we store in the checkpoint file is not scaled by dt,
    // whilst that stored in m_projection is.

    LevelData<FArrayBox>& pi = m_projection.Pi();
    m_scalarNew[m_pressure]->copyTo(pi);
    for (DataIterator dit = pi.dataIterator(); dit.ok(); ++dit)
    {
      pi[dit].divide(m_dt);
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

          amrMLcrse->fillScalars(*crsePressurePtr, m_time, m_pressure, true);


          // Make sure this is defined before trying to use it
          if (ml->m_quadCFInterpScalar.isDefined())
          {
            ml->m_quadCFInterpScalar.coarseFineInterp(ml->m_projection.Pi(), *crsePressurePtr);
            //ml->m_quadCFInterpScalar.coarseFineInterp(*ml->m_scalarOld[m_pressure], *crsePressurePtr);
          }
        }

        ml = ml->getFinerLevel();
        lev++;
      }


      // todo  - was this required?
      // Potentially do some composite initialisation here?
      /*
      int finest_level = lev-1; // because we've just looped over all levels
      Real dtInit = computeDtInit(finest_level);
      Real scale = 0.1;
      ParmParse ppMain("main");
      ppMain.query("init_dt_scale", scale);
      dtInit *= scale;

      initializeGlobalPressure(dtInit);
       */
    }

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

