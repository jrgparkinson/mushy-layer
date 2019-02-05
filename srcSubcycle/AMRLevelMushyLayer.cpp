
#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iomanip>

#include "parstream.H"

#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "BoxIterator.H"
#include "AMRIO.H"
#include "computeSum.H"
#include "computeNorm.H"
#include "AMRMultiGrid.H"
#include "AMRFASMultiGrid.H"

#include "AMRLevelMushyLayer.H"

#include "ParmParse.H"

#include "SetValLevel.H"
#include "Gradient.H"
#include "MushyLayerUtils.H"

#include "timeInterp.H"
#include "VCAMRPoissonOp2.H"
#include "DarcyBrinkmanOp.H"
#include "Divergence.H"
#include "ExtrapFillPatch.H"
#include "PiecewiseLinearFillPatchFace.H"
#include "PhysIBC.H"
#include "AMRNonLinearMultiCompOp.H"
#include "EnthalpyVariablesF_F.H"
#include "utils_F.H"


#include "BackwardEuler.H"

#include "Channel.h"


#include "NamespaceHeader.H"

BiCGStabSolver<LevelData<FArrayBox> > AMRLevelMushyLayer::s_botSolverUStar;
RelaxSolver<LevelData<FArrayBox> > AMRLevelMushyLayer::s_botSolverHC;

/*******/
AMRLevelMushyLayer::~AMRLevelMushyLayer()
{
}

/*******/
void AMRLevelMushyLayer::getHierarchyAndGrids(
    Vector<AMRLevelMushyLayer*>& a_hierarchy,
    Vector<DisjointBoxLayout>& a_grids, Vector<int>& a_refRat,
    ProblemDomain& a_lev0Dom, Real& a_lev0Dx)
{
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  int nlevels = hierarchy.size();

  a_hierarchy.resize(nlevels);
  a_refRat.resize(nlevels);
  a_grids.resize(nlevels);

  AMRLevelMushyLayer* coarsestLevel = (AMRLevelMushyLayer*) (hierarchy[0]);
  a_lev0Dx = coarsestLevel->m_dx;
  a_lev0Dom = coarsestLevel->m_problem_domain;

  for (int ilev = 0; ilev < nlevels; ilev++)
  {
    AMRLevelMushyLayer* adLevel = (AMRLevelMushyLayer*) (hierarchy[ilev]);

    a_hierarchy[ilev] = adLevel;
    a_grids[ilev] = adLevel->m_grids;
    a_refRat[ilev] = adLevel->m_ref_ratio;
  }
}

/*******/
void AMRLevelMushyLayer::defineSolvers(Real a_time)
{

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::defineSolvers" << endl;
  }
  CH_TIME("AMRLevelMushyLayer::defineSolvers");
  //  Real old_time = m_time-m_dt;
  //  Real new_time = m_time;

  bool a_homogeneous = false;
  Real alpha = 1.0;
  Real beta = 1.0;

  IntVect ivGhost = m_numGhost * IntVect::Unit;

  int numSmoothUp=4, numSmoothDown=1, numMG=1, maxIter=10, mgverb=0, bottomSolveIterations=40;
  Real tolerance=1e-10, hang=1e-10, normThresh=1e-10;
  int relaxMode = 1; // 1=GSRB, 4=jacobi
  bool useRelaxBottomSolverForHC = true;

  ParmParse ppAmrmultigrid("HCMultigrid");
  ppAmrmultigrid.query("num_smooth_up", numSmoothUp);
  ppAmrmultigrid.query("num_mg", numMG);
  ppAmrmultigrid.query("hang_eps", hang);
  ppAmrmultigrid.query("norm_thresh", normThresh);
  ppAmrmultigrid.query("tolerance", tolerance);
  ppAmrmultigrid.query("max_iter", maxIter);
  ppAmrmultigrid.query("verbosity", mgverb);
  ppAmrmultigrid.query("numSmoothDown", numSmoothDown);
  ppAmrmultigrid.query("relaxMode", relaxMode);
  ppAmrmultigrid.query("bottomSolveIterations", bottomSolveIterations);
  ppAmrmultigrid.query("useRelaxBottomSolver", useRelaxBottomSolverForHC);

  s_botSolverHC.m_imax = bottomSolveIterations;

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

  s_botSolverUStar.m_verbosity = max(mgverb - 2, 0);
  s_botSolverHC.m_verbosity = max(mgverb - 2, 0);

  {

    CH_TIME("AMRLevelMushyLayer::defineSolvers::defineViscousOp");
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    BCHolder viscousBCcomp = m_physBCPtr->velFuncBC(dir, m_viscousBCs);

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
    if (s_uStarAMRMG[idir] == NULL)
    {
      s_uStarAMRMG[idir] =
          RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >(
              new AMRMultiGrid<LevelData<FArrayBox> >());
      s_uStarAMRMG[idir]->setSolverParameters(numSmoothUp, numSmoothUp, numSmoothUp,
                                              numMG, maxIter, tolerance, hang, normThresh);
      s_uStarAMRMG[idir]->m_verbosity = mgverb;
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

    enthalpy[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], 1, ivGhost));
    bulkConcentration[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[lev], 1, ivGhost));

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

      amrML->fillScalarFace(*porosityFace[lev], a_time, m_porosity, true, secondOrder);

      //      amrML->fillScalars(*porosity[lev], m_time-m_dt/2, m_porosity, true, secondOrder);

      amrML->fillScalars(*enthalpy[lev], a_time, m_enthalpy, true , secondOrder);
      amrML->fillScalars(*bulkConcentration[lev], a_time, m_bulkConcentration, true, secondOrder);

      amrML->fillScalars(*enthalpySolidus[lev],a_time, m_enthalpySolidus, true, secondOrder);
      amrML->fillScalars(*enthalpyLiquidus[lev], a_time, m_enthalpyLiquidus, true, secondOrder);
      amrML->fillScalars(*enthalpyEutectic[lev],a_time, m_enthalpyEutectic, true ,secondOrder);
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
        (*bCoef[lev])[dit][dir].mult(-m_scalarDiffusionCoeffs[m_enthalpy], Hcomp);
        (*bCoef[lev])[dit][dir].mult(-m_scalarDiffusionCoeffs[m_bulkConcentration], Ccomp);
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
               relaxMode, porosityEdgeBC);

  HCOpFact = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(HCop);

  int maxAMRlevels = hierarchy.size();

  if (m_MGtype == m_standard)
  {
    MayDay::Error("Standard multigrid no longer supported");
  }
  else if (m_MGtype == m_FAS)
  {
    //FAS multigrid



      s_multiCompFASMG = RefCountedPtr<AMRFASMultiGrid<LevelData<FArrayBox> > >(
          new AMRFASMultiGrid<LevelData<FArrayBox> >());

      if (useRelaxBottomSolverForHC)
      {
        s_multiCompFASMG->define(lev0Dom, *HCOpFact, &s_botSolverHC,
                                 maxAMRlevels);
      }
      else
      {
        // Borrow the BiCGStab bottom solver from U star
        s_multiCompFASMG->define(lev0Dom, *HCOpFact, &s_botSolverUStar,
                                      maxAMRlevels);
      }

      s_multiCompFASMG->setSolverParameters(numSmoothDown, numSmoothUp, numSmoothUp, numMG,
                                            maxIter, tolerance, hang, normThresh);
      s_multiCompFASMG->m_verbosity = mgverb;



    // Not doing TGA yet
    s_enthalpySalinityTGA = RefCountedPtr<LevelTGA>(
        new LevelTGA(grids, refRat, lev0Dom, HCOpFact,
                     s_multiCompFASMG));


    s_enthalpySalinityBE = RefCountedPtr<LevelBackwardEuler>(
        new LevelBackwardEuler(grids, refRat, lev0Dom, HCOpFact,
                               s_multiCompFASMG));

  }
  else
  {
    MayDay::Error("Unknown multigrid type specified");
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

  int numSmooth=2, numMG=1, maxIter=10, mgverb=0;
  Real tolerance=1e-10, hang=1e-10, normThresh=1e-10;

  ParmParse ppAmrmultigrid("VelocityMultigrid");
  ppAmrmultigrid.query("num_smooth", numSmooth);
  ppAmrmultigrid.query("num_mg", numMG);
  ppAmrmultigrid.query("hang_eps", hang);
  ppAmrmultigrid.query("norm_thresh", normThresh);
  ppAmrmultigrid.query("tolerance", tolerance);
  ppAmrmultigrid.query("max_iter", maxIter);
  ppAmrmultigrid.query("verbosity", mgverb);


  //  Real half_time = m_time - m_dt / 2;
  //  Real new_time = m_time;
  Real old_time = m_time-m_dt;

  //    int nlevels = (m_level == 0) ? 1 : 2; // only 1 level if m_level = 0, as no coarser grids in this case
  //    int coarsestLevel = (m_level == 0) ? 0 : m_level - 1; // the coarsest level we care about
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

  for (int lev = coarsestLevel; lev < nlevels; lev++)
  {

    int relativeLev = lev - coarsestLevel;

    AMRLevelMushyLayer* thisLevel =
        (AMRLevelMushyLayer*) (hierarchy[lev]);
    DisjointBoxLayout& levelGrids = thisLevel->m_grids;

    solverGrids[relativeLev] = levelGrids;

    bCoef[relativeLev] = RefCountedPtr<LevelData<FluxBox> >(
        new LevelData<FluxBox>(levelGrids, 1, ivGhost)); // = 1
    aCoef[relativeLev] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(levelGrids, 1, ivGhost)); // = 1
    cCoef[relativeLev] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(levelGrids, 1, ivGhost)); // = porosity.pi_0/(pi)


    // Only actually fill levels if they're at m_level or coarser
    if (lev <= m_level)
    {
      LevelData<FArrayBox> porosity(levelGrids, 1, ivGhost);
      LevelData<FArrayBox> permeability(levelGrids, 1, ivGhost);

      // Fill these at the new time
      // this is correct for backward euler (where we evaluate the darcy term
      // at the new time) but we should modify this for TGA where we evaluate U at
      // different times within an update
      if (m_timeIntegrationOrder == 2 && m_time-m_dt == 0)
      {
        // Only display this warning at initial timestep
        //        MayDay::Warning("AMRlevelMushyLayer::defineUstarSolver - darcy term coefficient not time centred correctly for TGA");
      }

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
        if (explicitDarcyTerm)
        {
          (*cCoef[relativeLev])[dit].setVal(0.0);
        }
        else
        {
          (*cCoef[relativeLev])[dit].copy(porosity[dit]);
          (*cCoef[relativeLev])[dit].divide(permeability[dit]);
          (*cCoef[relativeLev])[dit].mult(m_parameters.m_darcyCoeff);
        }

        //                                                      (*cCoef[relativeLev])[dit].setVal(0.0); // testing

        // This is what multiplies du/dt
        (*aCoef[relativeLev])[dit].setVal(1.0);

        // this is what multiplies laplacian(u). Note the minus sign.
        (*bCoef[relativeLev])[dit].setVal(- m_parameters.m_viscosityCoeff );
      } // end loop over boxes

    } // end if level <= m_level
  } // end loop over levels

  for (int idir = 0; idir < SpaceDim; idir++)
  {

    BCHolder viscousBC = m_physBCPtr->velFuncBC(idir, m_viscousBCs);

    RefCountedPtr<DarcyBrinkmanOpFactory> vcamrpop = RefCountedPtr<DarcyBrinkmanOpFactory>(new DarcyBrinkmanOpFactory());
    vcamrpop->define(lev0Dom, allGrids, refRat, lev0Dx, viscousBC,
                     0.0, aCoef, -1.0, bCoef, cCoef); // Note that we should set m_dt*etc in bCoef, not beta!

    s_uStarOpFact[idir] = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(vcamrpop); // m_UstarVCAMRPOp[idir]);

    s_uStarAMRMG[idir]->define(lev0Dom, *s_uStarOpFact[idir],
                               &s_botSolverUStar, nlevels);

    s_uStarAMRMG[idir]->setSolverParameters(numSmooth, numSmooth, numSmooth,
                                            numMG, maxIter, tolerance, hang, normThresh);
  }
}

void AMRLevelMushyLayer::defineUstarSolver(	Vector<RefCountedPtr<LevelBackwardEuler> >& UstarBE,
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
    UstarBE[idir] = RefCountedPtr<LevelBackwardEuler> (new LevelBackwardEuler(allGrids, refRat, lev0Dom, s_uStarOpFact[idir], s_uStarAMRMG[idir]));
    UstarTGA[idir] = RefCountedPtr<LevelTGA> (new LevelTGA(allGrids, refRat, lev0Dom, s_uStarOpFact[idir], s_uStarAMRMG[idir]));
  }

}

/********/

/********/
/********/



void AMRLevelMushyLayer::getCoarseScalarDataPointers(const int a_scalarVar,
                                                     LevelData<FArrayBox>** a_coarserDataOldPtr,
                                                     LevelData<FArrayBox>** a_coarserDataNewPtr,
                                                     LevelFluxRegister** a_coarserFRPtr, LevelFluxRegister** a_finerFRPtr,
                                                     Real& a_tCoarserOld, Real& a_tCoarserNew)
{
  *a_coarserDataOldPtr = NULL;
  *a_coarserDataNewPtr = NULL;
  *a_coarserFRPtr = NULL;
  *a_finerFRPtr = NULL;

  a_tCoarserOld = 0.0;
  a_tCoarserNew = 0.0;

  // A coarser level exists
  if (m_hasCoarser)
  {
    AMRLevelMushyLayer* coarserPtr = getCoarserLevel();

    // Recall that my flux register goes between my level and the next
    // finer level
    *a_coarserFRPtr = &(*coarserPtr->m_fluxRegisters[a_scalarVar]);

    *a_coarserDataOldPtr = &(*coarserPtr->m_scalarOld[a_scalarVar]);
    *a_coarserDataNewPtr = &(*coarserPtr->m_scalarNew[a_scalarVar]);

    a_tCoarserNew = coarserPtr->m_time;
    a_tCoarserOld = a_tCoarserNew - coarserPtr->m_dt;
  }

  // A finer level exists
  if (m_hasFiner)
  {
    // Recall that my flux register goes between my level and the next
    // finer level

    // Check flux register is defined - if we are restarting,
    // then it might not be defined yet on finer (currently unused) levels
    if (m_fluxRegisters[a_scalarVar])
    {
      *a_finerFRPtr = &(*m_fluxRegisters[a_scalarVar]);
    }
  }
}

bool AMRLevelMushyLayer::convergedToSteadyState()
{
  // Check to see if we've crashed
  if(crashed())
  {
    pout() << "==============================" << endl;
    pout() << "Simulation crashed, stopping gracefully." << endl;
    pout() << "==============================" << endl;

    return true;
  }

  if (m_timestepFailed)
  {
    return false;
  }



  // Don't consider levels other than base level
  // Reflux corrections etc. can prevent us ever reaching a 'steady' state on these levels,
  // even when d/dt on the coarsest level is ~ 10^{-10}
  if (m_level > 0)
  {
    return true;
  }

  Real steadyStateCondition = 1e-3;
  ParmParse ppMain("main");
  ppMain.query("steady_state", steadyStateCondition);

  Real Tnorm = convergedToSteadyState(m_enthalpy);
  Real Cnorm = convergedToSteadyState(m_bulkConcentration);
  Real Unorm =  convergedToSteadyState(m_fluidVel, true);

  if (m_computeDiagnostics)
  {
    m_diagnostics.addDiagnostic(Diagnostics::m_dTdt, m_time, Tnorm);
    m_diagnostics.addDiagnostic(Diagnostics::m_dSdt, m_time, Cnorm);
    m_diagnostics.addDiagnostic(Diagnostics::m_dUdt, m_time, Unorm);

    // Can print diagnostics now if on processor 0
    bool printDiagnostics = (m_level == 0 && procID() ==0 ); // only print results on proc 0
    if (printDiagnostics)
    {
      m_diagnostics.addDiagnostic(Diagnostics::m_dt, m_time, m_dt);
      m_diagnostics.printDiagnostics(m_time);
    }
  }

  //  bool Tstalled = m_diagnostics.movingAverageHasConverged(Diagnostics::m_dTdt, m_time, m_dt)                                                                                                                                                                                                                                                                                                                                                                                       && Tnorm < 10*steadyStateCondition;
  //  bool Cstalled = m_diagnostics.movingAverageHasConverged(Diagnostics::m_dSdt, m_time, m_dt)                                                                                                                                                                                                                                                                                                                                                                                  && Cnorm < 10*steadyStateCondition;
  //  bool Ustalled = m_diagnostics.movingAverageHasConverged(Diagnostics::m_dUdt, m_time, m_dt)
  // Stop worrying about this for now
  bool Tstalled = false;
  bool Cstalled = false;
  bool Ustalled = false;

  bool metricConverged = false;
  // For some simulations, steady state can be defined by some other metric

  // Only test for this on level 0
  if (m_computeDiagnostics && !getCoarserLevel())
  {
    if ( (m_parameters.physicalProblem == MushyLayerParams::m_HRL
        || m_parameters.physicalProblem == MushyLayerParams::m_rayleighBenard))
    {
      metricConverged = metricConverged || m_diagnostics.movingAverageHasConverged(Diagnostics::m_Nu, m_time, m_dt);
    }

  }

  // Two ways we can 'converge':
  // 1) All fields have converged
  // 2) All fields have stalled (d/dt = constant) but metrics have converged


  bool hasConverged = false;

  bool ignoreVelocity = !solvingFullDarcyBrinkman();
  ppMain.query("ignoreVelocitySteadyState", ignoreVelocity);

  bool velConverged = ignoreVelocity || Unorm < steadyStateCondition;
  bool velStalled = ignoreVelocity || Ustalled;


  bool Tconverged = Tnorm < steadyStateCondition;
  bool Cconverged = Cnorm < steadyStateCondition;

  // Override to let us converge without salinity having actually converged
  ppMain.query("ignoreBulkConcentrationSteadyState", Cconverged);

  if (Tconverged
      && Cconverged
      && velConverged)
  {
    pout() << "Steady state reached (all fields converged)" << endl;
    hasConverged = true; // converged
  }
  else if (metricConverged)
  {
    pout() << "Steady state reached (metric converged)" << endl;
    hasConverged = true; // converged
  }
  else if (Tstalled
      && Cstalled
      && velStalled
      && metricConverged)
  {
    pout() << "Steady state reached (all fields nearly converged and metric converged)" << endl;
    hasConverged = true; // converged
  }
  else
  {
    //pout() << "Steady state not reached" << endl;
    hasConverged = false; // not converged
  }


  // One final thing to do
  // If we've converged, it may be that we just haven't kicked off the instability yet
  // Trying adding a small perturbation and keep going
  // Only do this once though
  Real maxVel = ::computeNorm(*m_vectorNew[m_advectionVel], NULL, -1, m_dx, Interval(0, SpaceDim-1), 0);
  if (hasConverged && m_doAutomaticRestart && abs(maxVel) < 1e-3)
  {
    pout() << "Max Vel = " << maxVel << ". Trying to restart with a small perturbation to kick off instability" << endl;
    addPerturbation(m_enthalpy, 1e-3);
    hasConverged = false;
    m_doAutomaticRestart = false;
  }

  if (hasConverged && m_dt < 1e-10)
  {
    pout() << "AMRLevelMushyLayer::convergedToSteadyState - all fields converged but dt < 1e-10 so keep solving" << endl;
    return false;
  }

  Real min_time = 0;
  ppMain.query("min_time", min_time);
  if (hasConverged && m_time < min_time)
  {
    pout() << "AMRLevelMushyLayer::convergedToSteadyState - converged but time < min_time (" << m_time << " < " << min_time << ")" << endl;
    return false;
  }


  return hasConverged;
}

bool AMRLevelMushyLayer::solvingFullDarcyBrinkman()
{
  if (m_parameters.darcy > 1e-10)
  {
    return true;
  }

  return false;
}

void AMRLevelMushyLayer::compute_d_dt(const int a_var, LevelData<FArrayBox>& diff, bool vector)
{



  if (vector)
  {
    diff.define(m_grids, SpaceDim);
    m_vectorNew[a_var]->copyTo(diff);

  }
  else
  {
    diff.define(m_grids, 1);
    m_scalarNew[a_var]->copyTo(diff);

  }



  DataIterator dit(m_grids);
  for (dit.reset(); dit.ok(); ++dit)
  {
    if (vector)
    {
      diff[dit] -= (*m_vectorOld[a_var])[dit];
    }
    else
    {
      diff[dit] -= (*m_scalarOld[a_var])[dit];
    }

    diff[dit] /= m_dt;
    //    diff[dit] /= max;
  }

  // Output some stuff maybe
  if (vector && a_var == m_fluidVel)
  {
    diff.copyTo(*m_vectorNew[m_dUdt]);
  }
  else if (!vector && a_var == m_enthalpy)
  {
    diff.copyTo(*m_scalarNew[m_dHdt]);
  }
  else if(!vector && a_var == m_bulkConcentration)
  {
    //    horizontallyAverage(*m_scalarNew[m_dSdt], diff);
    diff.copyTo(*m_scalarNew[m_dSdt]);
  }

}

Real AMRLevelMushyLayer::convergedToSteadyState(const int a_var, bool vector)
{

  ParmParse ppMain("main");

  LevelData<FArrayBox> diff;
  compute_d_dt(a_var, diff, vector);

  Real max;

  // Need this to ensure we only calculate sum over valid regions
  DisjointBoxLayout* fineGridsPtr = NULL;
  if (hasFinerLevel())
  {
    fineGridsPtr = &(getFinerLevel()->m_grids);
  }

  if (vector)
  {
    max = ::computeNorm(*m_vectorNew[a_var], fineGridsPtr, m_ref_ratio, m_dx, Interval(0, SpaceDim-1), 0);
  }
  else
  {
    max = ::computeNorm(*m_scalarNew[a_var], fineGridsPtr, m_ref_ratio, m_dx, Interval(0, 0), 0);
  }

  Real minAllowed = 1e-10;
  max = Max(max, minAllowed);

  // Don't consider sponge region in steady state calculations.
  //  ParmParse pp("main");
  Real spongeHeight = 0.0;
  ppMain.query("spongeHeight", spongeHeight);

  DataIterator dit(m_grids);

  if (spongeHeight > 0)
  {
    for (dit.reset(); dit.ok(); ++dit)
    {
      Box spongeBox;
      IntVect topCorner = IntVect::Zero;
      for (int dir=0; dir<SpaceDim; dir++)
      {
        if (dir == 0)
        {
          topCorner += m_numCells[0]-1 * BASISV(dir);
        }
        else
        {
          topCorner += round(m_numCells[1]*spongeHeight) * BASISV(dir);
        }
      }

      spongeBox.define(IntVect::Zero, topCorner);
      spongeBox &= diff[dit].box();
      diff[dit].setVal(0.0, spongeBox, 0);
    }
  }

  int largestDim = 0;
  if (vector)
  {
    largestDim = SpaceDim-1;
  }

  Real norm = ::computeNorm(diff, NULL, -1, m_dx, Interval(0, largestDim), m_steadyStateNormType);

  norm = norm/max;

  pout() << setiosflags(ios::scientific) << setprecision(5);

  char outString [100];
  string varName = vector ? m_vectorVarNames[a_var] : m_scalarVarNames[a_var] ;
  sprintf(outString, "d/dt (%-20s) = %e ", varName.c_str(), norm);

  if (s_verbosity > 3)
  {
    pout() << outString << endl;
  }


  return norm;
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

  ParmParse ppMain("main");

  // vel time is the time at which ths advection velocity is calculated
  // default is time centered, so at old_time + 0.5*dt
  //  ppMain.query("adv_vel_centering", advVelCentering);
  Real vel_time = old_time + advVelCentering*m_dt;
  Real half_dt = vel_time - old_time;

  LevelData<FArrayBox> porosityAv(m_grids, 1, ivGhost);
  RefCountedPtr<LevelData<FluxBox> > porosityFaceAvPtr = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids, 1, ivGhost));

  EdgeVelBCHolder edgeVelBC(m_physBCPtr->advectionVelFuncBC(m_viscousBCs));

  LevelData<FArrayBox> velOldGrown(m_grids, SpaceDim, ivGhost);
  fillVectorField(velOldGrown, old_time, m_fluidVel, true);

  Real advPorosityLimit = m_solidPorosity;
  ppMain.query("advPorosityLimit", advPorosityLimit);

  // Get initial guess at advection velocity from U^n
  bool useOldAdvVel = false;
  ppMain.query("useOldAdvVelForTracing", useOldAdvVel);
  if (useOldAdvVel)
  {
    // do nothing - we already have m_advVel from previous timestep
  }
  else
  {
    CellToEdge(velOldGrown, m_advVel);
  }

  edgeVelBC.applyBCs(m_advVel, m_grids, m_problem_domain, m_dx,
                     false); // inhomogeneous

  if (m_doEulerPart)
  {
    if (s_verbosity >= 5)
    {
      pout() << "AMRLevelMushyLayer::computeAdvectionVelocities - doing advection of velocities" << endl;
    }

    if (m_implicitAdvectionSolve)
    {

      bool doFRupdates = false;
      bool doProjection = false; // don't do CC projection
      bool compute_uDelU = true;
      bool MACProjection = true;

      ppMain.query("usePhiForImplicitAdvectionSolve", MACProjection);

      // do full timestep - this should now give the same velocity as the later solve?
      advectionSourceTerm.exchange();
      computeCCvelocity(advectionSourceTerm, old_time, m_dt, doFRupdates, doProjection,
                        compute_uDelU, MACProjection);

      m_vectorNew[m_viscousSolveSrc]->copyTo(*m_vectorNew[m_advectionImplicitSrc]);

      // New version:
      if (MACProjection)
      {
        fillVectorField(*m_vectorNew[m_advUpreProjection], vel_time, m_advUpreProjection, false, true);

      }
      else
      {
        for (DataIterator dit = m_vectorNew[m_advUpreProjection]->dataIterator(); dit.ok(); ++dit)
        {
          (*m_vectorNew[m_advUpreProjection])[dit].copy((*m_vectorNew[m_UpreProjection])[dit]);
        }
        fillVectorField(*m_vectorNew[m_advUpreProjection], vel_time, m_UpreProjection, false, true);

      }


      CellToEdge(*m_vectorNew[m_advUpreProjection], m_advVel);


    }
    else
    {
      if (s_verbosity >= 5)
      {
        pout() << "AMRLevelMushyLayer::computeAdvectionVelocities() - explicit tracing scheme" << endl;
      }

      // Explicit advection solve (tracing)

      int saveAdvectionMethod = m_advectionMethod;

      // Option to do a different method here
      int velPredMethod = m_advectionMethod;
      ppMain.query("velPredictionMethod", velPredMethod);
      m_advectionMethod = velPredMethod; //m_noPorosity;

      fillScalarFace(*porosityFaceAvPtr, vel_time, m_porosity, true);
      fillScalars(porosityAv, vel_time, m_porosity, true);

      // Construct the thing to advect
      LevelData<FArrayBox> porosityGrown(m_grids, 1, ivGhost);
      LevelData<FArrayBox> U_to_advect(m_grids, SpaceDim, ivGhost);
      fillScalars(porosityGrown, old_time, m_porosity, true);

      LevelData<FArrayBox> ccAdvVel(m_grids, SpaceDim, ivGhost);
      fillVectorField(ccAdvVel, old_time, m_fluidVel, true);

      setVelZero(ccAdvVel, advPorosityLimit);

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
        setVelZero(U_to_advect, advPorosityLimit);

        // Divide each component of velocity by porosity if we're advecting u/chi
        if ( m_advectionMethod == m_porosityInAdvection)
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
        if  (m_advectionMethod == m_porosityOutsideAdvection
            || m_advectionMethod == m_porosityInAdvection )
        {

          setVelZero(m_advVel, advPorosityLimit);

          m_advVel[dit].divide((*porosityFaceAvPtr)[dit], m_advVel[dit].box(), 0, 0);
          for (int dir = 0; dir < SpaceDim; dir++)
          {
            ccAdvVel[dit].divide(porosityAv[dit], ccAdvVel[dit].box(), 0, dir);
          }

          // this does BCs (doesn't fill interior)
          fillVectorField(ccAdvVel, old_time, m_fluidVel, false);
        }

        if ( m_advectionMethod == m_noPorosity)
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

      EdgeToCell(m_advVel, *m_vectorNew[m_advectionVel]);
      Real maxAdvU = ::computeNorm(*m_vectorNew[m_advectionVel], NULL, 1, m_dx, Interval(0,0));

      if (s_verbosity >= 1 && maxAdvU > 1e10)
      {
        pout() << "WARNING - max advection velocity (pre-projection) = " << maxAdvU << endl;
      }


      // Need to recover the advection velocity, as
      // traceAdvectionVel will replce m_advVel with the upwinded U_to_advect
      if ( m_advectionMethod == m_porosityInAdvection)
      {
        for (dit.reset(); dit.ok(); ++dit)
        {
          m_advVel[dit].mult((*porosityFaceAvPtr)[dit], m_advVel[dit].box(), 0, 0);
        }
      }


      m_advectionMethod = saveAdvectionMethod;

      // DO some post tracing smoothing in the x-direction
      // this is necessary to remove some \delta x scale instabilities which seem to arise

      m_advVel.exchange();

      horizontallySmooth(m_advVel);

      m_advVel.exchange();


    } // end if implicit/explicit solve

    setVelZero(m_advVel, advPorosityLimit);

    //m_advVel contains predicted face centred velocity at half time step
    //Now need to project it and do lambda corrections etc
    Real project_dt = 2*half_dt; // this usually just works out at m_dt, for time centred advection velocities
    if (m_implicitAdvectionSolve)
    {
      project_dt = m_dt;
    }

    correctEdgeCentredVelocity(m_advVel, project_dt);

    // Maybe replace m_advVel with time independent version in low porosity regions?
    Real chiLimit = 0.0;
    ppMain.query("porousAdvVelLimit", chiLimit);
    if (chiLimit > 0.0)
    {
      LevelData<FluxBox> mushyAdvVel(m_advVel.disjointBoxLayout(), 1, m_advVel.ghostVect());
      //        fillUnprojectedDarcyVelocity(mushyAdvVel, m_time-m_dt);
      calculateTimeIndAdvectionVel(m_time-m_dt, mushyAdvVel);

      int count = 0;
      for (DataIterator dit = m_advVel.dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& chi = (*m_scalarNew[m_porosity])[dit];
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
            if (chi(iv) < chiLimit)
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

    EdgeToCell(m_advVel, *m_vectorNew[m_advUstar]);
    correctEdgeCentredVelocity(m_advVel, m_dt);
  }


  bool enforceAnalyticVel = false;

  ppMain.query("analyticVel", enforceAnalyticVel);

  if (enforceAnalyticVel)
  {
    fillAnalyticVel(m_advVel);

    LevelData<FluxBox> velDiff(m_advVel.disjointBoxLayout(), 1, m_advVel.ghostVect());

    for(DataIterator dit = velDiff.dataIterator(); dit.ok(); ++dit)
    {
      velDiff[dit].copy(m_advVel[dit]);
    }

    bool correctVel = false;
    ppMain.query("correctAnalyticVel", correctVel);
    if (correctVel)
    {
      correctEdgeCentredVelocity(m_advVel, m_dt);
    }


    for(DataIterator dit = velDiff.dataIterator(); dit.ok(); ++dit)
    {
      velDiff[dit].minus(m_advVel[dit], m_advVel[dit].box(), 0,0,1);
    }


  } // end if enforce analytic vel

  Divergence::levelDivergenceMAC(*m_scalarNew[m_divUadv], m_advVel, m_dx);

  // Finally
  EdgeToCell(m_advVel, *m_vectorNew[m_advectionVel]);

}

void AMRLevelMushyLayer::horizontallySmooth(LevelData<FluxBox>& a_flux)
{
  ParmParse ppMain("main");
  DataIterator dit = a_flux.dataIterator();

  int dir = 0; // horizontally averaging

  Real alpha = 0.0; // default is no smoothing

  ppMain.query("postTraceSmoothing", alpha);

  if (alpha==0.0)
  {
    return;
  }

  // todo - write post trace smoothing in fortran
  for (dit.reset(); dit.ok(); ++dit)
  {
    for (int velDir=0; velDir < SpaceDim; velDir++)
    {
      FArrayBox& vel = a_flux[dit][velDir];
      Box b = vel.box();

      for (BoxIterator bit(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        IntVect ivUp = iv+BASISV(dir);

        if (b.contains(ivUp))
        {

          Real neighbour = vel(ivUp);

          // Make sure we don't change the velocity if the neighboroung value is much larger
          // this will prevent changing zero velocity cells, or including NaN type values
          if( abs(neighbour) < 1e100 ) //abs(vel(iv)) > 1e-15 &&
          {
            vel(iv) = (1-alpha)*vel(iv)+alpha*neighbour;
          }

        }
      }
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

  EdgeVelBCHolder edgeVelBC(m_physBCPtr->advectionVelFuncBC(m_viscousBCs));
  edgeVelBC.applyBCs(a_advVel, m_grids, m_problem_domain, m_dx,
                     false); // inhomogeneous
  if (s_verbosity >= 6)
  {
    pout() << "  AMRLevelMushyLayer::correctEdgeCentredVelocity() - applied init BCs" << endl;
  }

  if(m_scaleP_MAC)
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

  ParmParse ppMain("main");

  a_advVel.exchange();
  EdgeToCell(a_advVel, *m_vectorNew[m_advUstar]);

  if (m_doProjection)
  {
    if (s_verbosity >= 6)
    {
      pout() << "  AMRLevelMushyLayer::correctEdgeCentredVelocity() - do projection" << endl;
    }

    int projNum = 0;
    int maxNumProj = 1;
    ppMain.query("max_num_MAC_proj", maxNumProj);

    //    Divergence::levelDivergenceMAC(*m_scalarNew[m_divUadv], a_advVel, m_dx);
    Real maxDivU = 10*m_maxDivUFace ; //::computeNorm(*m_scalarNew[m_divUadv], NULL, 1, m_dx, Interval(0,0));

    while ( (maxDivU > m_maxDivUFace && projNum < maxNumProj) ||
        (ppMain.contains("analyticVel") && projNum < 1) ) // do at least one projection of analytic vel
    {
      projNum++;

      int exitStatus = m_projection.levelMacProject(a_advVel, old_time, a_dt, pressureScalePtr, crsePressureScalePtr,
                                   pressureScaleEdgePtrOneGhost, crsePressureScaleEdgePtr);


      Divergence::levelDivergenceMAC(*m_scalarNew[m_divUadv], a_advVel, m_dx);

      maxDivU = ::computeNorm(*m_scalarNew[m_divUadv], NULL, 1, m_dx, Interval(0,0), 0);
      pout() << "  MAC Projection (#" << projNum << " on level "<< m_level << "), exit status = " << exitStatus << ", max(div u) = " << maxDivU << endl;


    }

  }

  //  Divergence::levelDivergenceMAC(*m_scalarNew[m_divUadv], a_advVel, m_dx);

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

void AMRLevelMushyLayer::copyNewToOldStates()
{
  if (s_verbosity > 4)
  {
    pout() << "AMRLevelMushyLayer::copyNewToOldStates(" << m_time << " to " << m_time-m_dt << ") " << endl;
  }

  // Before we do this, populate d(porosity)/dt, which we need for computing advection update source terms
  // We actually don't use d(porosity)/dt any more, so stop doing this
  /*
  LevelData<FArrayBox> porosity_new(m_dPorosity_dt.disjointBoxLayout(), 1, m_dPorosity_dt.ghostVect());
  LevelData<FArrayBox> porosity_old(m_dPorosity_dt.disjointBoxLayout(), 1, m_dPorosity_dt.ghostVect());

  fillScalars(porosity_new, m_time, m_porosity, true, true);
  fillScalars(porosity_old, m_time-m_dt, m_porosity, true, true);

  for (DataIterator dit = m_dPorosity_dt.dataIterator(); dit.ok(); ++dit)
  {
    m_dPorosity_dt[dit].copy(porosity_new[dit]);
    m_dPorosity_dt[dit].minus(porosity_old[dit]);
    m_dPorosity_dt[dit].mult(1/m_dt);
  }
  */

  // Copy the new to the old
  // Old now contains values at n, new will contain values at n+1 eventually
  for (int a_scalarVar = 0; a_scalarVar < m_numScalarVars; a_scalarVar++)
  {
    m_scalarNew[a_scalarVar]->copyTo(m_scalarNew[a_scalarVar]->interval(),
                                     *m_scalarOld[a_scalarVar],
                                     m_scalarOld[a_scalarVar]->interval());

  }
  for (int vectorVar = 0; vectorVar < m_numVectorVars; vectorVar++)
  {
    m_vectorNew[vectorVar]->copyTo(m_vectorNew[vectorVar]->interval(),
                                   *m_vectorOld[vectorVar],
                                   m_vectorOld[vectorVar]->interval());
  }

  m_advVelNew.copyTo(m_advVelOld);

}




void AMRLevelMushyLayer::fillTCl(LevelData<FArrayBox>& a_phi, Real a_time,
                                 bool doInterior, bool quadInterp )
{
  fillMultiComp(a_phi, a_time, m_temperature, m_liquidConcentration, doInterior, quadInterp);
}


void AMRLevelMushyLayer::fillMultiComp(LevelData<FArrayBox>& a_phi, Real a_time, int scal1, int scal2,
                                       bool doInterior , bool quadInterp)
{
  CH_assert(a_phi.nComp() == 2);

  LevelData<FArrayBox> temp1(a_phi.disjointBoxLayout(), 1, a_phi.ghostVect());
  LevelData<FArrayBox> temp2(a_phi.disjointBoxLayout(), 1, a_phi.ghostVect());
  fillScalars(temp1, a_time, scal1, doInterior, quadInterp);
  fillScalars(temp2, a_time, scal2, doInterior, quadInterp);

  // This doesn't copy ghost cells!
  //    temp1.copyTo(Interval(0,0), a_phi, Interval(0,0));
  //    temp2.copyTo(Interval(0,0), a_phi, Interval(1, 1));

  // Need to do copying this way to transfer ghost cells
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
  {
    a_phi[dit].copy(temp1[dit], 0, 0, 1);
    a_phi[dit].copy(temp2[dit], 0, 1, 1);
  }
}


void AMRLevelMushyLayer::fillHC(LevelData<FArrayBox>& a_phi, Real a_time,
                                bool doInterior , bool quadInterp )
{
  fillMultiComp(a_phi, a_time, m_enthalpy, m_bulkConcentration, doInterior, quadInterp);
}
/*******/
Real AMRLevelMushyLayer::advance()
{


  ParmParse ppParams("parameters");
  ParmParse ppMain("main");

  if (m_dtReduction > 0)
  {
    //		m_dt = m_dt * m_dtReduction;
    m_dtReduction = -1;
  }


  AMRLevelMushyLayer *amrMLcrse = NULL;
  if (m_level > 0)
  {
    amrMLcrse = getCoarserLevel();
  }

  // Do some things on only the coarse level
  if (amrMLcrse == NULL)
  {
    // 1) set all levels reflux registers to zero
    AMRLevelMushyLayer *amrMLptr = this;

    // Loop over all levels, starting from 0
    Real oldLev0Dt = m_dt;
    Real oldTime = m_time;


    while(amrMLptr)
    {

      if (amrMLptr->hasFinerLevel())
      {
        amrMLptr->setFluxRegistersZero();
      }

      // 2) Set timestep failed flags to false, after halving dts
      if (amrMLptr->m_timestepFailed)
      {
        //        pout() << "PREVIOUS TIMESTEP FAILED - HALVING DT" << endl;
        Real prevDt = amrMLptr->dt();
        Real prevTime = amrMLptr->time();
        amrMLptr->dt(prevDt/2); // halve the timestep
        amrMLptr->time(oldTime - oldLev0Dt);

        pout() << "Reset on level " << amrMLptr->level() << " from (time=" << prevTime << ", dt = " << prevDt << ")";
        pout() << " to (time=" << amrMLptr->time() << ", dt = " << amrMLptr->dt() << ")" << endl;
      }
      amrMLptr->m_timestepFailed = false;


      // Check if we want to add a perturbation
      if (m_time < m_perturbationTime && (m_time + m_dt) > m_perturbationTime)
      {
        //pout() << "Adding perturbation with wavenumber " << m_perturbationWavenumber << endl;
        Real phaseShift = 0;
        ppMain.query("perturbationPhaseShift", phaseShift);
        addPerturbation(m_enthalpy, m_delayedPerturbation, m_perturbationWavenumber, phaseShift);
      }


      amrMLptr = amrMLptr->getFinerLevel();
    }

  }

  if (s_verbosity >= 1) // && m_level == 0
  {
    // Let's determine the CFL number we're running at
    Real maxAdvU = getMaxVelocity();
    m_computedCFL = m_dt*maxAdvU/m_dx;

    pout() << " AMRLevelMushyLayer::advance (level = "<< m_level << ", time=" << m_time << ", dt = " << m_dt << ", CFL=" << m_computedCFL << ")" << endl;
  }

  // Reset BCs if they change with time
  m_parameters.setTime(m_time); // BCs are stored in m_parameters
  setAdvectionBCs(); // Reset BCs on advection physics objects
  // Elliptic operators get the BC from m_parameters when they're defined (later)


  Real rampBuoyancy = 0.0;
  ppParams.query("rampBuoyancy", rampBuoyancy);
  Real maxRaC = 0.0;
  Real maxRaT = 0.0;

  ppParams.query("maxRaT", maxRaT);
  ppParams.query("maxRaC", maxRaC);

  Real initRaC = 0;
  Real initRaT = 0;
  ppParams.query("initRaT", initRaT);
  ppParams.query("initRaC", initRaC);
  if (rampBuoyancy > 0)
  {
    if (m_parameters.rayleighComposition == 0)
    {
      m_parameters.rayleighComposition = initRaC;
    }
    if (m_parameters.rayleighTemp == 0)
    {
      m_parameters.rayleighTemp = initRaT;
    }


    m_parameters.rayleighComposition = m_parameters.rayleighComposition *rampBuoyancy;
    m_parameters.rayleighTemp = m_parameters.rayleighTemp *rampBuoyancy;

    m_parameters.rayleighComposition = min(maxRaC, m_parameters.rayleighComposition);
    m_parameters.rayleighTemp = min(maxRaT, m_parameters.rayleighTemp);

    pout() << "RaC = " << m_parameters.rayleighComposition  << ", RaT = " <<  m_parameters.rayleighTemp  << endl;
  }

  // Updates 'new' variables
  updateEnthalpyVariables();

#ifdef CH_MPI
  {
      CH_TIME("MPI_Barrier exchange");
      MPI_Barrier(Chombo_MPI::comm);
  }
#endif


  Real old_time = m_time;
  Real new_time = m_time + m_dt;
  //  Real half_time = new_time - m_dt / 2;

  m_time = new_time;



  // Move 'new' variables to 'old' variables
  // do this after we've incremented m_time for consistency with coarser levels

  // Make sure we fill ghost cells before computing d(porosity)/dt
  LevelData<FArrayBox> oldPorosity(m_dPorosity_dt.disjointBoxLayout(), 1, m_dPorosity_dt.ghostVect());
  LevelData<FArrayBox> newPorosity(m_dPorosity_dt.disjointBoxLayout(), 1, m_dPorosity_dt.ghostVect());
  fillScalars(oldPorosity, m_time-m_dt, m_porosity, true, true);
  fillScalars(newPorosity, m_time, m_porosity, true, true);

  for (DataIterator dit = m_dPorosity_dt.dataIterator(); dit.ok(); ++dit)
  {
    m_dPorosity_dt[dit].copy(newPorosity[dit]);
    m_dPorosity_dt[dit].minus(oldPorosity[dit]);
    m_dPorosity_dt[dit].divide(m_dt);
  }
  copyNewToOldStates();

  //Gradually ramp up temperature forcing
  if (m_parameters.physicalProblem == MushyLayerParams::m_poiseuilleFlow)
  {
    stokesDarcyForcing(*m_scalarNew[m_temperature], new_time);
    stokesDarcyForcing(*m_scalarOld[m_temperature], old_time);
  }

  // Some useful things to have around
  DataIterator dit(m_grids);
  IntVect ivGhost = m_numGhost * IntVect::Unit;
  IntVect advectionGhost = m_numGhostAdvection * IntVect::Unit;

  /*
   * First step for momentum equation: predict face centred (u/chi)^{n+1/2} using (u/chi)^n
   */

  LevelData<FArrayBox> advectionSourceTerm(m_grids, SpaceDim, advectionGhost);

  calculatePermeability(); //make sure this is up to date



  LevelData<FArrayBox> zeroSrc(m_grids, 1, ivGhost);
  setValLevel(zeroSrc, 0.0);

  // Need to redefine solvers if variables have changed
  defineSolvers(m_time-m_dt); // define at old time

  if (solvingFullDarcyBrinkman())
  {
    // This fills all the ghost cells of advectionSourceTerm
     computeAdvectionVelSourceTerm(advectionSourceTerm);

    //    Real vel_centering = 0.5;
    //    ppMain.query("adv_vel_centering", vel_centering);
    computeAdvectionVelocities(advectionSourceTerm, m_adv_vel_centering);
  }
  else
  {
    //    this->finestLevel()
    //Calculate time centred advection velocity
    calculateTimeIndAdvectionVel(m_time-m_dt, m_advVel);
  }

  // Check we're still satisfying CFL condition
  // If not, skip scalar advection for this step
  bool doAdvectiveSrc = true;

  if (!currentCFLIsSafe(true))
  {
    doAdvectiveSrc = false;
  }

  // Another sanity check
//  Divergence::levelDivergenceMAC(*m_scalarNew[m_divUadv], m_advVel, m_dx);
//  Real  maxDivU = ::computeNorm(*m_scalarNew[m_divUadv], NULL, 1, m_dx, Interval(0,0));
//  pout() << "  Sanity check: max(div u) = " << maxDivU << endl;

  // always* advect lambda (and update flux registers)
  // do this as soon as we have advection velocities, in case we want to
  // correct them prior to HC advection
  // * skip if we're danger of violating the CFL condition
  if (doAdvectiveSrc)
  {
    advectLambda(true);
  }

  // Scalar advection diffusion
  bool doScalarAdvectionDiffusion = true;
  ppMain.query("doScalarAdvectionDiffusion", doScalarAdvectionDiffusion);

  bool skipNewLevelScalars = false;
  ppMain.query("skipNewLevelScalars", skipNewLevelScalars);
  if (m_newLevel && m_level > 0)
  {
    if (skipNewLevelScalars)
    {
      pout() << "First time step on a new level - skipping advection and diffusion" << endl;
      doScalarAdvectionDiffusion = false;
    }
    m_newLevel = false;
  }

  if (!(m_parameters.physicalProblem == MushyLayerParams::m_poiseuilleFlow ||
      m_parameters.physicalProblem == MushyLayerParams::m_soluteFluxTest ||
      m_parameters.physicalProblem == MushyLayerParams::m_zeroPorosityTest) && doScalarAdvectionDiffusion)
  {
    // Computation of advection velocity is done
    // Now do advection and diffusion of scalar fields
    CH_TIME("AMRLevelMushyLayer::advection-diffusion");
    if (s_verbosity >= 5)
    {
      pout() << "AMRLevelMushyLayer::advance - advect and diffuse scalars" << endl;
    }

    int exitStatus = 0;

    //bool doFRupdates = true;


    // Need to construct multi component object
    LevelData<FArrayBox> HC_new(m_grids, 2, IntVect::Unit);
    LevelData<FArrayBox> HC_old(m_grids, 2, IntVect::Unit);
    LevelData<FArrayBox> srcMultiComp(m_grids, 2, IntVect::Zero);

    fillHC(HC_new, m_time);
    fillHC(HC_old, m_time-m_dt);

    setValLevel(srcMultiComp, 0.0);

    bool doFRUpdates = true;

    exitStatus = multiCompAdvectDiffuse(HC_old, HC_new, srcMultiComp, doFRUpdates, doAdvectiveSrc);

    bool solverFailed = (exitStatus == 2 || exitStatus == 4 || exitStatus == 6);
//    bool solveSuccess = !solverFailed;

    //    if (solveSuccess)
    //    {
    // Get back the answer if solver was a success
    HC_new.copyTo(Interval(0,0), *m_scalarNew[m_enthalpy], Interval(0,0));
    HC_new.copyTo(Interval(1,1), *m_scalarNew[m_bulkConcentration], Interval(0,0));

    updateEnthalpyVariables();
    //    }

    if (solverFailed)
    {
      if (m_ignoreSolverFails)
      {
        pout() << "Ignoring all solver fails." << endl;
      }

      else
      {
        // Alternative way of restarting - set this to true to just do this timestep again
        // with half the dt. The trouble is we still output a bad file from this timestep which is annoying
        m_timestepFailed = true;

        Vector<string> failedReasons(9);

        failedReasons[4] = "Solver hang";
        failedReasons[8] = "Norm not reduced enough";
        failedReasons[2] = "Reached iter max";
        failedReasons[1] = "Initial norm not reduced enough";

        pout() << "Solver failed. Exit status: " << exitStatus << endl;

      }


    } // end if scalar diffusion solver failed

  } // end if doing scalar advection/diffusion



  if (solvingFullDarcyBrinkman())
  {
    // If we're skipping advective srcs for this timestep, skip this too
    if (!doAdvectiveSrc)
    {
      DataIterator dit = m_vectorNew[m_UdelU]->dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
      {
        (*m_vectorNew[m_UdelU])[dit].setVal(0.0);
      }
    }

    bool doFRupdates = true;
    bool compute_uDelU = !m_implicitAdvectionSolve && doAdvectiveSrc;
    bool doProjection = true;
    computeCCvelocity(advectionSourceTerm, m_time-m_dt, m_dt, doFRupdates, doProjection, compute_uDelU);
  }


  getExtraPlotFields();

  //supposed to return dt but the return value is never used, so don't.
  //if the return value is ever used we will know because this will break it
  return -1;
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

  VelBCHolder velBC(m_physBCPtr->uStarFuncBC(m_viscousBCs));
  RefCountedPtr<LevelData<FluxBox> > crsePressureScaleEdgePtr, pressureScaleEdgePtr;
  RefCountedPtr<LevelData<FArrayBox> > crsePressureScalePtr, pressureScalePtr;
  LevelData<FArrayBox> *crseVelPtr = NULL;

  AMRLevelMushyLayer* amrMLcrse = getCoarserLevel();

  DataIterator dit = m_scalarNew[m_porosity]->dataIterator();

  Real half_time = a_oldTime + a_dt/2;
  Real new_time = a_oldTime + a_dt; //m_time;

  // We seem to need to make this larger as the resolution decreases.
  // I don't know the actual scaling, so this is just a guess.
  Real ccVelPorosityLimit = m_solidPorosity; //1e-4;
  //  m_lowerPorosityLimit = 1e-10;
  ParmParse pp("main");
  pp.query("ccvel_porosity_cap", ccVelPorosityLimit);



  // Calculate Ustar
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMushyLayer::advance - computeUstar (level " << m_level << ")" << endl;
  }


  computeUstar(*m_vectorNew[m_UdelU], advectionSourceTerm, a_oldTime, a_dt, doFRupdates, a_MACprojection, compute_uDelU);

  Real maxUstar = ::computeNorm(*m_vectorNew[ustarVar], NULL, 1, m_dx, Interval(0,0), 0);
  if (s_verbosity >= 1 && maxUstar > 1e10)
  {
    pout() << "WARNING - max cell-centred velocity (pre-projection) = " << maxUstar << endl;
  }

  // Don't think this really helps
  Real porositySmoothing = 0.001;
  porositySmoothing = 0.0;

  if (m_scaleP_CC)
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

    if (m_scaleP_CC)
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
    //                      FArrayBox& U = (*m_vectorNew[m_fluidVel])[dit];
    //                      int temp=0;
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


  if (m_doProjection && a_doProjection)
    {


    CH_TIME("AMRLevelMushyLayer::levelProject");

    // Repeatedly apply projection - check divergence goes to 0
    //    int numProj = 2;

//    Real initMaxDivU = 1e300;
    QuadCFInterp interp;
    Divergence::levelDivergenceCC(*m_scalarNew[m_divU], *m_vectorNew[uvar], NULL, m_dx, true, interp);

    Real initMaxDivU = computeNorm(*m_scalarNew[m_divU], NULL, -1, m_dx, Interval(0,0), 0);
    if (s_verbosity >= 3)
    {
      pout() << "  CCProjection: init max(div(U)) = " << initMaxDivU << endl;
    }

    Real maxDivU = initMaxDivU;
    Real relTol = 1e-5;
    int i = 0;

    int maxNumProj = 1;
    pp.query("maxProjections", maxNumProj);

    //    for (int i = 0; i < numProj; i++)
    while(! (maxDivU < initMaxDivU*relTol || maxDivU < relTol || i >= maxNumProj))
    {

      setVelZero(*m_vectorNew[uvar], ccVelPorosityLimit);

      velBC.applyBCs(*m_vectorNew[uvar], m_grids, m_problem_domain,
                           m_dx, false); // inhomogeneous

      if (i > 0)
      {
        // todo - WARNING - THIS WILL BREAK AMR AS IT SETS THE BC WRONG FOR MAC SOLVE

        // this only works on level 0 - we don't have any CF boundaries
        if (m_level==0)
        {
        m_projection.AdditionalLevelProject(*m_vectorNew[uvar], crseVelPtr,
                                       new_time, a_dt,  m_isViscous);
        }
      }
      else
      {
      m_projection.LevelProject(*m_vectorNew[uvar], crseVelPtr,
                                new_time, a_dt, pressureScalePtr, crsePressureScalePtr, pressureScaleEdgePtr, crsePressureScaleEdgePtr,
                                m_isViscous);
      }
      // as things stand now, physical BC's are re-set in LevelProjection

      //Get a copy of the pressure so we can write it out
      m_projection.unscaledPi(*m_scalarNew[m_pressure], a_dt);

      // need to do physical boundary conditions and exchanges
      velBC.applyBCs(*m_vectorNew[uvar], m_grids, m_problem_domain,
                     m_dx, false); // inhomogeneous
      m_vectorNew[uvar]->exchange();

      QuadCFInterp interp;
      Divergence::levelDivergenceCC(*m_scalarNew[m_divU], *m_vectorNew[uvar], NULL, m_dx, true, interp);

      // Actually want norm (sum of absolute values)
      maxDivU = computeNorm(*m_scalarNew[m_divU], NULL, -1, m_dx, Interval(0,0), 0);

      pout() << "  CCProjection: max(div(U)) = " << maxDivU << endl;

      if (i == 0)
      {
        initMaxDivU = maxDivU;
      }

      i++;
    }



  } // end if doing projection


  // Make sure zero porosity regions have no velocity
  setVelZero(*m_vectorNew[uvar], ccVelPorosityLimit);

  //Testing:
//  LevelData<FArrayBox> U_chi(m_grids, SpaceDim);
//  fillVectorField(U_chi, m_time, m_U_porosity, true, true);
//  Real max = ::computeMax(U_chi, NULL, 1, Interval(0, SpaceDim-1));
//  if (max > 1e5)
//  {
//    pout() << "U_chi > 1e5!!! " << endl;
//  }


//  QuadCFInterp interp;
//  Divergence::levelDivergenceCC(*m_scalarNew[m_divU], *m_vectorNew[uvar], NULL, m_dx, true, interp);
//  Real maxDivU = computeNorm(*m_scalarNew[m_divU], NULL, -1, m_dx, Interval(0,0), 0);

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

void AMRLevelMushyLayer::fillPressureSrcTerm(LevelData<FArrayBox>& gradP,
                                             LevelData<FArrayBox>& pressureScale,
                                             Real a_time,
                                             bool a_MACprojection)
{

  computeGradP(gradP, a_time, a_MACprojection);

  DataIterator dit = gradP.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    if ((a_MACprojection && m_scaleP_MAC) || (!a_MACprojection && m_scaleP_CC))
    {
      for (int dir=0;dir<SpaceDim;dir++)
      {
        gradP[dit].mult(pressureScale[dit], 0, dir);
      }
    }
  }

}


void AMRLevelMushyLayer::setVelZero(FArrayBox& a_vel, const FArrayBox& a_porosity, const Real a_limit, const int a_radius)
{

  if (a_radius > 0)
  {

    IntVectSet porousCells;

    for (BoxIterator bit(a_porosity.box()); bit.ok(); ++bit)
    {
      IntVect iv = bit();

      if (a_porosity(iv) <= a_limit)
      {
        porousCells |= iv;
      }
    }

    porousCells.grow(a_radius);

    IntVectSet setZeroCells = porousCells;
    setZeroCells &= a_vel.box();

    for (IVSIterator its(setZeroCells); its.ok(); ++its)
    {
      for (int comp=0; comp<a_vel.nComp(); comp++)
      {
        a_vel(its(), comp) = 0;
      }
    }



  }
  else
  {

    Box region = a_vel.box();
    region &= a_porosity.box();

    FORT_SETZEROVELOCITY( CHF_FRA(a_vel),
                          CHF_CONST_FRA(a_porosity),
                          CHF_BOX(region),
                          CHF_CONST_REAL(a_limit));

  }
}

void AMRLevelMushyLayer::setCCVelZero(Real a_limit)
{
  setVelZero(*m_scalarNew[m_fluidVel], a_limit);
}
void AMRLevelMushyLayer::setVelZero(LevelData<FArrayBox>& a_vel, Real a_limit, int a_radius)
{
  if (a_limit < 0)
  {
    a_limit = m_solidPorosity;
  }

  for (DataIterator dit = a_vel.dataIterator(); dit.ok(); ++dit)
  {
    setVelZero(a_vel[dit], (*m_scalarNew[m_porosity])[dit], a_limit, a_radius);
  }
}

void AMRLevelMushyLayer::setVelZero(LevelData<FluxBox>& a_vel, Real a_limit)
{
  if (a_limit < 0)
  {
    a_limit = m_solidPorosity;
  }

  // Get the porosity
  LevelData<FluxBox> porosity(a_vel.disjointBoxLayout(), 1, a_vel.ghostVect());
  fillScalarFace(porosity, m_time, m_porosity, true);

  for (DataIterator dit = a_vel.dataIterator(); dit.ok(); ++dit)
  {

    for (int dir=0; dir<SpaceDim; dir++)
    {
      setVelZero(a_vel[dit][dir], porosity[dit][dir], a_limit);

    }
  }
}



void AMRLevelMushyLayer::computeScalDiffusion(LevelData<FArrayBox>& a_src, int a_var)
{

  BCHolder bc;
  getScalarBCs(bc, a_var, false);

  AMRPoissonOpFactory* op = new AMRPoissonOpFactory();
  op->define(m_problem_domain, m_grids, m_dx, bc);

  RefCountedPtr<AMRPoissonOpFactory> OpFact = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(op);

  RefCountedPtr<AMRPoissonOp> amrpop =  RefCountedPtr<AMRPoissonOp>(
      (AMRNonLinearMultiCompOp*) OpFact->AMRnewOp(m_problem_domain) );

  LevelData<FArrayBox> *crseVar = NULL;

  AMRLevelMushyLayer* crseML = getCoarserLevel();
  if (crseML)
  {
    crseVar = &(*crseML->m_scalarNew[a_var]);
  }

  amrpop->setAlphaAndBeta(0, 1);

  // This just calls applyOpI if crseHC = NULL, else does CF interpolation
  amrpop->applyOpMg(a_src, *m_scalarNew[a_var], crseVar, false);


}

void AMRLevelMushyLayer::advectLambda(bool doFRupdates)
{
  CH_TIME("AMRLevelMushyLayer::advectLambda");

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::advectLambda" << endl;
  }

  ParmParse ppMain("main");

  m_advVel.exchange();

  // Do lambda advection
  m_scalarNew[m_lambda]->copyTo(Interval(0,0), *m_scalarOld[m_lambda], Interval(0,0));

  // Want to get the lambda flux back so we can remove it later
  LevelData<FluxBox> lambdaFlux(m_grids, 1);
  advectScalar(m_lambda, m_lambda, m_advVel, true, lambdaFlux); // advect without a diffusive source term

  setValLevel(*m_vectorNew[m_advVelCorr], 0.0);

}


void AMRLevelMushyLayer::updateEnthalpyVariables()
{


  CH_TIME("AMRLevelMushyLayer::updateEnthalpyVariables");

  // Apply BCs
  fillScalars(*m_scalarNew[m_bulkConcentration], m_time, m_bulkConcentration);
  fillScalars(*m_scalarNew[m_enthalpy], m_time, m_enthalpy);

  ::updateEnthalpyVariables(*m_scalarNew[m_enthalpy], *m_scalarNew[m_bulkConcentration],
                            *m_scalarNew[m_temperature], *m_scalarNew[m_liquidConcentration], *m_scalarNew[m_solidConcentration],
                            *m_scalarNew[m_porosity],
                            *m_scalarNew[m_enthalpySolidus],*m_scalarNew[m_enthalpyLiquidus],*m_scalarNew[m_enthalpyEutectic],
                            m_parameters);


  doRegularisationOps(*m_scalarNew[m_liquidConcentration], m_liquidConcentration);
  doRegularisationOps(*m_scalarNew[m_porosity], m_porosity);

  // A few alterations for test problems
  if (m_parameters.physicalProblem == m_parameters.m_poiseuilleFlow)
  {
    initialDataPoiseuille();
  }
  else if(m_parameters.physicalProblem == m_parameters.m_convectionMixedPorous ||
      m_parameters.physicalProblem == m_parameters.m_zeroPorosityTest)
  {
    fillFixedPorosity(*m_scalarNew[m_porosity]);
    m_scalarNew[m_porosity]->copyTo(*m_scalarOld[m_porosity]);
  }


  computeLambdaPorosity();

}

void AMRLevelMushyLayer::computeLambdaPorosity()
{
  for (DataIterator dit = m_scalarNew[m_lambda]->dataIterator(); dit.ok(); ++dit)
  {
    (*m_scalarNew[m_lambda_porosity])[dit].copy((*m_scalarNew[m_lambda])[dit]);
    (*m_scalarNew[m_lambda_porosity])[dit].divide((*m_scalarNew[m_porosity])[dit]);

    (*m_scalarOld[m_lambda_porosity])[dit].copy((*m_scalarOld[m_lambda])[dit]);
    (*m_scalarOld[m_lambda_porosity])[dit].divide((*m_scalarOld[m_porosity])[dit]);
  }
}

void AMRLevelMushyLayer::fillFixedPorosity(LevelData<FArrayBox>& a_porosity)
{
  DataIterator dit = a_porosity.dataIterator();

  ParmParse pp("main");

  Real stdev = 0.005;
  Real maxChi = 1.05; // want to make this a bit more than 1, so we get porosity=1 region of finite size
  pp.query("maxChi", maxChi);
  pp.query("stdev", stdev);

  Real fractionalInnerRadius = 0.2;
  pp.query("innerRadius",fractionalInnerRadius);
  Real innerRadius = fractionalInnerRadius*m_domainWidth;

  Real porosityTimescale = 1/m_parameters.darcy;
  pp.query("porosityTimescale", porosityTimescale);

  Real porosityEndTime = -1;
  pp.query("porosityEndTime", porosityEndTime);

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

      if (m_parameters.m_porosityFunction == m_parameters.m_porosityLinear)
      {

        Real boundaryPorosity = m_parameters.bcValPorosityLo[0];

        Real xDistFromMiddle = abs(x-m_domainWidth/2);
        Real yDistFromMiddle = abs(y-m_domainWidth/2);



        Real gradient = (1-boundaryPorosity)/(0.5*m_domainWidth - innerRadius);

        Real maxDistFromMiddle = max(xDistFromMiddle, yDistFromMiddle);

        //          Real porosity = maxChi + boundaryPorosity - 2*maxDistFromMiddle;
        Real porosity = 1.0;
        if (maxDistFromMiddle >= innerRadius)
        {
          porosity = 1-gradient*(maxDistFromMiddle-innerRadius);
        }

        Real maxPorosity = 1.0;
        porosity = min(porosity, maxPorosity);

        //initialDataPoiseuille(x, y, valNew,  m_porosity, -1);
        a_porosity[dit](iv) = porosity;
        //            (*m_scalarOld[m_porosity])[dit](iv) = porosity;

      }
      else if (m_parameters.m_porosityFunction == m_parameters.m_porosityGaussian)
      {

        Real boundaryPorosity = m_parameters.bcValPorosityLo[0];

        Real xd = abs(x-m_domainWidth/2);
        Real yd = abs(y-m_domainWidth/2);

        //Real maxDistFromMiddle = max(xDistFromMiddle, yDistFromMiddle);

        //Real porosity = 1 + boundaryPorosity - 2*maxDistFromMiddle;

        Real porosity = boundaryPorosity + maxChi*exp(-(xd*xd + yd*yd)/(stdev*m_domainWidth));

        Real maxPorosity = 1.0;
        porosity = min(porosity, maxPorosity);
        porosity = max(porosity, 0.0);

        //initialDataPoiseuille(x, y, valNew,  m_porosity, -1);
        a_porosity[dit](iv) = porosity;
        //            (*m_scalarOld[m_porosity])[dit](iv) = porosity;

      }
      else if (m_parameters.m_porosityFunction == m_parameters.m_porosityEdge)
      {
        Real boundaryPorosity = m_parameters.bcValPorosityLo[0];

        Real xd = abs(x-m_domainWidth/2);
        Real yd = abs(y-m_domainWidth/2);

        //Real maxDistFromMiddle = max(xDistFromMiddle, yDistFromMiddle);

        //Real porosity = 1 + boundaryPorosity - 2*maxDistFromMiddle;

        Real porosity = boundaryPorosity* (1-maxChi*exp(-(xd*xd + yd*yd)/(stdev*m_domainWidth)) );

        Real maxPorosity = 1.0;
        porosity = min(porosity, maxPorosity);
        porosity = max(porosity, 0.0);

        //initialDataPoiseuille(x, y, valNew,  m_porosity, -1);
        a_porosity[dit](iv) = porosity;
        //            (*m_scalarOld[m_porosity])[dit](iv) = porosity;
      }
      else if (m_parameters.m_porosityFunction == m_parameters.m_porosityConstant)
      {
        // Assume all boundary values are equal, and just pick one
        a_porosity[dit](iv) = m_parameters.bcValPorosityLo[0];
        //            (*m_scalarOld[m_porosity])[dit](iv) = m_parameters.bcValPorosityLo[0];
      }
      else if (m_parameters.m_porosityFunction == m_parameters.m_porosityTimeDependent)
      {
        Real time = m_time;
        if (porosityEndTime >= 0 && m_time >= porosityEndTime)
        {
          time = porosityEndTime;
        }
        Real scale = min((1-m_parameters.bcValPorosityHi[1]), time/porosityTimescale);
        scale = time/porosityTimescale;

        scale = min(scale, 4.0);
        //        Real scale = min(0.9,
        Real chi =  1-scale*pow(y/m_domainHeight, 2)*0.3*sin(M_PI*x/m_domainWidth*2);
        chi = min(chi, 1.0);
        chi = max(chi, m_lowerPorosityLimit);
        a_porosity[dit](iv) = chi;
      }

    }
  }
}



void AMRLevelMushyLayer::updateEnthalpyVariablesOld()
{
  ::updateEnthalpyVariables(*m_scalarOld[m_enthalpy], *m_scalarOld[m_bulkConcentration],
                            *m_scalarOld[m_temperature], *m_scalarOld[m_liquidConcentration], *m_scalarOld[m_solidConcentration],
                            *m_scalarOld[m_porosity],
                            *m_scalarOld[m_enthalpySolidus],*m_scalarOld[m_enthalpyLiquidus],*m_scalarOld[m_enthalpyEutectic],
                            m_parameters);
}



void AMRLevelMushyLayer::computeUDelU(LevelData<FArrayBox>& U_adv_src, const LevelData<FArrayBox>& advectionSourceTerm, Real a_oldTime, Real a_dt)
{
  // May want to use pre existing value for this, rather than re computing.

  ParmParse pp("main");

  LevelData<FArrayBox> UdelU_porosity(m_grids, SpaceDim, m_numGhostAdvection*IntVect::Unit);
  DataIterator dit = UdelU_porosity.dataIterator();

  if (m_doEulerPart)
  {
    int uDeluMethod = 0;
    pp.query("uDeluMethod", uDeluMethod);
    if (uDeluMethod == 0)
    {
      bool doVelFRupdates = true;
      pp.query("doAdvVelFRupdates", doVelFRupdates);
      predictVelocities(UdelU_porosity, m_advVel, advectionSourceTerm, a_oldTime, a_dt, doVelFRupdates);

    }
    else if (uDeluMethod == 1)
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

    //      setVelZero(UdelU_porosity, m_solidPorosity);


  }
  else
  {
    // Don't want this term
    for (dit.reset(); dit.ok(); ++dit)
    {
      UdelU_porosity[dit].setVal(0.0);
    }
  }

  int uDelU_grow = 1;
  pp.query("uDelU_grow", uDelU_grow);

  // by default make this tiny (so essentially turned off)
  Real uDelU_porosityLimit = 10*m_lowerPorosityLimit; //1e-15;
  pp.query("uDelu_porosity", uDelU_porosityLimit);
  setVelZero(UdelU_porosity, uDelU_porosityLimit, uDelU_grow); // grow by 1

  UdelU_porosity.copyTo(U_adv_src);

  for (dit.reset(); dit.ok(); ++dit)
  {
    (*m_vectorNew[m_UdelU])[dit].copy(UdelU_porosity[dit], 0, 0, SpaceDim);
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
  //	LevelData<FArrayBox> U_porosity_crsePtr;
  int srcGhost = src.ghostVect()[0];
  //	int nRefCrse = -1;

  ParmParse pp("ccSrc");
  bool pressureSrc = m_addSubtractGradP;
  bool darcySrc = explicitDarcyTerm;
  bool viscousSrc = false;
  bool buoyancySrc = true;
  bool advSrc = true;

  pp.query("pressure", pressureSrc);
  pp.query("darcy", darcySrc);
  pp.query("viscous", viscousSrc);
  pp.query("buoyancy", buoyancySrc);
  pp.query("advection", advSrc);

  if (m_time <= m_skipTrickySourceTerm)
  {
    pout() << "AMRLevelMushyLayer::computeUstarSrc - time = " << m_time << " < " << m_skipTrickySourceTerm;
    pout() << ", skipping advection src term";
    advSrc = false;
  }

  Real oldTime = m_time - m_dt;
  //  Real half_time = m_time - (m_dt / 2); // m_time is new time

  // Trying to work out freestream preservation stability
  //  half_time = oldTime;

  if (m_level > 0)
  {
    amrMLcrse = getCoarserLevel();

    crseOldVel.define(amrMLcrse->m_grids, SpaceDim,
                      srcGhost * IntVect::Unit);
    amrMLcrse->fillVectorField(crseOldVel, oldTime, m_fluidVel, true);

    //		nRefCrse = amrMLcrse->m_ref_ratio;
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
    //                 porosity_centred[dit],  (*m_vectorNew[m_bodyForce])[dit]);


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

    if (buoyancySrc)
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

    if (explicitDarcyTerm || darcySrc)
    {
      src[dit] -= darcy_src[dit];
    }

    // This correction is for when using VCAMRPoissonOp2
    //		src[dit].plus(U_correction[dit]);


  }

  src.exchange();

  src.copyTo(*m_vectorNew[m_viscousSolveSrc]);

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

  Real src_time_centering = 0.5;
  ParmParse pp("main");
  pp.query("vel_src_centering", src_time_centering);

  Real src_time = a_oldTime + src_time_centering * a_dt;

  computeUstarSrc(src, advectionSourceTerm, src_time, a_MACprojection, compute_uDelU);

  if (!m_isViscous)
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
    Vector<RefCountedPtr<LevelBackwardEuler> > UstarBE;
    Vector<RefCountedPtr<LevelTGA> > UstarTGA;
    defineUstarSolver(UstarBE, UstarTGA);

    // Do the solve for each component separately (as the solver is hardwired for one component only)
    LevelData<FArrayBox> UstarCrse, UoldCrse, *UstarCrseComp, *UoldCrseComp, Uold;
    LevelData<FluxBox> diffusiveFlux(m_grids, SpaceDim);
    LevelFluxRegister *fineFluxRegPtr=NULL, *crseFluxRegPtr=NULL;

    if (m_vectorFluxRegisters[m_fluidVel])
    {
      fineFluxRegPtr = &(*m_vectorFluxRegisters[m_fluidVel]);
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

      crseFluxRegPtr = &(*(amrMLcrse->m_vectorFluxRegisters[m_fluidVel]));
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
    //todo - surely can solve both components at the same time? Wouldn't this be quicker with multigrid?
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

//      if (s_verbosity > 2)
//      {
//        pout() << "Viscous solve (level " << m_level << ") on component "
//            << comp << endl;
//      }

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

      //			s_uStarAMRMG[comp]->solve(UstarVectComp, rhsVectComp, maxLevel,
      //					maxLevel, false); // don't initialize to zero
      //			exitStatus += s_uStarAMRMG[comp]->m_exitStatus;

      if (!doFRupdates)
      {
        crseFluxRegPtr = NULL;
        fineFluxRegPtr = NULL;
      }

      int exitStatus = -1;
      Real resid = 0;

      if (m_timeIntegrationOrder == 1)
      {
        UstarBE[comp]->updateSoln(UstarComp, UoldComp, compSrc,
                                  &(*fineFluxRegPtr), &(*crseFluxRegPtr),
                                  UoldCrseComp, UstarCrseComp,
                                  old_time,  old_crseTime, new_crseTime,
                                  a_dt, m_level, false, comp); // False - don't zero phi


        exitStatus = UstarBE[comp]->exitStatus();
        resid = UstarBE[comp]->finalResidual();

      }
      else
      {
        UstarTGA[comp]->updateSoln(UstarComp, UoldComp, compSrc,
                                   &(*fineFluxRegPtr), &(*crseFluxRegPtr),
                                   UoldCrseComp, UstarCrseComp,
                                   old_time, old_crseTime, new_crseTime,
                                   a_dt, m_level, false, comp);  // False - don't zero phi

        exitStatus = UstarTGA[comp]->exitStatus();
        resid = UstarTGA[comp]->finalResidual();
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



  } // end if viscous

}


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

  // TODO - fill this with -dp/dz - Theta_l
  for (DataIterator dit = velocityBCVals->dataIterator(); dit.ok(); ++dit)
  {
    FluxBox& bcVel = (*velocityBCVals)[dit];

    FArrayBox& bc_u = bcVel[0];
    FArrayBox& bc_v = bcVel[1];

    FArrayBox& Sl = Theta_l_face[dit][1];
    FArrayBox& gradP_z = gradP[dit][1];

    bc_u.setVal(0.0);

    // vertical velocity boundary cells:
//    Box v_boundary = ::adjCellLo(m_problem_domain, 1, 1);
//    v_boundary.shift(1, 1);
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

    // Testing
//    bc_v.setVal(-20);

  }

  EdgeVelBCHolder edgeVelBC(m_physBCPtr->edgeVelFuncBC(m_viscousBCs, velocityBCVals));

  VelBCHolder velBC(m_physBCPtr->uStarFuncBC(m_viscousBCs));
  VelBCHolder velBCExtrap(m_physBCPtr->velExtrapBC());

  calculatePermeability();

  //  bool quadInterp = true;

  // if a coarser level exists, will need coarse-level data for proj
  if (m_level > 0)
  {
    const DisjointBoxLayout& crseGrids = amrMLcrse->m_grids;
    //    crseVelNewPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim);
    //    crseVelOldPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim);
    crsePressurePtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(crseGrids, 1));
    // coarse velocity BC data is interpolated in time

    // If we're not doing +/- grad(P) stuff then just pass in the velocity before it was projected
    // else the divergence at coarse-fine boundaries will be crazy
    //    amrMLcrse->fillVectorField(*crseVelNewPtr, new_time, m_Ustar, true);
    //    amrMLcrse->fillVectorField(*crseVelOldPtr, old_time, m_Ustar, true);

    amrMLcrse->fillScalars(*crsePressurePtr, time, m_pressure, true);

    if (m_scaleP_MAC)
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

  if(m_scaleP_MAC)
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

  if (m_isViscous)
  {
    if (s_verbosity >= 5)
    {
      pout() << "AMRLevelMushyLayer::calculateTimeIndAdvectionVel - with viscosity (level " << m_level << ")"    << endl;
    }

    // Average U^* to faces for projection
    CellToEdge(*m_vectorNew[m_Ustar], a_advVel);
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
    bool useIncrementalPressure = (getMaxLevel() == 0);
    ParmParse pp("projection");
    pp.query("useIncrementalPressure", useIncrementalPressure);

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
        Real basic_phiScale = 1;

        pp.get("phiScale", basic_phiScale);

        Real phiScale = m_projection.getScale(basic_phiScale, m_dt);

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
  Real smoothingCoeff = 0.0; //0.01;
  ParmParse ppProj("projection");
  ppProj.query("pre_smoothing", smoothingCoeff);

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

        velDir(iv) = velDir(iv) - smoothingCoeff*lapU;
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
    EdgeToCell(a_advVel, *m_vectorNew[m_advUstar]);

    a_advVel.exchange();
    int exitStatus = m_projection.levelMacProject(a_advVel, m_dt, crsePressurePtr, pressureScalePtr,
                                 crsePressureScalePtr, pressureScaleEdgePtr, crsePressureScaleEdgePtr,
                                 alreadyHasPressure, correctScale);

    correctScale = correctScale/10;
//    correctScale = 0.0;

    a_advVel.exchange();
    Divergence::levelDivergenceMAC(*m_scalarNew[m_divUadv], a_advVel, m_dx);
    maxDivU = ::computeNorm(*m_scalarNew[m_divUadv], NULL, 1, m_dx, Interval(0,0), 0);
    Real minPressure = ::computeMin(m_projection.phi(), NULL, 1, Interval(0,0));
    pout() << "  MAC Projection (level "<< m_level << "), exit status = " << exitStatus
        << ", max(div u) = " << maxDivU << ", min pressure = " << minPressure << endl;



    proj_i++;
  }


  fillAdvVel(time, a_advVel);
  a_advVel.exchange();

  m_projection.getPhi(*m_scalarNew[m_pressure]);

  // Test - overwrite with analytic advection velocity
  bool enforceAnalyticVel = (m_parameters.physicalProblem == MushyLayerParams::m_soluteFluxTest  );

  enforceAnalyticVel = false;

  ParmParse pp("main");
  pp.query("analyticVel", enforceAnalyticVel);

  if (enforceAnalyticVel)
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
  Divergence::levelDivergenceMAC(*m_scalarNew[m_divUadv], a_advVel, m_dx);

  // Need a m_fluidVel to be populated to some advection things the right way around
  EdgeToCell(m_advVel, *m_vectorNew[m_fluidVel]);
  this->fillVectorField(*m_vectorNew[m_fluidVel], time, m_fluidVel);

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

int AMRLevelMushyLayer::getMaxLevel()
{
  int maxLevel = -1;
  ParmParse pp("main");
  pp.get("max_level", maxLevel);

  return maxLevel;
}

void AMRLevelMushyLayer::fillUnprojectedDarcyVelocity(LevelData<FluxBox>& a_advVel, Real time)
{
  IntVect ghost = a_advVel.ghostVect();
  DataIterator dit = a_advVel.dataIterator();



  LevelData<FluxBox> permeability_face(m_grids, 1, ghost);
  LevelData<FluxBox> T_face(m_grids, 1, ghost);
  LevelData<FluxBox> C_face(m_grids, 1, ghost);
  fillScalarFace(permeability_face, time, m_permeability,        true, false);
  fillScalarFace(T_face,            time, m_temperature,         true, false);
  fillScalarFace(C_face,            time, m_liquidConcentration, true, false);

  for (dit.reset(); dit.ok(); ++dit)
  {
    a_advVel[dit].setVal(0.0);

    FArrayBox& fabVelz = a_advVel[dit][SpaceDim-1]; // for debugging

    //      fabVelz.plus(T_face[dit][SpaceDim-1],  m_parameters.rayleighTemp);
    //      fabVelz.plus(C_face[dit][SpaceDim-1], -m_parameters.rayleighComposition);
    fabVelz.plus(T_face[dit][SpaceDim-1],  m_parameters.m_buoyancyTCoeff);
    fabVelz.plus(C_face[dit][SpaceDim-1], -m_parameters.m_buoyancySCoeff);

    for (int idir=0; idir<SpaceDim; idir++)
    {
      a_advVel[dit][idir].mult(permeability_face[dit][idir]);
//      a_advVel[dit][idir].divide(m_parameters.m_darcyCoeff);
    }

  }
}

void AMRLevelMushyLayer::fillAnalyticVel(FArrayBox& velDir, int dir, int comp, bool project)
{
  const Box& b = velDir.box();

  ParmParse ppMain("main");

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

    Real freqx = (m_problem_domain.isPeriodic(0)) ? 2*M_PI/m_domainWidth : M_PI/m_domainWidth ;
    //    Real freqy = (m_problem_domain.isPeriodic(1)) ? 2*M_PI/m_domainWidth : M_PI/m_domainWidth ;

    Real Uscale;
    if (m_parameters.nonDimVel > 0)
    {
      Uscale = m_parameters.nonDimVel*100;
    }
    else if(m_parameters.darcy > 0)
    {
      Uscale = max(m_parameters.rayleighTemp, m_parameters.rayleighComposition)*m_parameters.darcy;
    }
    else
    {
      Uscale = sqrt(max(m_parameters.rayleighTemp, m_parameters.rayleighComposition));
    }

    Real magU = Uscale; //*0.15;

    int analyticVelType = m_parameters.physicalProblem;
    ppMain.query("analyticVelType", analyticVelType);

    if (analyticVelType == -1)
    {
      //      Real shift = 0;

      // This velocity should have zero divergence
      if (dir == 0)
      {
        //        velDir(iv, comp) = magU*sin(freqx*(x+shift))*sin(freqx*(x+shift))*sin(freqx*(y+shift))*cos(freqx*(y+shift));
        velDir(iv, comp) = 0;
      }
      else
      {
        velDir(iv, comp) = magU*sin(freqx*x)*(m_domainHeight-y);
      }

    }
    else if(analyticVelType == MushyLayerParams::m_soluteFluxTest)
    {

      if (dir == 0)
      {
        //              y = y + m_dx/2;
        velDir(iv, comp) = magU*(cos(freqx*x)+0.5*freqx*x)*exp(-freqx*y);
      }
      else
      {
        //              x = x + m_dx/2;
        velDir(iv, comp) = magU*(sin(freqx*x))*sin(freqx*y);

        //              fab(iv) = -magU*cos(freq*x)*sin(freq*y);

        //velDir(iv) = sin(freq*x);
      }

    }
    else if (analyticVelType == MushyLayerParams::m_mushyLayer)
    {
      if (project)
      {
        // This version is divergent, and needs to be projected
        if (dir == 0)
        {
          //              y = y + m_dx/2;
          velDir(iv, comp) = -magU*sin(freqx*x)*sin(freqx*x)*cos(freqx*y);
        }
        else
        {
          //              x = x + m_dx/2;
          velDir(iv, comp) = magU*cos(freqx*x)*sin(freqx*y)*sin(freqx*y);
        }
      }
      else
      {
        Real shift = 0.5*m_dx;
        shift = 0;
        // This velocity should have zero divergence
        if (dir == 0)
        {
          //          velDir(iv, comp) = magU*sin(freqx*x)*cos(freqx*y);
          //          velDir(iv, comp) = magU*sin(freqx*x)*cos(freqx*y);

          velDir(iv, comp) = magU*sin(freqx*(x+shift))*sin(freqx*(x+shift))*sin(freqx*(y+shift))*cos(freqx*(y+shift));
        }
        else
        {

          velDir(iv, comp) = -magU*sin(freqx*(x+shift))*cos(freqx*(x+shift))*sin(freqx*(y+shift))*sin(freqx*(y+shift));
          //          velDir(iv, comp) = -magU*cos(freqx*x)*sin(freqx*y);
        }
      }
    }
    else
    {
      if (dir == 0)
      {
        velDir(iv, comp) = -magU*sin(freqx*x)*cos(freqx*y);
      }
      else
      {
        velDir(iv, comp) = magU*cos(freqx*x)*sin(freqx*y);

      }
    }

  }
}


void AMRLevelMushyLayer::fillAnalyticVel(LevelData<FArrayBox>& a_advVel)
{
  DataIterator dit = m_grids.dataIterator();

  ParmParse ppMain("main");
  bool project = false;
  ppMain.query("analyticVelNonDivergent", project);

  for (dit.reset(); dit.ok(); ++dit)
  {
    FArrayBox& vel = a_advVel[dit];
    for (int dir=0; dir < SpaceDim; dir++)
    {
      fillAnalyticVel(vel, dir, dir, project);
    }
  }
}

void AMRLevelMushyLayer::fillAnalyticVel(LevelData<FluxBox>& a_advVel)
{
  DataIterator dit = m_grids.dataIterator();

  ParmParse ppMain("main");
  bool project = false;
  ppMain.query("analyticVelProject", project);

  for (dit.reset(); dit.ok(); ++dit)
  {
    for (int dir = 0; dir < SpaceDim; dir++)
    {
      FArrayBox& velDir = a_advVel[dit][dir];
      fillAnalyticVel(velDir, dir, 0, project);
    }
  }


}

void AMRLevelMushyLayer::fillAdvVel(Real time, LevelData<FluxBox>& a_advVel)
{
  // Fill interior ghost cells
  if (m_level > 0)
  {
    int interpRadius = a_advVel.ghostVect()[0];

    AMRLevelMushyLayer* amrMLcrse = getCoarserLevel();

    // Refinement ratio between this and coarser level
    int crseRefRatio = amrMLcrse->m_ref_ratio;

    ProblemDomain crseDomain = amrMLcrse->m_problem_domain;
    DisjointBoxLayout crseGrids(amrMLcrse->m_grids);

    //    bool secondOrderCorners = (CFinterpOrder_advection == 2);
    PiecewiseLinearFillPatchFace patcher(m_grids, crseGrids, 1,
                                         crseDomain, crseRefRatio,
                                         interpRadius);
    //                                         secondOrderCorners);

    Real timeInterpCoeff = 0;

    //    Real crseOldTime = amrMLcrse->time() - amrMLcrse->dt();
    //    Real crseNewTime = amrMLcrse->time();
    //
    //    if (abs(time - crseOldTime) < TIME_EPS)
    //    {
    //      timeInterpCoeff = 0;
    //    }
    //    else if (abs(time - crseNewTime) < TIME_EPS)
    //    {
    //      timeInterpCoeff = 1;
    //    }
    //    else
    //    {
    //      // do linear interpolation in time
    //      timeInterpCoeff = (time - crseOldTime)/(crseNewTime - crseOldTime);
    //    }

    // As advection velocities are stored at the half time step, we can't
    // necessarilly interpolate at exactly the right point in time, so accept
    // that we're making an O(dt) error here.
    patcher.fillInterp(a_advVel, amrMLcrse->m_advVel, amrMLcrse->m_advVel,
                       timeInterpCoeff, 0, 0, 1);


  }

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

  ParmParse ppMain("main");

  const DisjointBoxLayout& levelGrids = a_advVel.getBoxes();
  // for tracing, will need to fill in boundary values for
  // grown copy of velocity

  IntVect advectionGhostVect = m_numGhostAdvection * IntVect::Unit;
  //  LevelData<FArrayBox> viscousSource(levelGrids, SpaceDim, ghostVect);

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

    if ( m_advectionMethod == m_noPorosity )
    {
      // do nothing
    }
    else if ( m_advectionMethod == m_porosityOutsideAdvection
        || m_advectionMethod == m_porosityInAdvection)
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

  if (m_advectionMethod == m_porosityInAdvection)
  {
    fillVectorField(UtoAdvect_old, old_time, m_U_porosity, true);
  }
  else if (m_advectionMethod == m_porosityOutsideAdvection ||  m_advectionMethod == m_noPorosity)
  {
    fillVectorField(UtoAdvect_old, old_time, m_fluidVel, true);
  }

  Real advVelChiLimit = min(pow(10,5)*m_lowerPorosityLimit, pow(10,-10)) ; //was 1e-10

  ppMain.query("advPorosityLimit", advVelChiLimit);
  setVelZero(UtoAdvect_old, advVelChiLimit);
  setVelZero(advectionVelocity, advVelChiLimit);

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

    bool legacyCompute = false;

    ppMain.query("legacy_predict_vel", legacyCompute);

    // Get the true pressure correction (scaled with chi if appropriate)
    if (m_scaleP_MAC)
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
                               legacyCompute);

    U_chi_advected.exchange();

    for (DataIterator dit = U_chi_advected.dataIterator(); dit.ok(); ++dit)
    {
      for (int dir=0; dir<SpaceDim; dir++)
      {
//        U_chi_advected[dit].copy();
        ::EdgeToCell(U_chi_advected[dit], dir, (*m_vectorNew[m_U_porosity])[dit], 0, dir);
      }
    }
//    EdgeToCell(, *m_vectorNew[m_U_porosity]);

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


  int conservativeForm = false;
  ParmParse pp("main");
  pp.query("uDelU_conservativeForm", conservativeForm);

  if (conservativeForm)
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
          if (s_reflux_normal_momentum)
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
          if (s_reflux_normal_momentum)
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
  //	CellToEdge(a_old_vel, a_advVel);
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

    EdgeToCell(a_advVel, *m_vectorNew[m_advectionVel]);
    Real maxAdvU = ::computeNorm(*m_vectorNew[m_advectionVel], NULL, 1, m_dx, Interval(0,0), 0);

    if (s_verbosity >= 1 && maxAdvU > 1e10)
    {
      pout() << "WARNING - max advection velocity (post tracing, pre-projection) = " << maxAdvU << endl;
    }

  } // end loop over grids

  //  int temp = 0;
}

//Copied from AMRINS/AdvectUtil
void AMRLevelMushyLayer::computeGradP(LevelData<FArrayBox>& gradP,
                                      Real a_time,
                                      bool a_macProjection)
{
  DataIterator dit(m_grids);

  if (m_enforceGradP)
  {
    for (dit.reset(); dit.ok(); ++dit)
    {
      BoxIterator bit(gradP[dit].box());
      for (bit.reset(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        ::getLocation(iv, loc, m_dx);

        if (m_parameters.physicalProblem == MushyLayerParams::m_poiseuilleFlow)
        {
          gradP[dit](iv, 0) = 0;
          //					gradP[dit](iv, 1) = -stokesDarcyForcing(a_time);
          gradP[dit](iv, 1) = 0;
        }
        else
        {
          gradP[dit](iv, 0) = 0;
          gradP[dit](iv, 1) = 0;
        }
      }
    }

    return;
  }

  //note we use extrap BCs to try and make things a little smoother
  m_projection.gradPiBCs(gradP, true, a_macProjection);

  // CF Boundary conditions
  if (m_level > 0)
  {

    AMRLevelMushyLayer* amrMushyLayerCoarserPtr = getCoarserLevel();
    const DisjointBoxLayout& crseGrids = amrMushyLayerCoarserPtr->m_grids;

    LevelData<FArrayBox> crseGradP(crseGrids, SpaceDim, gradP.ghostVect());

    amrMushyLayerCoarserPtr->m_projection.gradPi(crseGradP);

    const ProblemDomain& crseDomain = amrMushyLayerCoarserPtr->problemDomain();
    int nRefCrse = amrMushyLayerCoarserPtr->refRatio();

    // This doesn't matter - new and old crse grad P are the same
    Real crse_time_interp_coeff = 1.0;
    int velGrow = gradP.ghostVect()[0];

    //    bool secondOrderCorners = (CFinterpOrder_advection==2);

    PiecewiseLinearFillPatch filpatcher(m_grids, crseGrids, SpaceDim,
                                        crseDomain, nRefCrse,
                                        velGrow,
                                        false); // not constant interp
    //                                        secondOrderCorners);

    filpatcher.fillInterp(gradP, crseGradP, crseGradP,
                          crse_time_interp_coeff, 0, 0, SpaceDim);
  }

}

void AMRLevelMushyLayer::upwind(LevelData<FluxBox>& a_edgeScal,
                                LevelData<FArrayBox>& a_old_scal,
                                LevelData<FluxBox>& a_adv_vel,
                                LevelData<FluxBox>& a_inflowOutflowVel,
                                LevelData<FArrayBox>& a_old_vel,
                                LevelData<FArrayBox>& a_diffusiveSrc,
                                PatchGodunov& a_patchGodScalar,
                                Real a_old_time, Real a_dt)
{
  // Get AdvectionPhysics object within the PatchGodunov object
  AdvectionPhysics* advectionPhysics =
      dynamic_cast<AdvectionPhysics*>(a_patchGodScalar.getGodunovPhysicsPtr());
  if (advectionPhysics == NULL)
  {
    MayDay::Error("AMRLevelMushyLayer::upwind - unable to upcast GodunovPhysics to AdvectionPhysics");
  }

  // set up patchGodunov for this problem
  a_patchGodScalar.setCurrentTime(a_old_time);

  // also need to build a grown advection velocity
  IntVect advectVelGrow(3*IntVect::Unit);
  const DisjointBoxLayout& levelGrids = a_old_vel.getBoxes();
  LevelData<FluxBox> grownAdvVel(levelGrids, 1, advectVelGrow);

  // now overwrite with advection velocities wherever possible
  a_adv_vel.copyTo(a_adv_vel.interval(), grownAdvVel,
                   grownAdvVel.interval());

  // now trace scalars to edges at time n+1/2
  DataIterator dit = levelGrids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    FluxBox& thisEdgeScal = a_edgeScal[dit()];
    FluxBox& thisAdvVel = grownAdvVel[dit()];
    FArrayBox& thisCellVel = a_old_vel[dit()];
    FArrayBox& thisOldScal = a_old_scal[dit()];
    FArrayBox& thisSrc = a_diffusiveSrc[dit()];

    FluxBox& thisInflowOutflowVel = a_inflowOutflowVel[dit()];

    // Init to zero
    thisEdgeScal.setVal(0.0);

    Box curBox(thisSrc.box());
    curBox.grow(-1);

    a_patchGodScalar.setCurrentBox(curBox);
    advectionPhysics->setCellVelPtr(&thisCellVel); // set cell-centered velocity field
    advectionPhysics->setAdvVelPtr(&thisAdvVel); // set advection velocity field
    advectionPhysics->setInflowOutflowVelPtr(&thisInflowOutflowVel);

    // compute face-centered, predicted scalars
    a_patchGodScalar.computeWHalf(thisEdgeScal, thisOldScal,
                                  thisSrc, a_dt, curBox); // was curBox


  }
}

void AMRLevelMushyLayer::
computeScalarAdvectiveFlux(LevelData<FluxBox>& a_edgeScal,
                           LevelData<FArrayBox>& a_old_scal,
                           LevelData<FluxBox>& a_adv_vel,
                           LevelData<FluxBox>& a_inflowOutflowVel,
                           LevelData<FArrayBox>& a_old_vel,
                           LevelData<FArrayBox>& a_diffusiveSrc,
                           PatchGodunov& a_patchGod,
                           Real a_old_time, Real a_dt)
{
  CH_TIME("AMRLevelMushyLayer::computeScalarAdvectiveFlux");

  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::computeScalarAdvectiveFlux, level " << m_level  << endl;
  }

  int numScal = a_edgeScal.nComp();

  // Predict half time face centered scalar components
  upwind(a_edgeScal, a_old_scal, a_adv_vel, a_inflowOutflowVel, a_old_vel, a_diffusiveSrc, a_patchGod, a_old_time, a_dt);

  DataIterator dit = a_old_vel.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    FluxBox& thisEdgeScal = a_edgeScal[dit()];

    for (int dir=0; dir<SpaceDim; dir++)
    {
      // multiply by edge velocity to get flux
      // do componentwise
      for (int comp=0; comp<numScal; comp++)
      {
        thisEdgeScal[dir].mult(a_adv_vel[dit][dir],0,comp,1);
        thisEdgeScal[dir].mult(m_parameters.m_advectionCoeff);
      }

    } // end loop over tracing directions

  } // end loop over grids for tracing



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
    if (m_advectionMethod == m_porosityOutsideAdvection)
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

    bool alterPredVel = true;
    ParmParse ppMain("main");
    ppMain.query("alter_predicted_vels", alterPredVel);

    if (alterPredVel)
    {
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

    } // end if altering our predicted velocities


  } // end loop over grids
}

void AMRLevelMushyLayer::calculatePermeability()
{
  if (s_verbosity >= 5)
  {
    pout() << "  AMRLevelMushyLayer::calculatePermeability" << endl;
  }
  LevelData<FArrayBox> porosityNew(m_grids, 1, m_numGhost*IntVect::Unit);
  LevelData<FArrayBox> porosityOld(m_grids, 1, m_numGhost*IntVect::Unit);

  fillScalars(porosityOld, m_time-m_dt, m_porosity, true);
  fillScalars(porosityNew, m_time, m_porosity, true);

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    FArrayBox solidFraction((*m_scalarNew[m_permeability])[dit].box(), 1);

    //First do new timestep
    solidFraction.setVal(1.0);
    solidFraction -= porosityNew[dit];
    ::calculatePermeability((*m_scalarNew[m_permeability])[dit],
                            solidFraction, m_parameters, m_dx);

    // Now do old timestep
    solidFraction.setVal(1.0);
    solidFraction -= porosityOld[dit];
    ::calculatePermeability((*m_scalarOld[m_permeability])[dit],
                            solidFraction, m_parameters, m_dx);
  }

  calculateCoarseFineBoundaries(m_permeability);
}

void AMRLevelMushyLayer::calculateCoarseFineBoundaries(int a_var, bool a_isVector)

{
  if (m_level == 0)
  {
    return;
  }

  AMRLevelMushyLayer *amrMLcrse = getCoarserLevel();

  if (a_isVector)
  {
    m_quadCFInterpVector.coarseFineInterp(*m_vectorNew[a_var],
                                          *amrMLcrse->m_vectorNew[a_var]);
    m_quadCFInterpVector.coarseFineInterp(*m_vectorOld[a_var],
                                          *amrMLcrse->m_vectorOld[a_var]);
  } else {
    m_quadCFInterpScalar.coarseFineInterp(*m_scalarNew[a_var],
                                          *amrMLcrse->m_scalarNew[a_var]);
    m_quadCFInterpScalar.coarseFineInterp(*m_scalarOld[a_var],
                                          *amrMLcrse->m_scalarOld[a_var]);
  }

}

void
AMRLevelMushyLayer::computeScalDiffusion(const int a_var,
                                         LevelData<FArrayBox>& a_lapScal,
                                         Real a_time)
{
  CH_TIME("AMRLevelMushyLayer::computeScalDiffusion");
  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::computeScalDiffusion for var " << a_var  << endl;
  }

  // Compute diffusion for both fields
  LevelData<FArrayBox> diffusiveSrc(m_grids, 2, IntVect::Zero);
  computeScalDiffusion(diffusiveSrc, a_time);

  // Extract the field we want
  //  for (dit.reset(); dit.ok(); ++dit)
  //    {
  //     // a_var = 0 (H, T) or 1 (S, S_l)
  //      a_lapScal[dit].copy(diffusiveSrc[dit], a_var, 0, 1);
  //
  //    }

  diffusiveSrc.copyTo(Interval(a_var, a_var), a_lapScal, Interval(0,0));

}

void
AMRLevelMushyLayer::computeScalDiffusion(LevelData<FArrayBox>& diffusiveSrc,
                                         Real a_time)
{

  CH_TIME("AMRLevelMushyLayer::computeScalDiffusion");

  bool isHomogeneous = false;

  int numComp = 2;
  CH_assert(diffusiveSrc.nComp() == numComp);

  RefCountedPtr<AMRNonLinearMultiCompOp> amrpop = RefCountedPtr<AMRNonLinearMultiCompOp>(
      (AMRNonLinearMultiCompOp*) HCOpFact->AMRnewOp(m_problem_domain));

  LevelData<FArrayBox> HC(m_grids, 2, IntVect::Unit);
  LevelData<FArrayBox> *crseHC = NULL;

  fillHC(HC, a_time, true, true);

  // Get coarse BC if there is one
  if (AMRLevelMushyLayer* mlCrse = getCoarserLevel())
  {
    crseHC = new LevelData<FArrayBox>(mlCrse->m_grids, numComp, IntVect::Unit);
    mlCrse->fillHC(*crseHC, a_time, true, true);
  }

  amrpop->setAlphaAndBeta(0, 1);

  // This just calls applyOpI if crseHC = NULL, else does CF interpolation
  amrpop->applyOpMg(diffusiveSrc, HC, crseHC, isHomogeneous);

  int order = 1;
  BCHolder bc = m_physBCPtr->extrapolationFuncBC(order);

  DataIterator dit = diffusiveSrc.dataIterator();
  const DisjointBoxLayout& grids = diffusiveSrc.getBoxes();
  for (dit.reset(); dit.ok(); ++dit)
  {
    bc(diffusiveSrc[dit],
       grids[dit],
       m_problem_domain, m_dx,
       false); // not homogeneous

    //    int temp=0;
  }


  // finally, do exchange
  diffusiveSrc.exchange(diffusiveSrc.interval());

  // Think we also need a corner copier (if we have ghost cells)
  if (diffusiveSrc.ghostVect()[0] > 0)
  {
    CornerCopier cornerCopy(grids, grids, m_problem_domain, diffusiveSrc.ghostVect(), true);
    diffusiveSrc.exchange(diffusiveSrc.interval(), cornerCopy);
  }

  // Clean up
  if (crseHC != NULL)
  {
    delete crseHC;
    crseHC = NULL;
  }

}

//This function taken from the Navier Stokes code
void AMRLevelMushyLayer::computeLapVel(LevelData<FArrayBox>& a_lapVel,
                                       LevelData<FArrayBox>& a_vel, const LevelData<FArrayBox>* a_crseVelPtr)
{

  // These are set in applyOpI anyway, but do it again here so we can see what happens
  VelBCHolder velBC(m_physBCPtr->uStarFuncBC(m_viscousBCs));
  velBC.applyBCs(a_vel, a_vel.disjointBoxLayout(), m_problem_domain, m_dx, false);

  // Copy velocity to new data holder before taking laplacian to avoid changing
  // original velocity
  LevelData<FArrayBox> a_vel_copy(m_grids, SpaceDim, 2*IntVect::Unit); // 2 ghost vectors for higher order laplace
  DataIterator dit(m_grids);
  for (dit.begin(); dit.ok(); ++dit)
  {
    //		a_vel_copy[dit].copy(a_vel[dit]);
    Box b = m_grids[dit];
    b&=m_problem_domain.domainBox();
//    a_vel_copy[dit].setVal(0.0);
//    a_vel_copy[dit] += a_vel[dit];
//    a_vel_copy[dit].plus(a_vel[dit], b, 0, 0, SpaceDim);
//    		int temp=0;
  }

  a_vel.copyTo(Interval(0, SpaceDim-1), a_vel_copy, Interval(0, SpaceDim-1));


  // Need 2nd order BCs for laplace operator
  if (a_crseVelPtr)
  {

    m_quadCFInterpVector.coarseFineInterp(a_vel_copy, *a_crseVelPtr);
    //		extrapBC.applyBCs(a_vel_copy, a_vel_copy.disjointBoxLayout(), m_problem_domain, m_dx, false);
  }

  // Apply domain BCs again?
  velBC.applyBCs(a_vel_copy, a_vel_copy.disjointBoxLayout(), m_problem_domain, m_dx, false);


  a_vel_copy.exchange(Interval(0, SpaceDim-1));

  ParmParse pp("computeLapVel");

  bool HOLapVel = false;
  pp.query("HO", HOLapVel);
  //todo - make HO work with AMR?


  for (int dir = 0; dir < SpaceDim; dir++)
  {
    // Fill outer ghost cells
    // Shouldn't need to do this any more - velBC should fill all ghost cells
//    VelBCHolder extrapBC = m_physBCPtr->velExtrapBC(Interval(1,1));
//    extrapBC.applyBCs(a_vel_copy, a_vel_copy.disjointBoxLayout(), m_problem_domain, m_dx, false);

    // todo: Try and replace this with the Darcy Brinkman op?
    if (HOLapVel)
    {

      MayDay::Error("Higher Order laplacian of velocity not implemented.");
//      m_viscousOp[dir]->applyOpI4(a_lapVel, a_vel_copy, false);

    }
    else
    {
      m_viscousOp[dir]->applyOpI(a_lapVel, a_vel_copy, false);
    }
    a_lapVel.exchange(a_lapVel.interval());

    // Do some smoothing

    int nsmooth = 0;
    Real a = 0.0;
    pp.query("num_passes", nsmooth);
    pp.query("scale", a);
    Box domBox = m_problem_domain.domainBox();
    for (int dir =0; dir <SpaceDim; dir ++)
    {
      if (!m_problem_domain.isPeriodic(dir))
      {
        domBox.grow(dir, -1);
      }
    }
//    domBox.grow(-1); // need to boundaries differently (one sided diffs)
    for (int i=0; i < nsmooth; i++)
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

          thisLap(iv, dir) = (1-a)*thisLap(iv, dir) + (a/4)*(thisLap(iv+BASISV(0), dir) +
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
    int order = 0;
    ParmParse ppMain("main");
    ppMain.query("lap_vel_bc_order", order);
    if (order >= 0)
    {
      BCHolder bc = m_physBCPtr->extrapolationFuncBC(order);
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

void AMRLevelMushyLayer::horizontallyAverage(LevelData<FArrayBox>& a_averaged, LevelData<FluxBox>& a_phi)
{
  Vector<Real> discardVector;
  horizontallyAverage(a_averaged, a_phi, discardVector);
}

void AMRLevelMushyLayer::horizontallyAverage(Vector<Real>& globalAveraged, LevelData<FluxBox>& a_phi  )
{
  LevelData<FArrayBox> discardLD(a_phi.disjointBoxLayout(), 1, IntVect::Zero);
  horizontallyAverage(discardLD, a_phi, globalAveraged);

}

void AMRLevelMushyLayer::horizontallyAverage(LevelData<FArrayBox>& a_averaged, LevelData<FluxBox>& a_phi,
                                             Vector<Real>& globalAveraged)
{
  //  a_phi.exchange();
  Box domBox = m_problem_domain.domainBox();

  IntVect smallEnd =domBox.smallEnd();
  int y_init = smallEnd[1];

  //  int nGhost = a_phi.ghostVect()[0];
  int idir = SpaceDim-1;

  // +2 because for N cells we have N+1 faces
  Real length = m_numCells[SpaceDim-1]+2;
  globalAveraged.resize(length, 0.0);

  for (int comp = 0; comp < a_phi.nComp(); comp++)
  {
    Vector<Real> averaged(length, 0.0);

    // Ensure we initialise this to 0
    for (int y_i = 0; y_i < globalAveraged.size(); y_i++)
    {
      globalAveraged[y_i] = 0;
    }

    Vector<Vector<Real> > allAveraged(numProc());
    for (int i = 0; i<numProc(); i++)
    {
      allAveraged[i].resize(length);
      allAveraged[i].assign(0.0); //= new Vector<Real>(length+2, 0.0);
    }

    DisjointBoxLayout dbl = a_phi.disjointBoxLayout();




    for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {

      FArrayBox& fluxDir = a_phi[dit][idir];
      Box box = fluxDir.box();

      // To stop double counting of interior faces
      // think this was actually causing us to not calculated fluxes at the top of boxes
      if (box.hiVect()[1] < domBox.hiVect()[1])
      {
        box.growDir(idir, Side::Hi, -1);
      }



      for (BoxIterator bit(box); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        int y_i = iv[SpaceDim-1];
        IntVect ivUp = iv + BASISV(idir);

        averaged[y_i-y_init] += fluxDir(iv, comp)*m_dx/m_domainWidth;
      }
    }

    //  pout() << endl;

    // Broadcast/gather to compute averages over whole domain
    int srcProc = 0;

    gather(allAveraged, averaged, srcProc);
    if (procID() == srcProc)
    {
      for (int ivec = 0; ivec<numProc(); ivec++)
      {
        for (int y_i = 0; y_i < globalAveraged.size(); y_i++)
        {
          globalAveraged[y_i] += (allAveraged[ivec])[y_i];
        }
      }
    }

    broadcast(globalAveraged, srcProc);

    //  pout() << averaged << endl;
    //    if (s_verbosity > 10)
    //    {
    //      pout() << "Processor ID: " << procID() << endl;
    //      pout() << "solute fluxes: " << globalAveraged << endl;
    //    }

    for (DataIterator dit = a_averaged.dataIterator(); dit.ok(); ++dit)
    {
      //      Box box = m_grids[dit];

      // To get the top face
      //      box.growDir(1, Side::Hi, 1);
      Box b = a_averaged[dit].box();

      for (BoxIterator bit(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        int y_i = iv[SpaceDim-1];
        y_i = y_i - y_init;
        a_averaged[dit](iv, comp) = globalAveraged[y_i];
      }
    }

  } // end loop over components

  a_averaged.exchange();

}

void AMRLevelMushyLayer::horizontallyAverage(Vector<Real>& averageVector, LevelData<FArrayBox>& a_phi)
{
  LevelData<FArrayBox> averageLD(a_phi.disjointBoxLayout(), a_phi.nComp(), IntVect::Zero); // no ghost cell as this causes issues in parallel
  horizontallyAverage(averageLD, a_phi, averageVector);
}

void AMRLevelMushyLayer::horizontallyAverage(LevelData<FArrayBox>& a_averaged, LevelData<FArrayBox>& a_phi)
{
  Vector<Real> discardVector;
  horizontallyAverage(a_averaged, a_phi, discardVector);
}

void AMRLevelMushyLayer::horizontallyAverage(LevelData<FArrayBox>& a_averaged, LevelData<FArrayBox>& a_phi, Vector<Real>& globalAveraged)
{
  Box domBox = m_problem_domain.domainBox();

  IntVect smallEnd = domBox.smallEnd();
  int z_low = smallEnd[1];

  //  Real length = domBox.hiVect()[1] + 2 - domBox.loVect()[1];
  Real length = m_numCells[SpaceDim-1]+1; // need an extra cell here
  //  Real width = (domBox.hiVect()[0] + 1 - domBox.loVect()[0])*m_dx;
  Real width = m_domainWidth;
  Vector<Real> averaged(length, 0.0);

  globalAveraged.resize(length, 0.0);
  Vector<Vector<Real> > allAveraged(numProc());

  for (int i = 0; i<numProc(); i++)
  {
    allAveraged[i].resize(length);
    allAveraged[i].assign(0.0); //= new Vector<Real>(length+2, 0.0);
  }

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
  {
    Box box = m_grids[dit];

    box &= m_problem_domain.domainBox();
    for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      int z_i = iv[SpaceDim-1];
      averaged[z_i-z_low] += a_phi[dit](iv)*(m_dx/width);
    }

  }

  //  pout() << "HorizontallyAverage - got local averages, now broadcast/gather" << endl;

  // Broadcast/gather to compute averages over whole domain
  int srcProc = 0;

  gather(allAveraged, averaged, srcProc);

  if (procID() == srcProc)
  {
    for (int ivec = 0; ivec<numProc(); ivec++)
    {
      for (int y_i = 0; y_i < globalAveraged.size(); y_i++)
      {
        globalAveraged[y_i] += (allAveraged[ivec])[y_i];
      }
    }


  }

  //  pout() << "HorizontallyAverage - completed gather" << endl;

  broadcast(globalAveraged, srcProc);

  //  pout() << "HorizontallyAverage - complete broadcast" << endl;

  for (DataIterator dit = m_scalarNew[m_bulkConcentration]->dataIterator(); dit.ok(); ++dit)
  {
    Box box = m_grids[dit];
    for (int idir=0; idir<SpaceDim; idir++)
    {
      box &= m_problem_domain.domainBox();
      for (BoxIterator bit(box); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        int z_i = iv[SpaceDim-1];
        z_i = z_i - z_low;
        a_averaged[dit](iv) = globalAveraged[z_i];
      }
    }
  }
}

void AMRLevelMushyLayer::computeAdvectionVelSourceTerm(LevelData<FArrayBox>& a_src)
{
  CH_TIME("AMRLevelMushyLayer::computeAdvectionVelSourceTerm");

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::computeAdvectionVelSourceTerm"            << endl;
  }

  ParmParse ppMain("main");

  IntVect ivGhost = a_src.ghostVect();
  LevelData<FArrayBox> darcy(m_grids, SpaceDim, ivGhost);
  LevelData<FArrayBox> viscous(m_grids, SpaceDim, ivGhost);
  LevelData<FArrayBox> buoyancy(m_grids, SpaceDim, ivGhost);
  LevelData<FArrayBox> pressure(m_grids, SpaceDim, ivGhost);
  LevelData<FArrayBox> temperature(m_grids, 1, ivGhost);
  LevelData<FArrayBox> liquidConc(m_grids, 1, ivGhost);
  LevelData<FArrayBox> porosity(m_grids, 1, ivGhost);
  LevelData<FArrayBox> pressureScale(m_grids, 1, ivGhost);

  //  Real old_time = m_time - m_dt;
  //  Real src_time = m_time - m_dt; //old_time;
  Real src_time = m_time - m_dt; //old_time;
  //	Real half_time = m_time - m_dt/2;

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

  //  Real minChi = ::computeMin(*m_scalarNew[m_porosity], NULL, 1);
  Real advVelsrcChiLimit = m_solidPorosity; //1e-10

  ppMain.query("advVelSrcChiLimit", advVelsrcChiLimit);
  setVelZero(velOld, advVelsrcChiLimit);

  //Apply BCs to grad(P) term - extrapolation
  DataIterator dit = m_grids.dataIterator();
  BCHolder bcExtrap = m_physBCPtr->extrapFuncBC();

  //  m_projection.gradPiBCs(pressure, true);
  fillPressureSrcTerm(pressure, pressureScale, src_time-m_dt/2,
                      false); // false - make sure we use Pi (true means use phi)

  ParmParse ppAdvsrc("advSrc");
  bool pressureSrc = false;
  bool darcySrc = true;
  bool viscousSrc = true;
  bool buoyancySrc = true;

  ppAdvsrc.query("pressure", pressureSrc);
  ppAdvsrc.query("darcy", darcySrc);
  ppAdvsrc.query("viscous", viscousSrc);
  ppAdvsrc.query("buoyancy", buoyancySrc);

  if (m_time <= m_skipTrickySourceTerm)
  {
    pout() << "AMRLevelMushyLayer::computeAdvectionVelSourceTerm - time = " << m_time << " < " << m_skipTrickySourceTerm;
    pout() << ", skipping darcy and viscous src terms";
    darcySrc = false;
    viscousSrc = false;
  }

  //Calculate Laplacian(U)
  if (m_isViscous)
  {
    bool recomputeLapVel = true;

    bool allowLaggedLapVel = false;
    ppAdvsrc.query("allow_lagged_lap_vel", allowLaggedLapVel);
    if (m_level > 0 && allowLaggedLapVel)
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
        (*m_vectorNew[m_advSrcLapU])[dit].copy(viscous[dit]);
      }

    }
    else
    {
      for (dit.reset(); dit.ok(); ++dit)
      {
        //      m_vectorNew[m_advSrcLapU]->copyTo(viscous);
        viscous[dit].copy((*m_vectorNew[m_advSrcLapU])[dit]);
        //        int temp=0;
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

      // This has already been done
      //      pressure[dit].mult(porosity[dit],
      //                         pressure[dit].box(), 0, dir);

    }

    darcy[dit].mult(m_parameters.m_darcyCoeff);
    viscous[dit].mult(m_parameters.m_viscosityCoeff);

  }

  // Finally, set darcy src term=0 near CF interface
  setVelZero(darcy, advVelsrcChiLimit, 2); // set to zero up to 2 cells from interface

  //    fillBuoyancy(buoyancy[dit], temperature[dit], liquidConc[dit], porosity[dit],  (*m_vectorNew[m_bodyForce])[dit]);

  for (dit.reset(); dit.ok(); ++dit)
  {
    a_src[dit].setVal(0.0);

    if (viscousSrc)
    {
      a_src[dit] += viscous[dit];
    }

    if (pressureSrc)
    {
      a_src[dit] -= pressure[dit];
    }

    if (buoyancySrc)
    {
      a_src[dit] += buoyancy[dit];
    }

    if (darcySrc)
    {
      a_src[dit] -= darcy[dit];
    }



  } // end dataiterator

  // Add in extra (u.grad(porosity)/porosity^2)u term if needed
  if (m_advectionMethod == m_porosityOutsideAdvection)
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
  else if (m_advectionMethod == m_porosityInAdvection)
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

  //Get CF BCs
  // Why weren't we doing this before???
  // This shouldn't be necessary - CF bcs should be set on fields as above
  bool fillGhost=false;
  ppAdvsrc.query("fillGhost", fillGhost);
  if (m_level > 0 && fillGhost)
  {
    fillVectorField(a_src, src_time, m_advectionSrc); // this should just fill ghost cells
  }



  setVelZero(a_src, advVelsrcChiLimit);

  a_src.exchange();

  for (dit.reset(); dit.ok(); ++dit)
  {
    // This source term is centered at step n, so store it in vectorOld.
    (*m_vectorOld[m_advectionSrc])[dit].copy(a_src[dit], 0, 0, SpaceDim);

    // Also need to store in vectorNew so that we can interpolate in time when subcycling
    (*m_vectorNew[m_advectionSrc])[dit].copy(a_src[dit], 0, 0, SpaceDim);
  }

//  if (s_verbosity >= 1)
//    {
//      DisjointBoxLayout* finerGridsPtr = NULL;
//      int nRefFine = -1;
//      Real sumRHS = computeSum(*m_vectorNew[m_advectionSrc], finerGridsPtr,
//                               nRefFine, m_dx, m_vectorNew[m_advectionSrc]->interval());
//      pout() << "  Advection Velocity Src (level " << m_level << ") -- sum(RHS) = " << sumRHS << endl;
//    }

  if (crseVelPtr != NULL)
  {
    delete crseVelPtr;
    crseVelPtr = NULL;
  }

}

void AMRLevelMushyLayer::fillBuoyancy(LevelData<FArrayBox>& a_buoyancy, Real a_time)
{
  IntVect ghost = a_buoyancy.ghostVect();

  if (m_level == 0)
  {
    DisjointBoxLayout dbl = a_buoyancy.disjointBoxLayout();


    LevelData<FArrayBox> temperature(dbl, 1, ghost);
    LevelData<FArrayBox> liquidConc(dbl, 1, ghost);
    LevelData<FArrayBox> porosity(dbl, 1, ghost);

    fillScalars(temperature,a_time,m_temperature, true, true);
    fillScalars(liquidConc,a_time,m_liquidConcentration, true, true);
    fillScalars(porosity,a_time,m_porosity, true, true);

    fillBuoyancy(a_buoyancy, temperature, liquidConc, porosity);
  }
  else
  {
    AMRLevelMushyLayer* mlCrse = getCoarserLevel();
    LevelData<FArrayBox> crseBuoyancy(mlCrse->m_grids, SpaceDim, ghost);
    mlCrse->fillBuoyancy(crseBuoyancy, a_time);

    // Refine back to this level
    m_fineInterpVector.interpToFine(a_buoyancy, crseBuoyancy);

    // Fill ghost cells
    PiecewiseLinearFillPatch filpatcher(m_grids, mlCrse->m_grids,
                                        a_buoyancy.nComp(), mlCrse->m_problem_domain, mlCrse->m_ref_ratio, ghost[0],
                                        false);

    Interval scalComps = a_buoyancy.interval();

    filpatcher.fillInterp(a_buoyancy, crseBuoyancy, crseBuoyancy,
                          0,
                          0, 0, scalComps.size());

  }
}

void AMRLevelMushyLayer::fillBuoyancy(LevelData<FArrayBox>& a_buoyancy,
                                      LevelData<FArrayBox>& a_temperature,
                                      LevelData<FArrayBox>& a_liquidConc,
                                      LevelData<FArrayBox>& a_porosity)
{
  ParmParse ppMain("main");
  Real zero_time = -1;
  ppMain.query("turn_off_buoyancy_time", zero_time);

  for (DataIterator dit = a_buoyancy.dataIterator(); dit.ok(); ++dit)
  {
    if (zero_time >= 0 && m_time > zero_time)
    {
      a_buoyancy[dit].setVal(0);
    }
    else
    {
      fillBuoyancy(a_buoyancy[dit], a_temperature[dit], a_liquidConc[dit], a_porosity[dit],
                   (*m_vectorNew[m_bodyForce])[dit]);
    }
  }
}

void AMRLevelMushyLayer::fillBuoyancy(FArrayBox& buoyancy,FArrayBox& temperature, FArrayBox& liquidConc,
                                      FArrayBox& porosity,
                                      FArrayBox& bodyForce)
{
  int gravityDir= SpaceDim-1;

  // buoyancy = porosity * (Ra_T * theta - Ra_c*Theta_l) e_z
  buoyancy.setVal(0.0);


  buoyancy.plus(temperature, m_parameters.m_buoyancyTCoeff, 0, gravityDir, 1);
  buoyancy.plus(liquidConc, -m_parameters.m_buoyancySCoeff, 0, gravityDir, 1);
  buoyancy.plus(bodyForce, 1, gravityDir, gravityDir, 1);

  buoyancy.mult(porosity, buoyancy.box(), 0, gravityDir, 1);
}

void AMRLevelMushyLayer::backupTimestep()
{
  // Make sure we copy ghost cells as well

  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::backupTimestep " << endl;
  }

  for (DataIterator dit = m_scalarNew[0]->dataIterator(); dit.ok(); ++dit)
  {
    for (int var = 0; var < m_numScalarVars; var++)
    {
      (*m_scalarRestart[var])[dit].copy((*m_scalarNew[var])[dit]);
    }

    for (int var = 0; var < m_numVectorVars; var++)
    {
      (*m_vectorRestart[var])[dit].copy((*m_vectorNew[var])[dit]);
    }
  }

  m_projectionBackup.copyPressure(m_projection);


}

void AMRLevelMushyLayer::restartTimestepFromBackup(bool ignorePressure)
{
  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::restartTimestepFromBackup"        << endl;
  }
  Vector<int> restartVars;

  // This should be all we need:
  restartVars.push_back(m_enthalpy);
  restartVars.push_back(m_bulkConcentration);
  restartVars.push_back(m_pressure);
  restartVars.push_back(m_lambda);

  Vector<int> vectRestartVars;
  vectRestartVars.push_back(m_fluidVel);

  // Make sure we copy ghost cells
  for (DataIterator dit = m_scalarNew[0]->dataIterator(); dit.ok(); ++dit)
  {

    for (int i = 0; i < restartVars.size(); i++)
    {
      int var = restartVars[i];
      if (! (var == m_pressure && ignorePressure))
      {
        (*m_scalarNew[var])[dit].copy((*m_scalarRestart[var])[dit]);
        (*m_scalarOld[var])[dit].copy((*m_scalarRestart[var])[dit]);
      }
    }

    for (int i = 0; i < vectRestartVars.size(); i++)
    {
      int var = vectRestartVars[i];
      (*m_vectorNew[var])[dit].copy((*m_vectorRestart[var])[dit]);
      (*m_vectorOld[var])[dit].copy((*m_vectorRestart[var])[dit]);

    }
  }

  if (!ignorePressure)
  {
    m_projection.copyPressure(m_projectionBackup);
  }

  // Go back to start of timestep
  //  m_time = m_time - m_dt;

  // Halve timestep
  //  m_dt = m_dt/2;
}

int AMRLevelMushyLayer::multiCompAdvectDiffuse(LevelData<FArrayBox>& a_phi_old, LevelData<FArrayBox>& a_phi_new,
                                               LevelData<FArrayBox>& a_src, bool doFRupdates, bool computeAdvectiveSrc)
{
  CH_TIME("AMRLevelMushyLayer::advectDiffuseScalar");
  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::advectDiffuseScalar: " << m_level << endl;
  }

  int exitStatus = 0;
  ParmParse ppMain("main");
  int numComp = a_src.nComp();

  LevelData<FArrayBox> full_src(m_grids, numComp, IntVect::Zero);
  LevelData<FluxBox> totalAdvectiveFlux(m_grids, numComp, IntVect::Unit); // Need ghost vector for dealing with patches
  DataIterator dit = m_grids.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    full_src[dit].copy(a_src[dit]);
    totalAdvectiveFlux[dit].setVal(0.0);
  }

  if (computeAdvectiveSrc)
  {
    // Always do multi comp advection...
    //    bool doMultiCompAdvection = true;

    computeScalarAdvectiveSrcHC(full_src, totalAdvectiveFlux, doFRupdates);


    bool skipSalt = false;
    ppMain.query("skipSaltUpdate", skipSalt);
    if (skipSalt)
    {
      for (dit.reset(); dit.ok(); ++dit)
      {
        full_src[dit].setVal(0.0, 1);
      }
      // this needs to be done before operators are defined - just set lewis=1e300 in inputs file
//      m_parameters.m_saltDiffusionCoeff = 0;
    }

  } // end if compute advective src

  full_src.copyTo(Interval(0,0), *m_scalarNew[m_enthalpySrc], Interval(0,0));
  full_src.copyTo(Interval(1,1), *m_scalarNew[m_saltEqnSrcGodunov], Interval(0,0));

  // Option here to not update enthalpy and salinity
  // (useful for debugging)
  bool skipUpdate = false;

  ppMain.query("skipHCUpdate", skipUpdate);
  if (skipUpdate)
  {
    pout() << "Skipping HC update" << endl;
    return 0;
  }

  // Set up coarse-fine boundary conditions

  LevelFluxRegister* coarserFRPtr = NULL;
  LevelFluxRegister* finerFRPtr = NULL;
  Real tCoarserOld, tCoarserNew;

  LevelData<FArrayBox>* coarserDataOldPtr = NULL;
  LevelData<FArrayBox>* coarserDataNewPtr = NULL;

  tCoarserOld = 0.0;
  tCoarserNew = 0.0;

  // If a coarser level exists
  if (m_hasCoarser)
  {
    AMRLevelMushyLayer* coarserPtr = getCoarserLevel();

    // Recall that my flux register goes between my level and the next
    // finer level
    coarserFRPtr = &(*coarserPtr->m_fluxRegHC);

    coarserDataOldPtr = new LevelData<FArrayBox>(coarserPtr->m_grids, 2, IntVect::Unit);
    coarserDataNewPtr = new LevelData<FArrayBox>(coarserPtr->m_grids, 2, IntVect::Unit);

    tCoarserNew = coarserPtr->m_time;
    tCoarserOld = tCoarserNew - coarserPtr->m_dt;

    coarserPtr->fillHC(*coarserDataOldPtr,tCoarserOld );
    coarserPtr->fillHC(*coarserDataNewPtr,tCoarserNew );
  }

  // If a finer level exists
  if (m_hasFiner)
  {
    // Recall that my flux register goes between my level and the next
    // finer level

    // Check flux register is defined - if we are restarting,
    // then it might not be defined yet on finer (currently unused) levels
    if (m_fluxRegHC)
    {
      finerFRPtr = &(*m_fluxRegHC);
    }
  }

  // Only update FRs if we've converged
  // Haven't sorted multi component flux registers yet
  if (doFRupdates)
  {
    incrementHCFluxRegisters(totalAdvectiveFlux, m_dt);
  }
  else
  {
    finerFRPtr = NULL;
    coarserFRPtr = NULL;
  }// end if do advective flux reg updates

  if (s_verbosity > 5)
  {
    pout() << "multiCompAdvectDiffuse - call backward euler on level " << m_level << endl;
    pout() << "multiCompAdvectDiffuse - has coarser? " << m_hasCoarser << ", has finer? " << m_hasFiner << endl;
  }

  Real old_time = m_time-m_dt;

  BaseLevelHeatSolver<LevelData<FArrayBox>, FluxBox, LevelFluxRegister>* baseLevBE = NULL;

  if (m_timeIntegrationOrder == 2)
  {
    //       MayDay::Error("multiCompAdvectDiffuse - TGA not implemented yet");
    //exitStatus = TGAUpdateScalar(a_var, a_src, converged);

    s_enthalpySalinityTGA->updateSoln(a_phi_new,
                                      a_phi_old, full_src, finerFRPtr, coarserFRPtr,
                                      coarserDataOldPtr, coarserDataNewPtr, old_time, tCoarserOld,
                                      tCoarserNew, m_dt, m_level, false); //false - don't zero phi

    baseLevBE = dynamic_cast<BaseLevelHeatSolver<LevelData<FArrayBox>, FluxBox, LevelFluxRegister> * > (&(*s_enthalpySalinityTGA));
  }
  else
  {

    if (ppMain.contains("noMultigrid"))
    {


      int num_iter=100;
      ppMain.query("noMultigridIter", num_iter);
      a_phi_old.copyTo(a_phi_new);


      LevelData<FArrayBox> thisSrc(m_grids, 2, IntVect::Unit);
      full_src.copyTo(thisSrc);
      for (DataIterator dit = thisSrc.dataIterator(); dit.ok(); ++dit)
      {
        thisSrc[dit].mult(m_dt);
        thisSrc[dit].plus(a_phi_old[dit]);
      }

      for (int i=0; i < num_iter; i++)
            {

        // Need to define solvers first and foremost
//        defineSolvers(m_time);

      // Try replacing with single level relax solver
      RefCountedPtr<AMRNonLinearMultiCompOp> amrpop = RefCountedPtr<AMRNonLinearMultiCompOp>(
            (AMRNonLinearMultiCompOp*) HCOpFact->AMRnewOp(m_problem_domain));

      // Solving (1-dt Op) HC^{n+1} = HC^{n}
      amrpop->setAlphaAndBeta(1.0, -m_dt);

      LevelData<FArrayBox> a_res(m_grids, 2, IntVect::Unit);

      amrpop->relax(a_phi_new, thisSrc, 1000);

      amrpop->residual(a_res, a_phi_new, thisSrc, false);
      Real maxRes = ::computeNorm(a_res, NULL, 1, m_dx, Interval(0, 1), 0);
      pout() << "  Max residual = " << maxRes << endl;

      if (maxRes < 1e-10)
      {
        pout() << "   Converged " << endl;
        break;

      }

      }
    }
    else
    {

    s_enthalpySalinityBE->updateSoln(a_phi_new,
                                     a_phi_old, full_src, finerFRPtr, coarserFRPtr,
                                     coarserDataOldPtr, coarserDataNewPtr, old_time, tCoarserOld,
                                     tCoarserNew, m_dt, m_level, false); //false - don't zero phi

    baseLevBE = dynamic_cast<BaseLevelHeatSolver<LevelData<FArrayBox>, FluxBox, LevelFluxRegister> * > (&(*s_enthalpySalinityBE));
    }

  }

  if (s_verbosity > 5)
  {
    pout() << "multiCompAdvectDiffuse -  updated solution" << endl;
  }

  a_phi_new.exchange();

  Real residual = 0;

#ifdef CH_FORK
  if (baseLevBE != NULL)
  {
    exitStatus = baseLevBE->exitStatus();
    residual = baseLevBE->finalResidual();
    int num_iter =  baseLevBE->numMGiterations();

    pout() << "  HC solve finished with exit status " << exitStatus << ", solver residual = " << residual << ", num MG iterations = " << num_iter << endl;

  }
#endif

  // Clean up
  if (coarserDataNewPtr != NULL)
  {
    delete coarserDataNewPtr;
    coarserDataNewPtr = NULL;
  }

  if (coarserDataOldPtr != NULL)
  {
    delete coarserDataOldPtr;
    coarserDataOldPtr = NULL;
  }


  return exitStatus;
}

void AMRLevelMushyLayer::incrementHCFluxRegisters(LevelData<FluxBox>& flux, Real fluxMult)
{
  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::incrementFluxRegisters on level " << m_level << endl;
  }
  // Do flux register updates here!

  Interval refluxComps = flux.interval();
  Interval scalComps = flux.interval();

  DataIterator dit = flux.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {

    FluxBox& localFlux = flux[dit];

    if (m_level > 0)
    {
      LevelFluxRegister& crseFR = (*getCoarserLevel()->m_fluxRegHC);
      CH_assert (crseFR.isDefined());
      for (int dir=0; dir<SpaceDim; dir++)
      {
        crseFR.incrementFine(localFlux[dir],
                             fluxMult, dit(),
                             scalComps, refluxComps,
                             dir);

      } // end loop over directions
    } // end if coarser level exists

    if (!finestLevel())
    {
      CH_assert(m_fluxRegHC->isDefined());
      for (int dir=0; dir<SpaceDim; ++dir)
      {
        m_fluxRegHC->incrementCoarse(localFlux[dir],
                                     fluxMult, dit(),
                                     scalComps, refluxComps,
                                     dir);

      } // end loop over directions
    } // end if finer level exists
  } // end loop over grids for flux register updates
}



void AMRLevelMushyLayer::computeTotalAdvectiveFluxes(LevelData<FluxBox>& edgeScalTotal)
{
  IntVect advect_grow = m_numGhostAdvection*IntVect::Unit;
  Real old_time = m_time - m_dt;
  int numComp = 2;

  CH_assert(edgeScalTotal.nComp() == numComp);

  ParmParse ppMain("main");

  IntVect ghostVect = IntVect::Unit;
  IntVect fluxGhostVect = edgeScalTotal.ghostVect();

  LevelData<FluxBox> edgeScalFluidAdv(m_grids, numComp, fluxGhostVect);
  LevelData<FluxBox> edgeScalFrameAdvection(m_grids, numComp, fluxGhostVect);

  // Need grown version of this
  LevelData<FArrayBox> TCl_old(m_grids, numComp, advect_grow);
  LevelData<FArrayBox> HC_old(m_grids, numComp, advect_grow);


  fillHC(HC_old, old_time,
         true,  // fill interior?
         (CFinterpOrder_advection==2)); // do quadratic interpolation at CF boundaries?
  fillTCl(TCl_old, old_time,
          true,  // fill interior?
          (CFinterpOrder_advection==2)); //  do quadratic interpolation at CF boundaries? - need this for corner cells that border CF and domain boundaries


  // Compute frame advection src term
  computeScalarAdvectiveFluxMultiComp(edgeScalFrameAdvection, m_frameAdvVel,
                                      m_patchGodHC, HC_old,
                                      old_time, m_dt);

  // todo - stop forcing the use of single component solvers
  bool allowMulticompAdvection = false;
  ppMain.query("allowMulticompAdvection", allowMulticompAdvection);

  if (s_reflux_enthalpy == s_reflux_concentration && allowMulticompAdvection)
  {
    computeScalarAdvectiveFluxMultiComp(edgeScalFluidAdv, m_advVel,
                                        m_patchGodTSl, TCl_old,
                                        old_time, m_dt);

  }
  else
  {
    // we have to consider enthalpy and bulk concentration src terms separately if we're doing
    // reflux differently between them. Specifically, if we're not refluxing then we shouldn't
    // add the freestream preservation term to the advection velocity.

    LevelData<FluxBox> advVel_noReflux(m_advVel.disjointBoxLayout(), 1, m_advVel.ghostVect());
    LevelData<FluxBox> advVel_reflux(m_advVel.disjointBoxLayout(), 1, m_advVel.ghostVect());
    LevelData<FArrayBox> vel(m_grids, SpaceDim, advect_grow);
    LevelData<FArrayBox> diffusiveSrcMultiComp(m_grids, 2, edgeScalFluidAdv.ghostVect()+IntVect::Unit); // this needs the ghost cell to cover slope boxes
    LevelData<FArrayBox> diffusiveSrcSingleComp(m_grids, 1, edgeScalFluidAdv.ghostVect()+IntVect::Unit); // this needs the ghost cell to cover slope boxes

    LevelData<FArrayBox> scalarOld(m_grids, 1, advect_grow);

    LevelData<FluxBox> a_HeatFlux(edgeScalFluidAdv.disjointBoxLayout(), 1, edgeScalFluidAdv.ghostVect());
    LevelData<FluxBox> a_ConcFlux(edgeScalFluidAdv.disjointBoxLayout(), 1, edgeScalFluidAdv.ghostVect());

    m_advVel.copyTo(advVel_noReflux);
    m_advVel.copyTo(advVel_reflux);

    for (DataIterator dit = edgeScalFluidAdv.dataIterator(); dit.ok(); ++dit)
    {
      advVel_noReflux[dit].copy(m_advVel[dit]);
      advVel_reflux[dit].copy(m_advVel[dit]);


    }

    m_projection.removeFreestreamCorrection(advVel_noReflux);

    LevelData<FluxBox> *enthalpy_advVel, *conc_advVel;

    computeScalDiffusion(diffusiveSrcMultiComp, old_time);

    if (s_reflux_enthalpy)
    {
      enthalpy_advVel = &advVel_reflux;
    }
    else
    {
      enthalpy_advVel = &advVel_noReflux;
    }

    if (s_reflux_concentration)
    {
      conc_advVel = &advVel_reflux;
    }
    else
    {
      conc_advVel = &advVel_noReflux;
    }

    EdgeToCell(*enthalpy_advVel, vel);
    diffusiveSrcMultiComp.copyTo(Interval(0,0), diffusiveSrcSingleComp, Interval(0,0));
    TCl_old.copyTo(Interval(0,0), scalarOld, Interval(0,0));

    for (DataIterator dit = edgeScalFluidAdv.dataIterator(); dit.ok(); ++dit)
    {
      diffusiveSrcSingleComp[dit].copy(diffusiveSrcMultiComp[dit], 0, 0, 1);
      scalarOld[dit].copy(TCl_old[dit], 0, 0, 1);
    }

    computeScalarAdvectiveFlux(a_HeatFlux, scalarOld, *enthalpy_advVel,
                               m_totalAdvVel,
                               vel,
                               diffusiveSrcSingleComp, m_patchGodT,
                               old_time, m_dt);


    diffusiveSrcMultiComp.copyTo(Interval(1,1), diffusiveSrcSingleComp, Interval(0,0));
    EdgeToCell(*conc_advVel, vel);
    TCl_old.copyTo(Interval(1,1), scalarOld, Interval(0,0));

    for (DataIterator dit = edgeScalFluidAdv.dataIterator(); dit.ok(); ++dit)
    {
      diffusiveSrcSingleComp[dit].copy(diffusiveSrcMultiComp[dit], 1, 0, 1);
      scalarOld[dit].copy(TCl_old[dit], 1, 0, 1);
    }

    computeScalarAdvectiveFlux(a_ConcFlux, scalarOld, *conc_advVel,
                               m_totalAdvVel,
                               vel,
                               diffusiveSrcSingleComp, m_patchGodSl,
                               old_time, m_dt);

    // Now copy predicted values to edgeScalFluidAdv
    for (DataIterator dit = edgeScalFluidAdv.dataIterator(); dit.ok(); ++dit)
    {
      edgeScalFluidAdv[dit].setVal(0.0);
      edgeScalFluidAdv[dit].plus(a_HeatFlux[dit], a_HeatFlux[dit].box(), 0, 0, 1);
      edgeScalFluidAdv[dit].plus(a_ConcFlux[dit], a_ConcFlux[dit].box(), 0, 1, 1);
    }

  }
  // Combine the two fluxes
  for (DataIterator dit = edgeScalFluidAdv.dataIterator(); dit.ok(); ++dit)
  {
    edgeScalTotal[dit].setVal(0.0);

    edgeScalTotal[dit].plus(edgeScalFrameAdvection[dit], edgeScalFrameAdvection[dit].box(),
                            0,0,numComp);
    edgeScalTotal[dit].plus(edgeScalFluidAdv[dit], edgeScalFluidAdv[dit].box(),
                            0,0,numComp);

    // Also copy the vertical component to an farray box
    //    (*m_scalarNew[m_FsVertFluid])[dit].copy(edgeScalFluidAdv[dit], 1, 0, 1);
    //    (*m_scalarNew[m_FsVertFrame])[dit].copy(edgeScalFrameAdvection[dit], 1, 0, 1);
    EdgeToCell(edgeScalFluidAdv[dit], 1, (*m_scalarNew[m_FsVertFluid])[dit], 0, 1);
    EdgeToCell(edgeScalFrameAdvection[dit], 1, (*m_scalarNew[m_FsVertFrame])[dit], 0, 1);

    for (int dir=0; dir <SpaceDim; dir++)
    {
      EdgeToCell(edgeScalFluidAdv[dit], 1, (*m_vectorNew[m_FsFluid])[dit], dir, dir);
    }

    (*m_scalarNew[m_FsVertFluid])[dit].mult(-1);
    (*m_scalarNew[m_FsVertFrame])[dit].mult(-1);
  }

  //  EdgeToCell()

  if (edgeScalTotal.ghostVect()[0] > 0)
  {
    CornerCopier cornerCopier;
    cornerCopier.define(m_grids, m_grids, m_problem_domain, edgeScalTotal.ghostVect(), true );
    edgeScalTotal.exchange(cornerCopier);
  }
  else
  {
    edgeScalTotal.exchange(edgeScalTotal.interval());
  }

  // Set some analytic source term if needed


  bool analyticSource = false;
  ppMain.query("analyticSourceTerm", analyticSource);
  if (analyticSource)
  {
    for (DataIterator dit = edgeScalTotal.dataIterator(); dit.ok(); ++dit)
    {
      FluxBox& thisFlux = edgeScalTotal[dit];

      for (BoxIterator bit = BoxIterator(thisFlux.box()); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, loc, m_dx);

        // Horizontal
        thisFlux[0](iv, 0) =  0; // H comp
        thisFlux[0](iv, 1) = 0.0; // C comp

        // Vertical
        //          thisFlux[1](iv, 0) =  -m_parameters.nonDimVel*loc[1]*loc[1]; // H comp
        thisFlux[1](iv, 0) =  -m_parameters.nonDimVel*loc[1]; // H comp
        thisFlux[1](iv, 1) = 0.0; // C comp
      }

    }
  }

  // These things are useful for debugging
  //  LevelData<FArrayBox> fluidAdvSrc(m_grids, numComp, IntVect::Zero); // want this for debugging

  //  LevelData<FArrayBox> frameAdvSrc(m_grids, numComp, IntVect::Zero); // want this for debugging

  //  Divergence::levelDivergenceMACMultiComp(fluidAdvSrc, edgeScalFluidAdv, m_dx);

  //  Divergence::levelDivergenceMACMultiComp(frameAdvSrc, edgeScalFrameAdvection, m_dx);
  //  int temp = 0;
}

void AMRLevelMushyLayer::computeScalarAdvectiveSrcHC(LevelData<FArrayBox>& a_src,
                                                     LevelData<FluxBox>& edgeScalTotal,
                                                     bool converged)
{
  int numComp = 2;
  CH_assert(a_src.nComp() == numComp);

  DataIterator dit = a_src.dataIterator();

  // Compute fluxes due to fluid advection and frame advection
  computeTotalAdvectiveFluxes(edgeScalTotal);
  //  edgeScalTotal.exchange();

  // The source term is the divergence of the fluxes
  LevelData<FArrayBox> advectiveSrc(m_grids, numComp, IntVect::Zero);
  Divergence::levelDivergenceMACMultiComp(advectiveSrc, edgeScalTotal, m_dx);

  for (dit.reset(); dit.ok(); ++dit)
  {
    a_src[dit].plus(advectiveSrc[dit],
                    -1.0,  // scale
                    0, // src comp
                    0, // dest comp
                    numComp); // num comps

    //    a_src[dit].mult(m_parameters.m_advectionCoeff);

  }

}


void AMRLevelMushyLayer::computeScalarAdvectiveFlux(LevelData<FluxBox>& a_edgeScal, int a_advectionVar, int a_diffusionVar,
                                                    LevelData<FluxBox>& a_advVel,
                                                    Real a_old_time, Real a_dt)
{
  // Need grown version of this
  IntVect advect_grow = m_numGhostAdvection*IntVect::Unit;
  LevelData<FArrayBox> scalar_advection_old(m_grids, 1, advect_grow);
  LevelData<FArrayBox> vel(m_grids, SpaceDim, advect_grow);
  LevelData<FArrayBox> diffusiveSrc(m_grids, 1, a_edgeScal.ghostVect()+IntVect::Unit); // this needs the ghost cell to cover slope boxes

  EdgeToCell(a_advVel, vel);
  fillScalars(scalar_advection_old, a_old_time, a_advectionVar,
              true, //do interior
              (CFinterpOrder_advection==2) // quad interp - this seems to fix previous issues at insulating side walls
  );

  // Make diffusive source
  bool doDiffusionSrc = true;
  if (a_diffusionVar > -1 && doDiffusionSrc)
  {
    LevelData<FArrayBox>* crseScalarDiffusion = NULL;

    if (m_level > 0)
    {
      // allocate crseBC info
      AMRLevelMushyLayer* mlCrse = getCoarserLevel();

      const DisjointBoxLayout& crseGrids = mlCrse->m_grids;
      crseScalarDiffusion = new LevelData<FArrayBox>(crseGrids,1);

      mlCrse->fillScalars(*crseScalarDiffusion, a_old_time, a_diffusionVar, true);
    }

    // Get something like grad^2(T) or div(chi dot grad S_l)
    computeScalDiffusion(a_diffusionVar, diffusiveSrc, a_old_time);
  }
  else
  {
    setValLevel(diffusiveSrc, 0.0);
  }

  // determine if we have inflow or outflow
  computeInflowOutflowAdvVel();

  // Compute advective flux
  computeScalarAdvectiveFlux(a_edgeScal, scalar_advection_old, a_advVel,
                             m_totalAdvVel,
                             vel,
                             diffusiveSrc, *m_patchGodScalars[a_advectionVar],
                             a_old_time, m_dt);
}

/// Compute temperature and salinity advective flux
void AMRLevelMushyLayer::computeScalarAdvectiveFluxMultiComp(LevelData<FluxBox>& a_edgeScal,
                                                             LevelData<FluxBox>& a_advVel,
                                                             PatchGodunov& a_patchGod,
                                                             LevelData<FArrayBox>& a_scalOld,
                                                             Real a_old_time, Real a_dt)
{
  // Need grown version of this
  int numComp = a_scalOld.nComp(); // 2 components - temperature and liquid salinity (or enthalpy and bulk concentration)

  CH_assert(a_edgeScal.nComp() == numComp);
  //  CH_assert(a_scalOld.nComp() == numComp);

  IntVect advect_grow = m_numGhostAdvection*IntVect::Unit;
  LevelData<FArrayBox> vel(m_grids, SpaceDim, advect_grow);
  LevelData<FArrayBox> diffusiveSrc(m_grids, numComp, a_edgeScal.ghostVect()+IntVect::Unit); // this needs the ghost cell to cover slope boxes

  EdgeToCell(a_advVel, vel);
  vel.exchange();

  // Make diffusive source
  bool doDiffusionSrc = true; // hard code this for now
  ParmParse ppMain("main");
  ppMain.query("diffusiveSrcForAdvection", doDiffusionSrc);
  if (doDiffusionSrc)
  {
    // Get something like grad^2(T) or div(chi dot grad S_l)
    computeScalDiffusion(diffusiveSrc, a_old_time);
  }
  else
  {
    setValLevel(diffusiveSrc, 0.0);
  }

  // determine if we have inflow or outflow
  computeInflowOutflowAdvVel();

  // Compute advective flux
  computeScalarAdvectiveFlux(a_edgeScal, a_scalOld, a_advVel,
                             m_totalAdvVel,
                             vel,
                             diffusiveSrc, a_patchGod,
                             a_old_time, m_dt);


}



void AMRLevelMushyLayer::advectScalar(const int a_scalarVar, const int a_advectionVar,
                                      LevelData<FluxBox>& a_advVel, bool doFRupdates)
{
  LevelData<FluxBox> flux(m_grids, 1);
  advectScalar(a_scalarVar,a_advectionVar,
               a_advVel, doFRupdates,
               flux);
}


void AMRLevelMushyLayer::advectScalar(const int a_scalarVar, const int a_advectionVar,
                                      LevelData<FluxBox>& a_advVel, bool doFRupdates,
                                      LevelData<FluxBox>& flux)
{
  LevelFluxRegister* coarserFRPtr = NULL;
  LevelFluxRegister* finerFRPtr = NULL;
  LevelData<FArrayBox>* coarserDataOldPtr = NULL;
  LevelData<FArrayBox>* coarserDataNewPtr = NULL;
  Real tCoarserOld, tCoarserNew;

  getCoarseScalarDataPointers(a_scalarVar,
                              &coarserDataOldPtr,  &coarserDataNewPtr,  // we don't use these two, but need them as dummy arguments
                              &coarserFRPtr, &finerFRPtr, // get the flux registers for the thing we're updating, a_scalarVar
                              tCoarserOld, tCoarserNew); // don't need these either, they're just dummy arguments

  if (!doFRupdates)
  {
    coarserFRPtr = NULL;
    finerFRPtr = NULL;
  }

  DataIterator dit(m_grids);

  // Get the flux of a_advectionVar, i.e. u*a_advectionVar
  computeScalarAdvectiveFlux(flux, a_advectionVar, -1, a_advVel, m_time-m_dt, m_dt); // -1 means no diffusive src

  // Make the source term, div(u*a_advectionVar)
  LevelData<FArrayBox> update(m_grids, 1);
  Divergence::levelDivergenceMAC(update, flux, m_dx);

  // Add the source term to the old time solution
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    update[dit].mult(m_dt);
    (*m_scalarNew[a_scalarVar])[dit] -= update[dit];

  }

  // Flux register updates
  if (doFRupdates)
  {
    Real scale = m_dt;
    updateScalarFluxRegister(a_scalarVar, flux, scale);

  }


}

void AMRLevelMushyLayer::updateScalarFluxRegister(int a_scalarVar, LevelData<FluxBox>& flux, Real scale)
{
  //  Real scale = m_dt;

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {


    FluxBox& f = flux[dit];

    for (int dir = 0; dir < SpaceDim; dir++)
    {
      FArrayBox& fluxDir = f[dir];
      if (this->hasCoarserLevel())
      {
        getCoarserLevel()->m_fluxRegisters[a_scalarVar]->incrementFine(fluxDir, scale, dit(),
                                                                       Interval(0,0), Interval(0,0), dir);
      }

      if (hasFinerLevel())
      {
        m_fluxRegisters[a_scalarVar]->incrementCoarse(fluxDir, scale, dit(),
                                                      Interval(0,0), Interval(0,0), dir);
      }

    }
  }
}




bool AMRLevelMushyLayer::finestLevel()
{
  return !hasFinerLevel();
}


void AMRLevelMushyLayer::computeDiagnostics()
{
  if (!m_computeDiagnostics)
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
  if ((m_parameters.physicalProblem == MushyLayerParams::m_HRL
      || m_parameters.physicalProblem == MushyLayerParams::m_rayleighBenard)
      && m_level == 0)
  {
    fillScalars(*m_scalarNew[m_temperature], m_time, m_temperature);

    Real Nu = ::computeNusselt(*m_scalarNew[m_temperature], *m_vectorNew[m_fluidVel],
                               m_dx, m_parameters,
                               m_domainWidth, m_domainHeight);

    m_diagnostics.addDiagnostic(Diagnostics::m_Nu, m_time, Nu);
  }
  else if (m_parameters.physicalProblem == MushyLayerParams::m_convectionMixedPorous
      && m_level == 0)
  {

    // Calculate average dT/dx at left and right boundaries

    //    Vector<LevelData<FArrayBox>* > T(nLevels);
    Vector<LevelData<FArrayBox>* > gradT(nLevels);
    Vector<LevelData<FArrayBox>* > AMRNu(nLevels);
    Vector<LevelData<FluxBox>* > gradEdgeT(nLevels);
    LevelData<FArrayBox>* crseT = NULL;
    //    LevelData<FArrayBox>* fineT = NULL;

    // Use domain flux registers to calculate Nu
    Vector<LevelDomainFluxRegister*> Tflux(nLevels);

    LevelDomainFluxRegister* fineFR = NULL;
    LevelDomainFluxRegister* coarseFR = NULL;

    for (int lev=0; lev<nLevels; lev++)
    {
      Tflux[lev] = new LevelDomainFluxRegister();
    }

    //    Real Nu = 0;
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
        coarseFR = NULL;
      }

      if (lev < (nLevels-1))
      {
        fineFR = Tflux[lev+1];
      }
      else
      {
        fineFR = NULL;
      }


      Tflux[lev]->define(ml->m_problem_domain, ml->m_grids,refRat[lev], ml->m_dx, fineFR, coarseFR, SpaceDim);

      AMRLevelMushyLayer* mlCrse = ml->getCoarserLevel();
      if (mlCrse)
      {
        crseT = mlCrse->m_scalarNew[m_temperature];
      }

      Box domainBox = ml->m_problem_domain.domainBox();

      // Ensure BCs are set
      ml->fillScalars(*ml->m_scalarNew[m_temperature], ml->m_time, m_temperature, false, true);

      Gradient::levelGradientMAC(*gradEdgeT[lev], (*ml->m_scalarNew[m_temperature]),
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

    m_diagnostics.addDiagnostic(Diagnostics::m_Nu, m_time, Nu);
    m_diagnostics.addDiagnostic(Diagnostics::m_NuLeft, m_time, Nuleft);
    m_diagnostics.addDiagnostic(Diagnostics::m_NuRight, m_time, Nuright);

    // Alternative calculation
//    fillScalars(*m_scalarNew[m_temperature], m_time, m_temperature);
//
//    Nu = ::computeNusselt(*m_scalarNew[m_temperature], *m_vectorNew[m_fluidVel],
//                                   m_dx, m_parameters,
//                                   m_domainWidth, m_domainHeight);

    // Compute max velocities at midpoints
    Box horizBox = ::adjCellHi(m_problem_domain.domainBox(), 1, 1);
    horizBox.shift(1, -int(m_numCells[1]/2));

    Real maxVertVel = 0.0;

    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& vel = (*m_vectorNew[m_fluidVel])[dit];

      Box b = m_grids[dit];
      b &= horizBox;

      if (b.size() >= IntVect::Unit)
      {
        // TODO: this won't work in parallel
        maxVertVel = max(maxVertVel, vel.max(b, 1));
      }
    }


    Box vertBox = ::adjCellHi(m_problem_domain.domainBox(), 0, 1);
    vertBox.shift(0, -int(m_numCells[0]/2));

    Real maxHorizVel = 0.0;
    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          Box b = m_grids[dit];
          b &= vertBox;

          FArrayBox& vel = (*m_vectorNew[m_fluidVel])[dit];
          if (b.size() >= IntVect::Unit)
               {
          // TODO: this won't work in parallel
          maxHorizVel = max(maxHorizVel, vel.max(b, 0));
               }
        }

    m_diagnostics.addDiagnostic(Diagnostics::m_maxVhalf, m_time, maxHorizVel);
    m_diagnostics.addDiagnostic(Diagnostics::m_maxUhalf, m_time, maxVertVel);

    for (int lev=0; lev<nLevels; lev++)
    {
      if (Tflux[lev] != NULL)
      {
        delete Tflux[lev];
        Tflux[lev] = NULL;
      }

      if (gradT[lev] != NULL)
      {
        delete gradT[lev];
        gradT[lev] = NULL;
      }

      if ( gradEdgeT[lev] != NULL)
      {
        delete gradEdgeT[lev] ;
        gradEdgeT[lev]  = NULL;
      }

      if ( AMRNu[lev]  != NULL)
      {
        delete  AMRNu[lev] ;
        AMRNu[lev]  = NULL;
      }

    }

  }

  // Composite Nusselt calculation

  // Compute increments to solute flux on each level
  int numComp = 2;

  IntVect fluxGhost = IntVect::Zero;

  LevelData<FluxBox> totalFlux(m_grids, numComp, IntVect::Unit);
  getTotalFlux(totalFlux);

  LevelData<FluxBox> totalFluxNoGhost(m_grids, numComp, IntVect::Zero);
  totalFlux.copyTo(totalFluxNoGhost);

  // new way of keeping track of fluxes
  Real scale = (1/m_dx)*m_dt;

  m_heatDomainFluxRegister.incrFlux(totalFluxNoGhost, scale, 0);
  m_saltDomainFluxRegister.incrFlux(totalFluxNoGhost, scale, 1);

  // Calculate horizontally averaged fluxes
  LevelData<FArrayBox> averageFrameFlux(m_grids, numComp);
  LevelData<FArrayBox> averageAdvectiveFlux(m_grids, numComp);
  LevelData<FArrayBox> averageDiffusiveFlux(m_grids, numComp);
  LevelData<FArrayBox> averageVerticalFlux(m_grids, numComp);

  horizontallyAverage(averageVerticalFlux, totalFluxNoGhost);

  // Copy across for writing to plot files
  averageVerticalFlux.copyTo(Interval(0, 0), *m_scalarNew[m_averageHeatFlux], Interval(0,0));
  averageVerticalFlux.copyTo(Interval(1, 1), *m_scalarNew[m_averageVerticalFlux], Interval(0,0));

  // Also want solute flux at each point in space written out
  int soluteComp = 1;
  for (DataIterator dit = totalFlux.dataIterator(); dit.ok(); ++dit)
  {
    FluxBox& flux = totalFlux[dit];
    FArrayBox& fab = flux[SpaceDim-1];
    Box b = flux.box();

    b.growDir(1, Side::Hi, -1); // need this because we also grab the flux in the cell one up

    Box b2 = (*m_scalarNew[m_verticalFlux])[dit].box();
    b &= b2; // ensure we don't try and fill cells which don't exist
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      IntVect ivUp = iv + BASISV(1);
      (*m_scalarNew[m_verticalFlux])[dit](iv) = 0.5*(fab(iv, soluteComp) + fab(ivUp, soluteComp));
    }
  }



  if (m_level == 0)
  {
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

    AMRLevelMushyLayer* ml = this;

    for(int lev = 0; lev < nLevels; lev++)
    {
      S_new[lev] = ml->m_scalarNew[m_bulkConcentration];
      H_new[lev] = ml->m_scalarNew[m_enthalpy];
      horizAvFs[lev] = new LevelData<FArrayBox>(ml->m_grids, 1, IntVect::Zero); // no ghost vector
      Fs_vert_diffusion[lev] = new LevelData<FArrayBox>(ml->m_grids, 1, IntVect::Zero); // no ghost vector
      Fs_vert_fluid[lev] = new LevelData<FArrayBox>(ml->m_grids, 1, IntVect::Zero); // no ghost vector
      Fs_vert_frame[lev] = new LevelData<FArrayBox>(ml->m_grids, 1, IntVect::Zero); // no ghost vector

      ml->m_scalarNew[m_averageVerticalFlux]->copyTo(*horizAvFs[lev]);
      ml->m_scalarNew[m_FsVertDiffusion]->copyTo(*Fs_vert_diffusion[lev]);
      ml->m_scalarNew[m_FsVertFluid]->copyTo(*Fs_vert_fluid[lev]);
      ml->m_scalarNew[m_FsVertFrame]->copyTo(*Fs_vert_frame[lev]);

      ml = ml->getFinerLevel();

    }

    // Another diagnostic - average vertical solute flux across the whole domain
    //    Real domainSize = m_domainWidth*m_domainHeight;
    Real volume = 0.0;
    Real domainAverageFs = computeSum(volume, horizAvFs, refRat, lev0Dx, Interval(0,0), 0);
    //    domainAverageFs = domainAverageFs / domainSize;
    domainAverageFs = domainAverageFs/volume;
    m_diagnostics.addDiagnostic(Diagnostics::m_averageVerticalSaltFlux, m_time, domainAverageFs);

    Real L2FsVertDiffusion = computeNorm(Fs_vert_diffusion, refRat, lev0Dx, Interval(0,0), 2);
    Real L2FsVertFluid = computeNorm(Fs_vert_fluid, refRat, lev0Dx, Interval(0,0), 2);
    Real L2FsVertFrame = computeNorm(Fs_vert_frame, refRat, lev0Dx, Interval(0,0), 2);

    m_diagnostics.addDiagnostic(Diagnostics::m_L2FsVertDiffusion, m_time, L2FsVertDiffusion);
    m_diagnostics.addDiagnostic(Diagnostics::m_L2FsVertFluid, m_time, L2FsVertFluid);
    m_diagnostics.addDiagnostic(Diagnostics::m_L2FsVertFrame, m_time, L2FsVertFrame);


    Real L1FsVertDiffusion = computeNorm(Fs_vert_diffusion, refRat, lev0Dx, Interval(0,0), 1);
    Real L1FsVertFluid = computeNorm(Fs_vert_fluid, refRat, lev0Dx, Interval(0,0), 1);
    Real L1FsVertFrame = computeNorm(Fs_vert_frame, refRat, lev0Dx, Interval(0,0), 1);

    m_diagnostics.addDiagnostic(Diagnostics::m_L1FsVertDiffusion, m_time, L1FsVertDiffusion);
    m_diagnostics.addDiagnostic(Diagnostics::m_L1FsVertFluid, m_time, L1FsVertFluid);
    m_diagnostics.addDiagnostic(Diagnostics::m_L1FsVertFrame, m_time, L1FsVertFrame);

    Real L0FsVertDiffusion = computeNorm(Fs_vert_diffusion, refRat, lev0Dx, Interval(0,0), 0);
    Real L0FsVertFluid = computeNorm(Fs_vert_fluid, refRat, lev0Dx, Interval(0,0), 0);
    Real L0FsVertFrame = computeNorm(Fs_vert_frame, refRat, lev0Dx, Interval(0,0), 0);

    m_diagnostics.addDiagnostic(Diagnostics::m_L0FsVertDiffusion, m_time, L0FsVertDiffusion);
    m_diagnostics.addDiagnostic(Diagnostics::m_L0FsVertFluid, m_time, L0FsVertFluid);
    m_diagnostics.addDiagnostic(Diagnostics::m_L0FsVertFrame, m_time, L0FsVertFrame);

    //Clean up
    for(int lev = 0; lev < nLevels; lev++)
    {
      delete horizAvFs[lev];
      horizAvFs[lev] = NULL;

      delete Fs_vert_diffusion[lev];
      Fs_vert_diffusion[lev] = NULL;

      delete Fs_vert_fluid[lev];
      Fs_vert_fluid[lev]= NULL;

      delete Fs_vert_frame[lev];
      Fs_vert_frame[lev]= NULL;
    }

    // Calculate sums over valid regions

    // Volume averaged sum i.e. sum[H*dx^(SpaceDim)]
    AMRSaltSum_new = ::computeSum(S_new, refRat, lev0Dx, Interval(0,0), 0);
    AMREnthalpySum_new = ::computeSum(H_new, refRat, lev0Dx, Interval(0,0), 0);

    //    Real Fs_bottom = ::computeSum(fluxBottom, refRat, lev0Dx, Interval(1, 1), 0)/scale;
    //    Real Fs_top = ::computeSum(fluxTop, refRat, lev0Dx, Interval(1, 1), 0)/scale;

    Real Fs_bottom = m_saltDomainFluxRegister.getFluxHierarchy(1, Side::Lo, 1.0);
    Real Fs_top = m_saltDomainFluxRegister.getFluxHierarchy(1, Side::Hi, 1.0);
    //    Real Fs_left = m_saltDomainFluxRegister.getFluxHierarchy(0, Side::Lo, 1.0);
    //    Real Fs_right = m_saltDomainFluxRegister.getFluxHierarchy(0, Side::Hi, 1.0);

    Real Fh_bottom = m_heatDomainFluxRegister.getFluxHierarchy(1, Side::Lo, 1.0);
    Real Fh_top =m_heatDomainFluxRegister.getFluxHierarchy(1, Side::Hi, 1.0);
    //    Real Fh_left = m_heatDomainFluxRegister.getFluxHierarchy(0, Side::Lo, 1.0);
    //    Real Fh_right = m_heatDomainFluxRegister.getFluxHierarchy(0, Side::Hi, 1.0);


    Real scale = 1/(m_dt*m_domainWidth);
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

    //    Real salt_mismatch = (ds_diff- flux_diff)/Fs_top;
    //    Real heat_mismatch = (H_diff - heat_flux_diff)/Fh_top;

    Real salt_mismatch = (ds_diff + salt_flux_diff);
    Real heat_mismatch = (H_diff + heat_flux_diff);

    Real rel_salt_mismatch = salt_mismatch/salt_flux_diff;
    Real rel_heat_mismatch = heat_mismatch/heat_flux_diff;


    m_diagnostics.addDiagnostic(Diagnostics::m_soluteFluxBottom, m_time, totalF_bottom);
    m_diagnostics.addDiagnostic(Diagnostics::m_soluteFluxTop, m_time, totalF_top);

    m_diagnostics.addDiagnostic(Diagnostics::m_heatFluxAbsMismatch, m_time, heat_mismatch);
    m_diagnostics.addDiagnostic(Diagnostics::m_saltFluxAbsMismatch, m_time, salt_mismatch);

    m_diagnostics.addDiagnostic(Diagnostics::m_heatFluxRelMismatch, m_time, rel_heat_mismatch);
    m_diagnostics.addDiagnostic(Diagnostics::m_saltFluxRelMismatch, m_time, rel_salt_mismatch);

    m_diagnostics.addDiagnostic(Diagnostics::m_heatFluxBottom, m_time, totalFh_bottom);
    m_diagnostics.addDiagnostic(Diagnostics::m_heatFluxTop, m_time, totalFh_top);

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

  } // end if level 0


  if (calcDiagnostics)
  {
    Real maxVel = getMaxVelocity();
    m_diagnostics.addDiagnostic(Diagnostics::m_maxVel, m_time, maxVel);

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

    Real sumLambda = ::computeNorm(amrLambda,nRef,m_dx, Interval(0,0), 2, m_level) ;
    Real maxLambda = ::computeMax(amrLambda, nRef, Interval(0,0), m_level); //::computeSum(amrLambda,nRef,m_dx, Interval(0,0), m_level);

    m_diagnostics.addDiagnostic(m_diagnostics.m_maxLambda, m_time, maxLambda);
    m_diagnostics.addDiagnostic(m_diagnostics.m_sumLambda, m_time, sumLambda);

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
  if (calcDiagnostics)
  {
    //    pout() << "AMRLevelMushyLayer::computeDiagnostics - compute average salinity" << endl;
    Vector<Real> averagedSalinity;
    horizontallyAverage(averagedSalinity, *m_scalarNew[m_liquidConcentration]);

    // Also do average of liquid salinity over liquid regions
    Real averageSl = averageOverLiquidRegion(m_liquidConcentration);

    //    if (printDiagnostics)
    //    {
    //      diagnosticsFile << "Horizontally averaged salinity at (0, 20%, 40%, 60%):  ";
    //      for (int i = 1; i < 5; i++)
    //      {
    //        int y_i = round(averagedSalinity.size()*0.2*(i-1));
    //        diagnosticsFile << averagedSalinity[y_i] << ", ";
    //      }


    int y_0 = round(averagedSalinity.size()*0.2*(1-1));
    int y_1 = round(averagedSalinity.size()*0.2*(2-1));
    int y_2 = round(averagedSalinity.size()*0.2*(3-1));
    int y_3 = round(averagedSalinity.size()*0.2*(4-1));

    m_diagnostics.addDiagnostic(Diagnostics::m_HorizAvSalinity0, m_time, averagedSalinity[y_0]);
    m_diagnostics.addDiagnostic(Diagnostics::m_HorizAvSalinity20, m_time, averagedSalinity[y_1]);
    m_diagnostics.addDiagnostic(Diagnostics::m_HorizAvSalinity40, m_time, averagedSalinity[y_2]);
    m_diagnostics.addDiagnostic(Diagnostics::m_HorizAvSalinity60, m_time, averagedSalinity[y_3]);

    m_diagnostics.addDiagnostic(Diagnostics::m_avSalinity, m_time, averageSl);


    //      diagnosticsFile << endl;

    //      diagnosticsFile << "Average salinity over all liquid regions: " << averageSl << endl;
    //    }



  }

  // Work out mushy layer depth
  if (calcDiagnostics)
  {
    Vector<Real> averagedPorosity;
    horizontallyAverage(averagedPorosity, *m_scalarNew[m_porosity]);



    //    if (printDiagnostics)
    //    {
    //      diagnosticsFile << "Mushy layer depth = ";
    int depth_i = 0;
    Real depth = -1.0;
    for (int i = averagedPorosity.size()-1; i >= 0 ; i--)
    {
      depth_i++;
      if (averagedPorosity[i] > 0.999)
      {
        depth = depth_i*m_dx;
        break;
      }

    }

    m_diagnostics.addDiagnostic(Diagnostics::m_mushDepth, m_time, depth);


    //      diagnosticsFile << endl;
    //    }

  }



  // Now lets work out some chimney geometry:
  // - how big
  // - what spacing
  bool doChimneyDiagnostics = false; // turn this off for now as it seems to be very slow
  if (m_level == 0 && doChimneyDiagnostics)
  {
    //    AMRLevelMushyLayer* ml = this;

    computeChimneyDiagnostics();
  }// end if level = 0


  // Deferred this to convergedToSteadyState(), so we also print out d/dt info
  //  bool printDiagnostics = (m_level == 0 && procID() ==0 ); // only print results on proc 0
  //  if (printDiagnostics)
  //  {
  //    m_diagnostics.printDiagnostics(m_time);
  //  }

}

void AMRLevelMushyLayer::getTotalFlux(LevelData<FluxBox>& totalFlux)
{
  int numComp = totalFlux.nComp();
  IntVect fluxGhost = totalFlux.ghostVect();

  LevelData<FluxBox> diffusiveTSlFlux(m_grids, numComp, fluxGhost);

  // Get diffusive salt and heat fluxes
  RefCountedPtr<AMRNonLinearMultiCompOp> HCOp = RefCountedPtr<AMRNonLinearMultiCompOp>(
      (AMRNonLinearMultiCompOp*)this->HCOpFact->AMRnewOp(m_problem_domain));


  LevelData<FArrayBox> HC(m_grids, numComp, fluxGhost + IntVect::Unit); // Need one more ghost cell than the flux
  fillHC(HC, m_time);

  // diffusive flux
  for (DataIterator dit = m_scalarNew[m_enthalpy]->dataIterator(); dit.ok(); ++dit)
  {
    HCOp->getFlux(diffusiveTSlFlux[dit], HC, diffusiveTSlFlux[dit].box(), dit(), 1.0);

    // Copy the vertical component to an farray box
    //    (*m_scalarNew[m_FsVertDiffusion])[dit].copy(diffusiveTSlFlux[dit], 1, 0, 1);

    EdgeToCell(diffusiveTSlFlux[dit], 1, (*m_scalarNew[m_FsVertDiffusion])[dit], 0, 1);
    (*m_scalarNew[m_FsVertDiffusion])[dit].mult(-1);

    EdgeToCell(diffusiveTSlFlux[dit], 1, (*m_vectorNew[m_FsDiffusion])[dit], 0, 0);
    EdgeToCell(diffusiveTSlFlux[dit], 1, (*m_vectorNew[m_FsDiffusion])[dit], 1, 1);


  }




  // Set F = F_{fluid} + F_{frame}
  computeTotalAdvectiveFluxes(totalFlux);



  for (DataIterator dit = m_scalarNew[m_bulkConcentration]->dataIterator(); dit.ok(); ++dit)
  {
    for (int idir=0; idir<SpaceDim; idir++)
    {
      Box b = totalFlux[dit][idir].box();

      totalFlux[dit][idir].minus(diffusiveTSlFlux[dit][idir], 0, 0, numComp);


      EdgeToCell(totalFlux[dit], 1, (*m_vectorNew[m_Fs])[dit], idir, idir);


      //      if (m_parameters.rayleighComposition != 0.0 || m_parameters.rayleighTemp != 0.0)
      //      {
      //        totalFlux[dit][idir].minus(fluidAdvFlux[dit][idir], b, 0, 0, numComp);
      //      }
      //
      //      if (m_parameters.nonDimVel != 0.0)
      //      {
      //        totalFlux[dit][idir].minus(frameAdvFlux[dit][idir], b, 0, 0, numComp);
      //      }
    }
  }

  totalFlux.exchange();

}

Real AMRLevelMushyLayer::averageOverLiquidRegion(int a_var)
{
  Real average = 0;
  Real vol = 0;

  // Need to be careful to do this right in parallel

  Box domBox = m_problem_domain.domainBox();

  Vector<Real> allAveraged(numProc(), 0.0);

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    Box box = m_grids[dit];

    box &= domBox;
    for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if ((*m_scalarNew[m_porosity])[dit](iv) > 0.999)
      {
        average += (*m_scalarNew[a_var])[dit](iv);
        vol += 1;
      }
    }
  }

  average /= vol;

  //  pout() << "HorizontallyAverage - got local averages, now broadcast/gather" << endl;

  // Broadcast/gather to compute averages over whole domain
  int srcProc = 0;

  gather(allAveraged, average, srcProc);

  Real globalAverage = 0;

  if (procID() == srcProc)
  {
    for (int ivec = 0; ivec<numProc(); ivec++)
    {
      globalAverage += allAveraged[ivec]/numProc();
    }


  }

  broadcast(globalAverage, srcProc);

  return globalAverage;
}

void AMRLevelMushyLayer::computeChimneyDiagnostics()
{
  IntVectSet liquidIVs, mushyIVs;
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    for (BoxIterator bit(m_grids[dit]); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      Real porosity = (*m_scalarNew[m_porosity])[dit](iv);

      // we use >= rather == because due to rounding errors we can get porosities very slightly above 1.0
      if (porosity >= 1.0)
      {
        liquidIVs |= iv;
      }
      else
      {
        mushyIVs |= iv;
      }
    }
  }

  mushyIVs.grow(1);

  Box solidBox = mushyIVs.minBox();
  int lowestSolidRegion  =solidBox.smallEnd(1);
  Box domainBox = m_problem_domain.domainBox();

  IntVect topCorner;
  for (int dir = 0; dir<SpaceDim; dir++)
  {
    if (dir == SpaceDim -1)
    {
      // Vertical direction
      topCorner += BASISV(dir)*lowestSolidRegion;
    }
    else
    {
      topCorner += BASISV(dir)*domainBox.bigEnd(dir);
    }
  }

  Box purelyLiquidBox(domainBox.smallEnd(), topCorner);

  // get rid off all IVs in an entirely liquid region, below the mushy layer
  liquidIVs -= purelyLiquidBox;
  //    purelyLiquidBox -= solidBox;

  // Start iterating in strips from the top of the domain.
  // Any liquid regions are either added to an existing IVS (if contiguous) or start their own

  //    Vector<IntVectSet*> channels;
  Vector<Channel*> channels;

  Box domBox = m_problem_domain.domainBox();
  Box b = adjCellHi(m_problem_domain.domainBox(), 1, 1);
  b.shift(1, -1);


  while (domBox.contains(b))
  {
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (liquidIVs.contains(iv))
      {
        // Two options: 1) new channel, 2) existing channel
        // Check if existing channel:
        bool existingChannel = false;
        for (int channel_i = 0; channel_i < channels.size(); channel_i++)
        {
          Channel& thisChannel = *channels[channel_i];

          if (thisChannel.borders(iv))
          {
            // Only add the IntVect to the channel if it isn't finished
            // i.e. if it hasn't reached the liquid yet.
            if (!thisChannel.isFinished())
            {
              thisChannel |= iv;
              //                pout() << iv << endl;
            }
            existingChannel = true;
            continue;
          }
        }

        // If we couldn't add IV to existing channel, make a new channel
        // a new channel must border a mushy region!
        bool bordersMush = mushyIVs.contains(iv);
        if (!existingChannel && bordersMush)
        {
          //            pout() << "new channel - " << iv << endl;
          Channel *newChannel = new Channel(iv);
          channels.push_back(newChannel);
        }
      } // end if intvect is liquid
    } // end loop over this box

    // Check if we've entered the actual liquid - have the channel widths massively increased?
    for (int channel_i = 0; channel_i < channels.size(); channel_i++)
    {
      Channel& thisChannel = *channels[channel_i];

      if (thisChannel.isFinished())
      {
        continue;
      }

      // Compare the width at the lowest point and second lowest point
      Real bottomWidth = thisChannel.width(Side::Lo, 0, m_dx);
      Real penultimateBottomWidth = thisChannel.width(Side::Lo, 1, m_dx);

      if (penultimateBottomWidth > 0.0
          //            && bottomWidth/a_dx > 10 // require a certain number of cells, as small channels can indeed grow quickly
          && bottomWidth > 3*penultimateBottomWidth)
      {
        // Want to end this channel at this point (and stop adding IVS), but can't do this until we create a channel class
        thisChannel.removeBottomCells();
        thisChannel.setFinished();
      }
    }

    b.shift(1, -1);
  } // end loop over boxes

  // Only keep channels which are higher than wide?
  Vector<Channel*> actualChannels;

  for (int i = 0; i < channels.size(); i++)
  {
    Channel& thisChannel = *channels[i];
    if (thisChannel.averageWidth(m_dx) < thisChannel.height(m_dx))
    {
      actualChannels.push_back(channels[i]);
    }

  }

  // Finally - may want to merge channels if we had wierd branching going on


  // Now we have the channels, can do statistics on them
  Real averageChannelWidth = 0.0, averageChannelSpacing = 0.0;

  Vector<Real> channelSpacing;
  Channel::channelSpacing(channelSpacing, actualChannels, m_dx, m_problem_domain);
  pout() << "Channel spacing: " << channelSpacing << endl;
  for (int i = 0; i < channelSpacing.size(); i++)
  {
    averageChannelSpacing += channelSpacing[i];
  }
  averageChannelSpacing /= channelSpacing.size();

  for (int i = 0; i < actualChannels.size(); i++)
  {

    Channel& thisChannel = *actualChannels[i];

    Real averageWidth = thisChannel.averageWidth(m_dx);
    pout() << "Channel average width = " << averageWidth << endl;

    averageChannelWidth += averageWidth;

  }
  averageChannelWidth /= actualChannels.size();


  if (m_level == 0)
  {
    m_diagnostics.addDiagnostic(Diagnostics::m_chimneySpacing, m_time, averageChannelSpacing);
    m_diagnostics.addDiagnostic(Diagnostics::m_chimneyWidth, m_time, averageChannelWidth);
  }

}

/*******/

Real AMRLevelMushyLayer::computeDt()
{

  if (m_timestepFailed)
  {
    return computeDt(false); // false means don't grow dt
  }
  else
  {
    return computeDt(true);
  }
}
Real AMRLevelMushyLayer::computeDt(bool growdt)
{
  // If we have a fixed dt  and no subcycling, we must make sure we use the dt on all levels
  if (m_fixedDt > 0 && !m_useSubcycling)
  {
    m_dt = m_fixedDt;
    return m_dt;
  }

  // If we've just reduced our timestep, keep it.
  if (m_timestepReduced)
  {
    m_timestepReduced = false; // Make sure we go back to normal next time
    m_dt = m_dt/m_dtReduction;
  }

  // Shouldn't be growing this dt - amr handles this itself
  Real grownDt = m_dt;
  if (growdt)
  {
    grownDt = m_dt*m_max_dt_growth;
  }

  Real solnDt = computeDt(m_cfl);

  Real newDt = min(grownDt, solnDt);

  return newDt;
}

Real AMRLevelMushyLayer::getMaxAdvVel()
{
  Real maxAdvULocal = 0.0;

  if (SpaceDim==3)
  {

    LevelData<FArrayBox> U(m_grids, SpaceDim);
    EdgeToCell(m_advVel, U);
    maxAdvULocal = ::computeNorm(U, NULL, 1, m_dx, Interval(0,SpaceDim-1), 0);
    return maxAdvULocal;
  }

  // alternate method
  Box domBox = m_problem_domain.domainBox();



  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    for (int dir=0; dir < SpaceDim-1; dir++)
    {
      Box faceBox = domBox.surroundingNodes(dir);

      FArrayBox& velDir = m_advVel[dit][dir];
      Box b = velDir.box();
      b.grow(-m_advVel.ghostVect());
      if (SpaceDim < 3)
      {
        //TODO - work out how to make this work in 3d
        b &= faceBox;
      }
      Real thisMax = velDir.norm(b, 0, 0);

      maxAdvULocal = max(maxAdvULocal, thisMax);
    }
  }

  // Need to gather max values across all processors
#ifdef CH_MPI
  Real recv;
  int result = MPI_Allreduce(&maxAdvULocal, &recv, 1, MPI_CH_REAL,
                             MPI_MAX, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("Sorry, but I had a communication error in getMaxAdvVel");
  }

  maxAdvULocal = recv;
#endif

  if (maxAdvULocal > 1e100)
  {
    pout() << "  WARNING: max advection velocity (level " << m_level << ") = " << maxAdvULocal << endl;
  }

  return maxAdvULocal;
}

Real AMRLevelMushyLayer::getMaxVelocity()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRlevelMushyLayer::getMaxVelocity" << endl;
  }

  Real maxAdvU = getMaxAdvVel();

  if (s_verbosity >= 3)
  {
    pout() << "AMRlevelMushyLayer::getMaxVelocity - max (face centered U) = " << maxAdvU << endl;
  }

  if (maxAdvU > 1e100)
  {
    maxAdvU = ::computeNorm(*m_vectorNew[m_fluidVel], NULL, 1, m_dx, Interval(0,SpaceDim-1), 0);
    if (s_verbosity >= 3)
    {
      pout() << "AMRlevelMushyLayer::getMaxVelocity - max (cell centered U) = " << maxAdvU << endl;
    }

    if (maxAdvU > 1e100)
    {
      maxAdvU = 0;
    }
  }

  if (s_verbosity > 3)
  {
    pout() << "AMRLevelMushyLayer::getMaxVelocity() - Max velocity = " << maxAdvU << endl;
  }
  maxAdvU = abs(maxAdvU);

  maxAdvU = max(maxAdvU, abs(m_parameters.nonDimVel));

  return maxAdvU;
}

Real AMRLevelMushyLayer::computeDt(Real cfl)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRlevelMushyLayer::computeDt" << endl;
  }

  if (m_grids.size() == 0)
  {
    return 0;
  }

  ParmParse ppMain("main");
  Real max_dt = -1;
  ppMain.query("max_dt", max_dt);
  Real maxAdvU = getMaxVelocity();


  // If we're doing advection with u/chi as the advection velocity,
  // then we need to use it for our cfl condition
  bool considerUChi = (m_advectionMethod == m_porosityInAdvection ||
  m_advectionMethod == m_porosityOutsideAdvection ||  m_advectionMethod == m_noPorosity);
  ppMain.query("consider_u_chi_dt", considerUChi);
  Real maxUChi = 0.0;

  if (considerUChi)
  {
    LevelData<FArrayBox> U_chi(m_grids, SpaceDim, IntVect::Zero);
    fillVectorField(U_chi, m_time, m_U_porosity, true);

    maxUChi = ::computeMax(U_chi, NULL, -1, Interval(0,SpaceDim-1));
    if (maxUChi > 1e10)
    {
      maxUChi = 1e-100;
    }
  }

  if (s_verbosity > 4)
  {
    pout() << "  Max(U) = " << maxAdvU << ", max(U/chi) = " << maxUChi << endl;
  }

  if (maxUChi < 1000*maxAdvU)
  {
    maxAdvU = max(maxAdvU, maxUChi);
  }

  Real newDT = cfl * m_dx  / maxAdvU;

  AMRLevelMushyLayer* ml = this;
  while (ml->getFinerLevel())
  {
    ml = ml->getFinerLevel();
  }

  Real finest_dx = ml->m_dx;

//  Real maxTemp = ::computeMax(*m_scalarNew[m_temperature], NULL, -1, Interval(0,0));
  Real maxPorosity = ::computeMax(*m_scalarNew[m_porosity], NULL, -1, Interval(0,0));
  //  Real maxUChi = ::computeMax(*m_vectorNew[m_U_porosity], NULL, -1, Interval(0,SpaceDim-1));
  //   Real minPerm = ::computeMin(*m_scalarNew[m_permeability], NULL, -1, Interval(0,0));

  Real buoyancy_acceleration = abs(max(m_parameters.m_buoyancyTCoeff, m_parameters.m_buoyancySCoeff) * maxPorosity);
  Real darcy_acceleration = abs(m_parameters.m_darcyCoeff*maxUChi); //*maxPorosity/minPerm
  Real viscous_acceleration = abs(m_parameters.m_viscosityCoeff*maxAdvU/(m_domainHeight));
//  Real acceleration = max(buoyancy_acceleration, darcy_acceleration);
//    acceleration = max(acceleration, viscous_acceleration);

   Real acceleration = computeNorm(*m_vectorNew[m_advectionSrc], NULL, 1, m_dx, Interval(0, SpaceDim-1), 0);

   // Ignore bogus values
   if (abs(acceleration) > 1e100)
   {
     acceleration = 0.0;
   }

    // Could just compute norm of advection velocity solve source term?

//  acceleration = max(acceleration);

  // Factor of 10 is fairly arbitrary
   Real accelCFL = cfl;
   ppMain.query("accelCFL", accelCFL);
  Real accelDt = sqrt(accelCFL*finest_dx/acceleration);

  bool printAccelDt = false;
  ppMain.query("printAccelDt", printAccelDt);
  if (printAccelDt && s_verbosity >= 2)
  {
    pout() << "  Max dt computed from acceleration = " << accelDt << endl;
    pout() << "  Accleration terms: buoyancy = " << buoyancy_acceleration << ", viscous = " << viscous_acceleration << ", darcy = " << darcy_acceleration << endl;

    //    Real accelCFL = accelDt*m_dt/m_dx;
    //    if (accelCFL > 1.0)
    //      {
    //     pout() << "  WARNING - CFL computed due to advection terms  = " << accelCFL << endl;
    //      }
  }

  // make sure we consider acceleration for determinig dt at time 0.
  // The initialisation may have generated a small initial velocity, with a corresponding CFL timestep,
  // which we want to override.
  bool useAccelDt = false;
  ppMain.query("useAccelDt", useAccelDt);

  bool useInitAccelDt = true;
  ppMain.query("useInitAccelDt", useInitAccelDt);

  if (maxAdvU == 0 || (m_time == 0 && useInitAccelDt) || useAccelDt)
  {
    newDT = min(newDT, accelDt);
  }

  if (max_dt > 0)
  {
    newDT = min(max_dt, newDT);
  }

  return newDT;
}

/*******/
Real AMRLevelMushyLayer::computeInitialDt()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRlevelMushyLayer::computeInitialDt" << endl;
  }

  Real dt = computeDt(m_initial_dt_multiplier);

  Real max_dt = 1e100;
  ParmParse ppMain("main");
  ppMain.query("max_dt", max_dt);
  Real max_init_dt = max_dt*m_initial_dt_multiplier;

  ppMain.query("max_init_dt", max_init_dt);

  dt = min(dt, max_init_dt);

  return dt;
}

void AMRLevelMushyLayer::setFluxRegistersZero()
{
  if (s_verbosity > 5)
  {
    pout() << "AMRLevelMushyLayer::setFluxRegistersZero on level " << m_level << endl;
  }

  if (m_grids.size() == 0)
  {
    return;
  }

  for (int a_scalarVar = 0; a_scalarVar < m_numScalarVars; a_scalarVar++)
  {
    if (m_makeFluxRegForScalarVar[a_scalarVar])
    {
      if (m_fluxRegisters[a_scalarVar] == NULL)
      {
        return;
      }
      m_fluxRegisters[a_scalarVar]->setToZero();
    }
  }

  for (int a_vectorVar = 0; a_vectorVar < m_numVectorVars; a_vectorVar++)
  {
    if (m_makeFluxRegForVectorVar[a_vectorVar])
    {
      m_vectorFluxRegisters[a_vectorVar]->setToZero();
    }
  }

  m_fluxRegHC->setToZero();

  m_heatDomainFluxRegister.setToZero();
  m_saltDomainFluxRegister.setToZero();

}
/*******/



/*******/

AMRLevelMushyLayer*
AMRLevelMushyLayer::getCoarsestLevel()
{
  AMRLevelMushyLayer* amrML = this;
  while(amrML->getCoarserLevel())
  {
    amrML = amrML->getCoarserLevel();
  }
  return amrML;
}


AMRLevelMushyLayer*
AMRLevelMushyLayer::getCoarserLevel() const {
  AMRLevelMushyLayer* amrADCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
  {
    amrADCoarserPtr =  dynamic_cast<AMRLevelMushyLayer*>(m_coarser_level_ptr);

    if (amrADCoarserPtr == NULL)
    {
      MayDay::Error("AMRLevelMushyLayer::getCoarserLevel: dynamic cast failed");
    }
  }

  return amrADCoarserPtr;
}

/*******/
AMRLevelMushyLayer*
AMRLevelMushyLayer::getFinerLevel() const {
  AMRLevelMushyLayer* amrADFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
  {
    amrADFinerPtr = dynamic_cast<AMRLevelMushyLayer*>(m_finer_level_ptr);

    if (amrADFinerPtr == NULL)
    {
      MayDay::Error(
          "AMRLevelMushyLayer::getFinerLevel: dynamic cast failed");
    }
  }

  return amrADFinerPtr;
}

int AMRLevelMushyLayer::getFinestLevel()
{
  //  int finest_level = m_level;
  // index through levels to find out what finest level is
  AMRLevelMushyLayer* thisMLPtr = this;
  while (!thisMLPtr->finestLevel())
  {
    thisMLPtr = thisMLPtr->getFinerLevel();
  }
  CH_assert(thisMLPtr->finestLevel());
  return thisMLPtr->m_level;
}



void AMRLevelMushyLayer::fillVectorField(LevelData<FArrayBox>& a_vector,
                                         Real a_time, int a_var, bool doInterior, bool quadInterp)
{
  Interval vectorComps(0, SpaceDim - 1);

  Real old_time = m_time - m_dt;

  LevelData<FArrayBox> porosityOld(m_grids, 1, a_vector.ghostVect());
  LevelData<FArrayBox> porosityNew(m_grids, 1, a_vector.ghostVect());

  fillScalars(porosityOld, old_time, m_porosity, true);
  fillScalars(porosityNew, m_time, m_porosity, true);

  if (a_var == m_U_porosity)
  {

    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      (*m_vectorOld[m_U_porosity])[dit].copy((*m_vectorOld[m_fluidVel])[dit]);
      (*m_vectorNew[m_U_porosity])[dit].copy((*m_vectorNew[m_fluidVel])[dit]);

      for (int dir=0; dir<SpaceDim; dir++)
      {
        (*m_vectorOld[m_U_porosity])[dit].divide(porosityOld[dit], 0, dir, 1);
        (*m_vectorNew[m_U_porosity])[dit].divide(porosityNew[dit], 0, dir, 1);
      }
    }

//    m_vectorOld[m_U_porosity]->exchange();
//    m_vectorNew[m_U_porosity]->exchange();
  }


  if (doInterior)
  {
    if (abs(a_time - old_time) < TIME_EPS)
    {
      m_vectorOld[a_var]->copyTo(vectorComps, a_vector, vectorComps);
//      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
//      {
//        a_vector[dit].copy((*m_vectorOld[a_var])[dit]);
//      }
    }
    else if (abs(a_time - m_time) < TIME_EPS)
    {
            m_vectorNew[a_var]->copyTo(vectorComps, a_vector, vectorComps);
//      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
//      {
//        a_vector[dit].copy((*m_vectorNew[a_var])[dit]);
//      }
    } else {
      // do linear interpolation in time
      timeInterp(a_vector, a_time, *m_vectorOld[a_var], old_time,
                 *m_vectorNew[a_var], m_time, vectorComps);

    }
  }

  int nGhost = a_vector.ghostVect()[0];
  if (nGhost > 0)
  {

    // if necessary, do interpolation from coarser levels
    if (m_level > 0)
    {
      const DisjointBoxLayout& levelGrids = a_vector.getBoxes();
      const DisjointBoxLayout& thisLevelsGrids = m_grids;

      const IntVect& growVect = a_vector.ghostVect();
      int velGrow = growVect[0];

      // if grids for a_vel are the same as those for this level, and
      // there are no ghost cells, we don't need to do this at all
      if (!((velGrow == 0) && (levelGrids == thisLevelsGrids)))
      {
        // fill in coarse-fine BC data by conservative linear interp
        AMRLevelMushyLayer& crseLevel = *getCoarserLevel();

        LevelData<FArrayBox>& oldCrseVector = *crseLevel.m_vectorOld[a_var];
        LevelData<FArrayBox>& newCrseVector = *crseLevel.m_vectorNew[a_var];
        const DisjointBoxLayout& crseGrids = oldCrseVector.getBoxes();
        const ProblemDomain& crseDomain = crseLevel.problemDomain();
        int nRefCrse = crseLevel.refRatio();

        Real crse_new_time = crseLevel.m_time;
        Real crse_dt = crseLevel.dt();
        Real crse_old_time = crse_new_time - crse_dt;
        Real crse_time_interp_coeff;

        // check for "essentially 0 or 1"
        if (abs(a_time - crse_old_time) < TIME_EPS)
        {
          crse_time_interp_coeff = 0.0;
        } else if (abs(a_time - crse_new_time) < TIME_EPS)
        {
          crse_time_interp_coeff = 1.0;
        } else {
          crse_time_interp_coeff = (a_time - crse_old_time) / crse_dt;
        }

        //        bool doSecondOrderCorners = (CFinterpOrder_advection == 2);

        PiecewiseLinearFillPatch filpatcher(levelGrids, crseGrids, SpaceDim,
                                            crseDomain, nRefCrse, velGrow,
                                            false); //doSecondOrderCorners);

        filpatcher.fillInterp(a_vector, oldCrseVector, newCrseVector,
                              crse_time_interp_coeff, 0, 0, SpaceDim);

        if (quadInterp)
        {
          LevelData<FArrayBox> avCrseVect(crseGrids, SpaceDim, oldCrseVector.ghostVect());
          ::timeInterp(avCrseVect, a_time, oldCrseVector, crse_old_time, newCrseVector, crse_new_time, Interval(0,SpaceDim-1));

          m_quadCFInterpVector.coarseFineInterp(a_vector, avCrseVect);
        }

      }
    }

    // need to set physical boundary conditions here
    const DisjointBoxLayout& levelGrids = a_vector.getBoxes();

    // Changing U and U/chi to other BCs, as these ones were filling interior CF ghost cells
    if (1==0)
    {
      if (m_isViscous)
      {

        // This BC fills all ghost cells, which is what we want
        BCHolder bc = m_physBCPtr->extrapolationFuncBC();
        DataIterator dit = a_vector.dataIterator();
        for (dit.reset(); dit.ok(); ++dit)
        {
          bc(a_vector[dit], levelGrids[dit], m_problem_domain, m_dx, false); // not homogeneous
        }
      } else {
        VelBCHolder velBC(m_physBCPtr->tracingVelFuncBC());
        velBC.applyBCs(a_vector, levelGrids, m_problem_domain, m_dx, false); // inhomogeneous
      }
    }
    else if (a_var == m_Ustar || a_var==m_UpreProjection || a_var == m_advUstar || a_var == m_advUpreProjection
        || a_var == m_fluidVel or a_var == m_U_porosity)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        BCHolder viscousBC = m_physBCPtr->velFuncBC(idir, m_viscousBCs, Interval(idir, idir) );

        DataIterator dit = m_grids.dataIterator();

        for (dit.reset(); dit.ok(); ++dit)
        {
          viscousBC(a_vector[dit],
                    m_grids[dit],
                    m_problem_domain, m_dx,
                    false); // not homogeneous
        }
      }

    } else {
      // Don't know how to set BCs for this variable, don't bother
    }

  }

  a_vector.exchange();
}




void AMRLevelMushyLayer::smoothScalarField(LevelData<FArrayBox>& a_phi, int a_var, Real a_smoothing)
{
  // Just solve on one level!

  int finest_level = m_level;
  // index through levels to find out what finest level is
  AMRLevelMushyLayer* thisMLPtr = this;
  while (!thisMLPtr->finestLevel())
  {
    thisMLPtr = thisMLPtr->getFinerLevel();
  }
  CH_assert(thisMLPtr->finestLevel());
  finest_level = thisMLPtr->m_level;

  // amr grid info for solvers
  Vector<DisjointBoxLayout> AmrGrids(finest_level + 1);
  Vector<int> AmrRefRatios(finest_level + 1);
  Vector<Real> AmrDx(finest_level + 1);
  ProblemDomain baseDomain;

  // loop over levels, set up for AMRMultiGrid solve
  thisMLPtr = this;
  int startLev = m_level;
  // if crser level exists, define it as well for BC's
  if (startLev > 0)
  {
    startLev = startLev - 1;
    thisMLPtr = thisMLPtr->getCoarserLevel();
  }
  AMRLevelMushyLayer* startLevelPtr = thisMLPtr;

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


  int numLevels = finest_level + 1;

  // This is a Helmholtz operator
  Real alpha = 1.0;
  Real beta = -a_smoothing;

  BCHolder bc;
  // Get inhomogeneous form of BCs
  getScalarBCs(bc, a_var, false);


  RelaxSolver<LevelData<FArrayBox> > bottomSolver;
  bottomSolver.m_verbosity = s_verbosity;

  AMRMultiGrid<LevelData<FArrayBox> > diffusionSolver;

  AMRPoissonOpFactory diffusiveOpFactory;
  diffusiveOpFactory.define(baseDomain, AmrGrids,
                            AmrRefRatios, AmrDx[0],
                            bc,
                            alpha, beta);

  AMRLevelOpFactory<LevelData<FArrayBox> >& castFact =
      (AMRLevelOpFactory<LevelData<FArrayBox> >&) diffusiveOpFactory;

  diffusionSolver.define(baseDomain, castFact,
                         &bottomSolver, numLevels);


  diffusionSolver.m_verbosity = m_verbosity_multigrid;
  diffusionSolver.m_eps = s_viscous_solver_tol;
  diffusionSolver.m_normThresh = s_viscous_solver_tol;
  diffusionSolver.m_hang = s_viscous_solver_tol;

  // Just solve on this level
  Vector<LevelData<FArrayBox>*> correction(finest_level + 1,
                                           NULL);
  Vector<LevelData<FArrayBox>*> rhs(finest_level + 1,
                                    NULL);
  //
  thisMLPtr = startLevelPtr;
  // Only solve on one level!
  for (int lev = startLev; lev <= startLev; lev++)
  {
    const DisjointBoxLayout& levelGrids = thisMLPtr->m_grids;
    // recall that AMRMultiGrid can only do one component.
    // rhs has no ghost cells
    rhs[lev] = &(a_phi);

    //soln has one layer of ghost cells
    IntVect ghostVect(D_DECL(1, 1, 1));
    correction[lev] = new LevelData<FArrayBox>(levelGrids,
                                               1, ghostVect);
    // initialize corr to 0
    DataIterator levelDit = correction[lev]->dataIterator();

    LevelData<FArrayBox>& levelCorr =
        *(correction[lev]);
    setValLevel(levelCorr, 0.0);
    thisMLPtr = thisMLPtr->getFinerLevel();
  }

  diffusionSolver.solve(correction, rhs, m_level, m_level, true, false);

  // Apply smoothing correction

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
  {
    a_phi[dit].copy((*correction[startLev])[dit]);
  }

}

void AMRLevelMushyLayer::fillScalarFace(LevelData<FluxBox>& a_scal,
                                        Real a_time, const int a_var, bool doInterior , bool quadInterp)
{
  CH_TIME("AMRLevelMushyLayer::fillScalarFace");
  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::fillScalarFace - field: " << m_scalarVarNames[a_var] << ", level: " << m_level << endl;
  }

  CellToEdgeAveragingMethod method = arithmeticAveraging; //Arithmetic mean
  // Geometric mean for permeability was making pressure projection impossible
  // || a_var == m_pressureScaleVar
  if (a_var == m_porosity || a_var == m_permeability) // || a_var == m_permeability
  {
    method = geometricAveraging; // Geometric mean
    //    method = arithmeticAveraging;
  }

  fillScalarFace(a_scal, a_time, a_var, method, doInterior, quadInterp);

}

void AMRLevelMushyLayer::fillScalarFace(LevelData<FluxBox>& a_scal, Real a_time,
                                        const int a_var, CellToEdgeAveragingMethod method, bool doInterior,
                                        bool quadInterp, Real smoothing)
{
  CH_TIME("AMRLevelMushyLayer::fillScalarFace");

  // Need a ghost vector here to get domain faces correct
  LevelData<FArrayBox> temp(m_grids, a_scal.nComp(), IntVect::Unit);
  fillScalars(temp, a_time, a_var, doInterior, quadInterp); //NB - this includes exchanges (and corner copiers)

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
    bc.applyBCs(a_scal, a_scal.disjointBoxLayout(), m_problem_domain, m_dx, false);
  }
  else if (a_var == m_permeability)
  {
    EdgeVelBCHolder bc(m_physBCPtr->permeabilityFaceBC());
    bc.applyBCs(a_scal, a_scal.disjointBoxLayout(), m_problem_domain, m_dx, false);
  }

  // This is important
  a_scal.exchange();

  doRegularisationOps(a_scal, a_var);

  // Try a corner copier?
  if (a_scal.ghostVect()[0] > 0)
  {
    CornerCopier cornerCopy(a_scal.disjointBoxLayout(), a_scal.disjointBoxLayout(),
                            m_problem_domain, a_scal.ghostVect(), true);
    a_scal.exchange(a_scal.interval(), cornerCopy);
  }

}
// Fill a single component of a scalar field
void AMRLevelMushyLayer::fillScalars(LevelData<FArrayBox>& a_scal, Real a_time,
                                     const int a_var, bool doInterior, bool quadInterp)
{

  CH_TIME("AMRLevelMushyLayer::fillScalars");
  if (s_verbosity >= 5)
  {
    pout() << "AMRLevelMushyLayer::fillScalars - field: " << m_scalarVarNames[a_var] << ", level: " << m_level << endl;
  }

  const DisjointBoxLayout& levelGrids = m_grids;

  Interval scalComps = Interval(0,0); // only works with this for now
  CH_assert(a_scal.nComp() == 1);

  Real old_time = m_time - m_dt;

  if (doInterior)
  {
    CH_TIME("AMRLevelMushyLayer::fillScalars::interior");
    if (s_verbosity >= 6)
    {
      pout() << "  AMRLevelMushyLayer::fillScalars - start interior filling, a_time:" << a_time << ", old_time: " << old_time << endl;
    }

    if (abs(a_time - old_time) < TIME_EPS)
    {
      m_scalarOld[a_var]->copyTo(scalComps, a_scal, scalComps);
    } else if (abs(a_time - m_time) < TIME_EPS)
    {
      m_scalarNew[a_var]->copyTo(scalComps, a_scal, scalComps);
    } else {
      // do linear interpolation in time
      timeInterp(a_scal, a_time, *m_scalarOld[a_var], old_time,
                 *m_scalarNew[a_var], m_time, scalComps);
    }

  }

  if (s_verbosity >= 6)
  {
    pout() << "  AMRLevelMushyLayer::fillScalars - done interior filling" << endl;
  }

  // Only do CF interp if we have ghost vectors
  if (m_level > 0  && a_scal.ghostVect() >= IntVect::Unit )
  {
    CH_TIME("AMRLevelMushyLayer::fillScalars::coarseFineInterp");

    // Make sure fill objects are defined
    if (!m_piecewiseLinearFillPatchScalarOne.isDefined())
    {
      defineCFInterp();
    }

    // fill in coarse-fine BC data by conservative linear interp
    AMRLevelMushyLayer& crseLevel = *getCoarserLevel();

    LevelData<FArrayBox>& oldCrseScal = *(crseLevel.m_scalarOld[a_var]);
    LevelData<FArrayBox>& newCrseScal = *(crseLevel.m_scalarNew[a_var]);

    const DisjointBoxLayout& crseGrids = oldCrseScal.getBoxes();
    const ProblemDomain& crseDomain = crseLevel.problemDomain();
    int nRefCrse = crseLevel.refRatio();

    Real crse_new_time = crseLevel.m_time;
    Real crse_dt = crseLevel.dt();
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

    if (s_verbosity >= 5)
    {
      pout() << "  AMRLevelMushyLayer::fillScalars - crse_time_inter_coeff = " << crse_time_interp_coeff << endl;
    }

    const IntVect& scalGrowVect = a_scal.ghostVect();
    int scalGrow = scalGrowVect[0];

    {
      CH_TIME("AMRLevelMushyLayer::fillScalars::linearFillPatch");

      if (scalGrow == 1)
      {
        m_piecewiseLinearFillPatchScalarOne.fillInterp(a_scal, oldCrseScal, newCrseScal,
                                                       crse_time_interp_coeff,
                                                       scalComps.begin(), scalComps.end(), scalComps.size());
      }
      else if (scalGrow == 2)
      {
        m_piecewiseLinearFillPatchScalarTwo.fillInterp(a_scal, oldCrseScal, newCrseScal,
                                                       crse_time_interp_coeff,
                                                       scalComps.begin(), scalComps.end(), scalComps.size());
      }
      else if (scalGrow == 3)
      {
        m_piecewiseLinearFillPatchScalarThree.fillInterp(a_scal, oldCrseScal, newCrseScal,
                                                         crse_time_interp_coeff,
                                                         scalComps.begin(), scalComps.end(), scalComps.size());
      }
      else if (scalGrow == 4)
      {
        m_piecewiseLinearFillPatchScalarFour.fillInterp(a_scal, oldCrseScal, newCrseScal,
                                                        crse_time_interp_coeff,
                                                        scalComps.begin(), scalComps.end(), scalComps.size());
      }
      else
      {
        pout() << "ERROR: No fill patch object for " << scalGrow << "ghost cells" << endl;
        MayDay::Error("No fill patch object for this number of ghost cells");
      }

    }

    if (quadInterp && m_quadCFInterpScalar.isDefined())
    {
      CH_TIME("AMRLevelMushyLayer::fillScalars::quadCFinterp");
      if (s_verbosity >= 6)
      {
        pout() << "  AMRLevelMushyLayer::fillScalars - do quad interp" << endl;
      }

      LevelData<FArrayBox> avCrseScal(crseGrids, 1, oldCrseScal.ghostVect());
      ::timeInterp(avCrseScal, a_time, oldCrseScal, crse_old_time, newCrseScal, crse_new_time, scalComps);

      if (s_verbosity >= 6)
      {
        pout() << "  AMRLevelMushyLayer::fillScalars - made quadInterp CF BC" << endl;
      }

      // Fill BCs on coarse scalar - doesn't seem to make a difference
      //      crseLevel.fillScalars(avCrseScal, a_time, a_var, false, true);

      m_quadCFInterpScalar.coarseFineInterp(a_scal, avCrseScal);
    }

  }

  const ProblemDomain& physDomain = problemDomain();

  int numGhost = a_scal.ghostVect()[0];
  if (s_verbosity >= 6)
  {
    pout() << "  AMRLevelMushyLayer::fillScalars - num ghost: " << numGhost << endl;
  }

  // Do domain BCs if we have ghost cells
  if (numGhost > 0)
  {
    CH_TIME("AMRLevelMushyLayer::fillScalars::domainBCs");

    if (s_verbosity >= 5)
    {
      pout() << "  AMRLevelMushyLayer::fillScalars - do  BCs" << endl;
    }

    BCHolder thisBC;
    getScalarBCs(thisBC, a_var, false); // inhomogeneous

    // loop over boxes
    DataIterator dit = a_scal.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      thisBC(a_scal[dit], m_grids[dit], physDomain, m_dx, false); // inhomogeneous
    }

  }

  if (s_verbosity >= 5)
  {
    pout() << "  AMRLevelMushyLayer::fillScalars - regularisation ops" << endl;
  }

  ParmParse ppMain("main");

  doRegularisationOps(a_scal, a_var);

  {
    CH_TIME("AMRLevelMushyLayer::fillScalars::exchange");

    a_scal.exchange();

    // Try a corner copier?
    bool doCorners = true;
    ppMain.query("scalarExchangeCorners", doCorners);
    if (a_scal.ghostVect()[0] > 0 && doCorners)
    {
      if (s_verbosity >= 5)
      {
        pout() << "  AMRLevelMushyLayer::fillScalars - corner copier" << endl;
      }

      CornerCopier cornerCopy(m_grids, m_grids, m_problem_domain, a_scal.ghostVect(), true);
      a_scal.exchange(a_scal.interval(), cornerCopy);
    }

  }

  if (s_verbosity >= 5)
  {
    pout() << "  AMRLevelMushyLayer::fillScalars - finished" << endl;
  }

}

void AMRLevelMushyLayer::doRegularisationOps(LevelData<FluxBox>& a_scal, int a_var)
{
  if (a_var == m_porosity || a_var == m_permeability || a_var == m_bulkConcentration)
  {
    DataIterator dit2 = a_scal.dataIterator();

    ParmParse ppMain("main");
      bool use_new_version = false;
      ppMain.query("fortranRegularisationFace", use_new_version);

    for (dit2.reset(); dit2.ok(); ++dit2)
    {
      Box b = a_scal[dit2].box();

      for (int dir=0; dir<SpaceDim; dir++)
      {
        if (use_new_version)
        {
          doRegularisationOpsNew(a_var, a_scal[dit2][dir]);
        }
        else
        {
          doRegularisationOps(a_var, a_scal[dit2][dir]);
        }
      }
    }
  }

}

void AMRLevelMushyLayer::doRegularisationOps(int a_var, FArrayBox& a_state)
{
  CH_TIME("AMRLevelMushyLayer::doRegularisationOpsOld");

  Box b = a_state.box();

  if (a_var == m_porosity)
  {
    // Ensure porosity is greater than 0 and less than or equal to 1
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      a_state(iv) = max(m_lowerPorosityLimit, a_state(iv));
      a_state(iv) = min(1.0, a_state(iv));
    }

  }
  else if (a_var == m_permeability)
  {
    //Real minPermeability = m_parameters.calculatePermeability(m_lowerPorosityLimit);
    Real minPermeability = pow(m_lowerPorosityLimit,3);
    // Ensure porosity is greater than 0
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();

      a_state(iv) = max(minPermeability, a_state(iv));
    }

  }
  else if (a_var == m_bulkConcentration)
  {
    // Was hoping the enthalpy-concentration update would guarantee this, but maybe not

    //Real minPermeability = m_parameters.calculatePermeability(m_lowerPorosityLimit);
    Real minVal = -m_parameters.compositionRatio;
    Real maxVal = 0;


    // Ensure porosity is greater than 0
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();

      a_state(iv) = max(minVal, a_state(iv));
      a_state(iv) = min(maxVal, a_state(iv));
    }


  }


}

void AMRLevelMushyLayer::doRegularisationOpsNew(int a_var, FArrayBox& a_state)
{
  CH_TIME("AMRLevelMushyLayer::doRegularisationOpsNew");
  Box region = a_state.box();

  if (a_var == m_porosity)
  {

    Real maxVal = 1.0;
    Real minVal = 0.01; //m_lowerPorosityLimit;

    FORT_SETMINMAXVAL( CHF_FRA(a_state),
                       CHF_BOX(region),
                       CHF_CONST_REAL(minVal),
                       CHF_CONST_REAL(maxVal));
  }
  else if (a_var == m_permeability)
  {
    //Real minPermeability = m_parameters.calculatePermeability(m_lowerPorosityLimit);
    Real minPermeability = pow(m_lowerPorosityLimit,3);

    FORT_SETMINVAL( CHF_FRA(a_state),
                          CHF_BOX(region),
                          CHF_CONST_REAL(minPermeability));

  }
  else if (a_var == m_bulkConcentration)
  {
    // Was hoping the enthalpy-concentration update would guarantee this, but maybe not

    Real minVal = -m_parameters.compositionRatio;
    Real maxVal = 0;

    FORT_SETMINMAXVAL( CHF_FRA(a_state),
                       CHF_BOX(region),
                       CHF_CONST_REAL(minVal),
                       CHF_CONST_REAL(maxVal));

  }


}

void AMRLevelMushyLayer::doRegularisationOps(LevelData<FArrayBox>& a_scal,
                                             int a_var)
{
  CH_TIME("AMRLevelMushyLayer::doRegularisationOps");
  DataIterator dit2 = a_scal.dataIterator();

  ParmParse ppMain("main");
  bool use_new_version = true;
  ppMain.query("fortranRegularisation", use_new_version);

  for (dit2.reset(); dit2.ok(); ++dit2)
  {

    if (use_new_version)
    {
      doRegularisationOpsNew(a_var, a_scal[dit2]);
    }
    else
    {
      doRegularisationOps(a_var, a_scal[dit2]);
    }

  }
}

// Different interface for backward compatability
int AMRLevelMushyLayer::convertBCType(const int a_implicitBC)
{
  int a_explicitBC=0;
  Real temp=0.0;
  convertBCType(a_implicitBC, temp, a_explicitBC, temp);
  return a_explicitBC;
}

void AMRLevelMushyLayer::convertBCType(const int a_implicitBC,  const Real a_implicitVal,
                                       int a_explicitBC, Real a_explicitVal)
{
  // Default BCs
  //  int ibcType = AdvectIBC::m_dirichlet;

  m_physBCPtr->convertBCType(a_implicitBC, a_implicitVal, a_explicitBC, a_explicitVal);

}



PhysIBC* AMRLevelMushyLayer::getScalarIBCs(int a_var)
{


  IntVect bcTypeHi, bcTypeLo;
  RealVect bcValHi, bcValLo;

  for (int dir=0; dir<SpaceDim; dir++)
  {
    bcTypeHi[dir] = convertBCType(m_parameters.bcTypeScalarHi[dir]);
    bcTypeLo[dir] = convertBCType(m_parameters.bcTypeScalarLo[dir]);
  }

  m_physBCPtr->applyFrameAdvectionBC(bcTypeHi, bcTypeLo);

  // Do lambda differently as I don't think it has a plume
  if (a_var == m_lambda || a_var == m_lambda_porosity)
  {

    Real plumeVal = 1.0;

    Vector<Real> plumeBounds;
    plumeBounds.resize(2);
    plumeBounds[0] = 0.0; plumeBounds[1] = 0.0;

    for (int dir=0; dir<SpaceDim; dir++)
    {
      bcValHi[dir] = 1.0;
      bcValLo[dir] = 1.0;

      if (a_var == m_lambda_porosity)
      {
        bcValHi[dir] /= m_parameters.bcValPorosityHi[dir];
        bcValLo[dir] /= m_parameters.bcValPorosityLo[dir];
      }

      bcTypeHi[dir] = AdvectIBC::m_dirichlet;
      bcTypeLo[dir] = AdvectIBC::m_dirichlet;
    }
    return m_physBCPtr->scalarTraceIBC(bcValHi, bcValLo,
                                       bcTypeLo, bcTypeHi,
                                       plumeVal, plumeBounds);

  }
  else
  {
    pout() << "WARNING, NO IBCs for scalar var " << a_var << endl;
    return  m_physBCPtr->scalarTraceIBC(bcValLo, bcValHi,
                                        bcTypeLo, bcTypeHi);
  }
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

void AMRLevelMushyLayer::getScalarBCs(BCHolder& thisBC, int a_var, bool a_homogeneous)
{

  computeInflowOutflowAdvVel();
  // Reset this
  m_physBCPtr->setAdvVel(&m_totalAdvVel);


  if (a_var == m_temperature)
  {
    thisBC = m_physBCPtr->BasicthetaFuncBC(a_homogeneous, &m_totalAdvVel);
  }
  else if (a_var == m_enthalpy)
  {
    thisBC = m_physBCPtr->BasicEnthalpyFuncBC(a_homogeneous, &m_totalAdvVel);
  }
  else if(a_var == m_enthalpySolidus)
  {
    //    thisBC = m_physBCPtr->BasicScalarFuncBC(m_parameters.HSolidusTop,
    //                                            m_parameters.HSolidusBottom,
    //                                            m_parameters.HSolidusPlume,
    //                                            a_homogeneous);
    thisBC = m_physBCPtr->BasicSolidusBC(a_homogeneous, &m_totalAdvVel);
  }
  else if(a_var == m_enthalpyEutectic)
  {
    //    thisBC = m_physBCPtr->BasicScalarFuncBC(m_parameters.HEutecticTop,
    //                                            m_parameters.HEutecticBottom,
    //                                            m_parameters.HEutecticPlume,
    //                                            a_homogeneous);
    thisBC = m_physBCPtr->BasicEutecticBC(a_homogeneous, &m_totalAdvVel);
  }
  else if(a_var == m_enthalpyLiquidus)
  {
    //    thisBC = m_physBCPtr->BasicScalarFuncBC(m_parameters.HLiquidusTop,
    //                                            m_parameters.HLiquidusBottom,
    //                                            m_parameters.HLiquidusPlume,
    //                                            a_homogeneous);
    thisBC = m_physBCPtr->BasicLiquidusBC(a_homogeneous, &m_totalAdvVel);
  }
  else if (a_var == m_porosity)
  {
    thisBC = m_physBCPtr->BasicPorosityFuncBC(a_homogeneous);
  }
  else if (a_var == m_pressure)
  {
    thisBC = m_physBCPtr->BasicPressureFuncBC(a_homogeneous);
  }
  else if (a_var == m_permeability)
  {
    thisBC = m_physBCPtr->BasicPermeabilityFuncBC(a_homogeneous);
  }
  else if (a_var == m_liquidConcentration)
  {
    thisBC = m_physBCPtr->ThetaLFuncBC(a_homogeneous, &m_totalAdvVel);
  }
  else if (a_var == m_bulkConcentration)
  {
    thisBC = m_physBCPtr->ThetaFuncBC(a_homogeneous, &m_totalAdvVel);
  }
  else if (a_var == m_lambda)
  {
    thisBC = m_physBCPtr->basicLambdaFuncBC(a_homogeneous);
  }
  else if (a_var == m_lambda_porosity)
  {
    thisBC = m_physBCPtr->basicLambdaFuncBC(a_homogeneous,
                                            true); // true - scale with porosity
  }
  else
  {
    pout() << "WARNING No BCs for " << m_scalarVarNames[a_var] << endl;
    MayDay::Warning("AMRLevelMushyLayer::getScalarBCs - No BCs for scalar variable, using extrapolation ");

    thisBC = m_physBCPtr->extrapFuncBC();
  }
}

void AMRLevelMushyLayer::stokesDarcyForcing(LevelData<FArrayBox>& T, Real time)
{
  Real timescale = 0.5;
  ParmParse ppParameters("parameters");
  ppParameters.query("forcing_timescale", timescale);

  if (timescale > 0)
  {

    for (DataIterator dit=T.dataIterator(); dit.ok(); ++dit)
    {
      Box b = T[dit].box();
      for (BoxIterator bit(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        RealVect loc;
        getLocation(iv, loc, m_dx);

        Real one = 1.0;
        Real zero = 0.0;


        Real sinPart = sin(M_PI*loc[0]);
        Real thisT = min(one,  (time/timescale)-sinPart);

        Real powPart = pow(loc[0]-2.5,2);
        thisT = exp(-powPart/(time/timescale) ) ;
        //			thisT = log(1+time/timescale)*exp(-pow(loc[0]-2.5,2)) ;

        //			thisT = 0.01;

        thisT = max(zero, thisT);
        thisT = min(one, thisT);
        T[dit](iv) = thisT;
      }
    }
  }

}

bool AMRLevelMushyLayer::crashed()
{

  Vector<int> vars;
  vars.append(m_enthalpy);
  vars.append(m_vorticity);
  vars.append(m_enthalpySrc);
  vars.append(m_lambda);

  for (int i=1; i<vars.size(); i++)
  {
    Real max = computeMax(*m_scalarNew[vars[i]], NULL, m_ref_ratio, Interval(0,0));
    Real min = computeMin(*m_scalarNew[vars[i]], NULL, m_ref_ratio, Interval(0,0));

    Real limit = 1e50;
    if (max > limit || min < -limit)
    {
      return true;
    }
  }


  return false;
}

bool AMRLevelMushyLayer::currentCFLIsSafe(bool printWarning)
{

  Real maxAdvU = getMaxVelocity();
  Real newCFL = maxAdvU * m_dt / m_dx;

  // CFL above 1 definitely unstable. CFL > 2*user limit means the velocity has increased rapidly
  if (newCFL > 1.0 || newCFL > m_cfl*2)
  {
    if (printWarning)
    {
      pout() << "WARNING: new max(U) = " << maxAdvU << " means that CFL = " << newCFL << " which may be unstable" << endl;
//      pout() << "Therefore, we will be skipping fluid advection during this timestep" << endl;
    }
//    return false;
    //always return true for now
    return true;
  }

  return true;
}


/*******/

#include "NamespaceFooter.H"


