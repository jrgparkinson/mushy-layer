
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

#include "SetValLevel.H"
#include "BackwardEuler.H"
#include "VCAMRPoissonOp2.H"
#include "Divergence.H"
#include "ExtrapFillPatch.H"
#include "PhysIBC.H"

#include "EnthalpyVariablesF_F.H"
#include "utils_F.H"

#include "AMRLevelMushyLayer.H"
#include "Channel.h"
#include "MushyLayerUtils.H"


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

  Real Tnorm = convergedToSteadyState(ScalarVars::m_enthalpy);
  Real Cnorm = convergedToSteadyState(ScalarVars::m_bulkConcentration);
  Real Unorm =  convergedToSteadyState(m_fluidVel, true);

  if (m_opt.computeDiagnostics)
  {
    // If a diagnostic period has been declared, check this time has passed since we last produced diagnostics
      if (m_opt.diagnostics_period > 0 &&
          m_time - m_prev_diag_output < m_opt.diagnostics_period)
      {
        // do nothing
      }
      else
      {

        m_diagnostics.addDiagnostic(DiagnosticNames::diag_dTdt, m_time, Tnorm);
        m_diagnostics.addDiagnostic(DiagnosticNames::diag_dSdt, m_time, Cnorm);
        m_diagnostics.addDiagnostic(DiagnosticNames::diag_dUdt, m_time, Unorm);

        // Can print diagnostics now if on processor 0
        bool printDiagnostics = (m_level == 0 && procID() ==0 ); // only print results on proc 0
        if (printDiagnostics)
        {
          m_diagnostics.addDiagnostic(DiagnosticNames::diag_dt, m_time, m_dt);
          m_diagnostics.printDiagnostics(m_time);
        }

        m_prev_diag_output = m_time;
      }
  }

  // Stop worrying about this for now
  bool Tstalled = false;
  bool Cstalled = false;
  bool Ustalled = false;

  bool metricConverged = false;
  // For some simulations, steady state can be defined by some other metric

  // Only test for this on level 0
  if (m_opt.computeDiagnostics && !getCoarserLevel())
  {
    if ( (m_parameters.physicalProblem == PhysicalProblems::m_HRL
        || m_parameters.physicalProblem == PhysicalProblems::m_rayleighBenard))
    {
      metricConverged = metricConverged || m_diagnostics.movingAverageHasConverged(DiagnosticNames::diag_Nu, m_time, m_dt);
    }

  }

  // Two ways we can 'converge':
  // 1) All fields have converged
  // 2) All fields have stalled (d/dt = constant) but metrics have converged
  bool hasConverged = false;

  bool velConverged = m_opt.ignoreVelocitySteadyState || Unorm < m_opt.steadyStateCondition;
  bool velStalled = m_opt.ignoreVelocitySteadyState || Ustalled;

  bool Tconverged = Tnorm < m_opt.steadyStateCondition;
  bool Cconverged = Cnorm < m_opt.steadyStateCondition;

  // Override to let us converge without salinity having actually converged
  if (m_opt.ignoreBulkConcSteadyState)
  {
    Cconverged = true;
  }

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
//  Real maxVel = ::computeNorm(*m_vectorNew[VectorVars::m_advectionVel], NULL, -1, m_dx, Interval(0, SpaceDim-1), 0);
//  if (hasConverged && m_doAutomaticRestart && abs(maxVel) < 1e-3)
//  {
//    pout() << "Max Vel = " << maxVel << ". Trying to restart with a small perturbation to kick off instability" << endl;
//    addPerturbation(ScalarVars::m_enthalpy, 1e-3);
//    hasConverged = false;
//    m_doAutomaticRestart = false;
//  }

  if (hasConverged && m_dt < 1e-10)
  {
    pout() << "AMRLevelMushyLayer::convergedToSteadyState - all fields converged but dt < 1e-10 so keep solving" << endl;
    return false;
  }


  if (hasConverged && m_time < m_opt.min_time)
  {
    pout() << "AMRLevelMushyLayer::convergedToSteadyState - converged but time < min_time (" << m_time << " < " << m_opt.min_time << ")" << endl;
    return false;
  }


  return hasConverged;
}

bool AMRLevelMushyLayer::solvingFullDarcyBrinkman()
{
  return m_parameters.isDarcyBrinkman();
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
    diff.copyTo(*m_vectorNew[VectorVars::m_dUdt]);
  }
  else if (!vector && a_var == ScalarVars::m_enthalpy)
  {
    diff.copyTo(*m_scalarNew[ScalarVars::m_dHdt]);
  }
  else if(!vector && a_var == ScalarVars::m_bulkConcentration)
  {
    //    horizontallyAverage(*m_scalarNew[ScalarVars::m_dSdt], diff);
    diff.copyTo(*m_scalarNew[ScalarVars::m_dSdt]);
  }

}

Real AMRLevelMushyLayer::convergedToSteadyState(const int a_var, bool vector)
{

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

  DataIterator dit(m_grids);

  if (m_opt.spongeHeight > 0)
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
          topCorner += round(m_numCells[1]*m_opt.spongeHeight) * BASISV(dir);
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

  Real norm = ::computeNorm(diff, NULL, -1, m_dx, Interval(0, largestDim), m_opt.steadyStateNormType);

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




//void AMRLevelMushyLayer::horizontallySmooth(LevelData<FluxBox>& a_flux)
//{
//  CH_TIME("AMRLevelMushyLayer::horizontallySmooth");
//
//  DataIterator dit = a_flux.dataIterator();
//
//  int dir = 0; // horizontally averaging
//
//
//
//  for (dit.reset(); dit.ok(); ++dit)
//  {
//    for (int velDir=0; velDir < SpaceDim; velDir++)
//    {
//      FArrayBox& vel = a_flux[dit][velDir];
//      Box b = vel.box();
//
//      for (BoxIterator bit(b); bit.ok(); ++bit)
//      {
//        IntVect iv = bit();
//        IntVect ivUp = iv+BASISV(dir);
//
//        if (b.contains(ivUp))
//        {
//
//          Real neighbour = vel(ivUp);
//
//          // Make sure we don't change the velocity if the neighboroung value is much larger
//          // this will prevent changing zero velocity cells, or including NaN type values
//          if( abs(neighbour) < 1e100 ) //abs(vel(iv)) > 1e-15 &&
//          {
//            vel(iv) = (1-m_opt.postTraceSmoothing)*vel(iv)+m_opt.postTraceSmoothing*neighbour;
//          }
//
//        }
//      }
//    }
//
//  }
//}

bool AMRLevelMushyLayer::doVelocityAdvection()
{
  if (m_parameters.darcy==0)
  {
    return false;
  }

  if (m_opt.doEulerPart)
  {
    return true;
  }
  else
  {
    return false;
  }

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




/*******/
Real AMRLevelMushyLayer::advance()
{

  ///////////////////////////////////////////////////////////////
  // This method advances the solution by one timestep on a level.
  // Before solving the equations, we need to set some things up...
  ///////////////////////////////////////////////////////////////

  // If a timestep fails, we sometimes try reducing the timestep by setting m_dtReduction to some factor > 0
  // If we've just reduced the timestep, set this < 0 to ensure we don't do it again
  if (m_dtReduction > 0)
  {
    m_dtReduction = -1;
  }

  // Get the coarser level, so we can work out later if this is in fact the coarsest level
  AMRLevelMushyLayer *amrMLcrse = NULL;
  if (m_level > 0)
  {
    amrMLcrse = getCoarserLevel();
  }

  // Do some setup operations on only the coarsest level
  // If the coarser level pointer is null, it means that this is the coarsest level
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
      if (m_time < m_opt.perturbationTime && (m_time + m_dt) > m_opt.perturbationTime)
      {
        addPerturbation(ScalarVars::m_enthalpy, m_opt.delayedPerturbation, m_opt.perturbationWavenumber, m_opt.perturbationPhaseShift);
      }

      amrMLptr = amrMLptr->getFinerLevel();
    }
  }


  if (s_verbosity >= 1)
  {
    // Let's determine the CFL number we're running at
    Real maxAdvU = getMaxVelocityForCFL();
    m_computedCFL = m_dt*maxAdvU/m_dx;

    pout() << " AMRLevelMushyLayer::advance (level = "<< m_level << ", time=" << m_time << ", dt = " << m_dt << ", CFL=" << m_computedCFL << ")" << endl;
  }

  // Reset BCs in case they change with time
  m_parameters.setTime(m_time); // BCs are stored in m_parameters
  this->m_physBCPtr->Time(m_time);
  setAdvectionBCs(); // Reset BCs on advection physics objects

  // Elliptic operators get the BC from m_parameters when they're defined (later)
  if (m_opt.rampBuoyancy > 0)
  {
    if (m_parameters.rayleighComposition == 0)
    {
      m_parameters.rayleighComposition = m_opt.initRaC;
    }
    if (m_parameters.rayleighTemp == 0)
    {
      m_parameters.rayleighTemp = m_opt.initRaT;
    }

    m_parameters.rayleighComposition = m_parameters.rayleighComposition *m_opt.rampBuoyancy;
    m_parameters.rayleighTemp = m_parameters.rayleighTemp *m_opt.rampBuoyancy;

    m_parameters.rayleighComposition = min(m_opt.maxRaC, m_parameters.rayleighComposition);
    m_parameters.rayleighTemp = min(m_opt.maxRaT, m_parameters.rayleighTemp);

    pout() << "RaC = " << m_parameters.rayleighComposition  << ", RaT = " <<  m_parameters.rayleighTemp  << endl;
  }

  // Compute temperature, porosity, liquid concentration and sold concentration
  // from the enthalpy and bulk concentration to ensure they're consistent with each other
  updateEnthalpyVariables();

#ifdef CH_MPI
  {
    CH_TIME("MPI_Barrier exchange");
    MPI_Barrier(Chombo_MPI::comm);
  }
#endif

  // Increment time by dt, so m_time represent the time at the end of this timestep,
  // and m_time-m_dt is the time at the start of this timestep

//  Real old_time = m_time;
//  Real half_time = new_time - m_dt / 2;
  Real new_time = m_time + m_dt;
  m_time = new_time;


  // Compute d(porosity)/dt from the old timestep in case we need it at some point
  // First create data structures
  LevelData<FArrayBox> oldPorosity(m_dPorosity_dt.disjointBoxLayout(), 1, m_dPorosity_dt.ghostVect());
  LevelData<FArrayBox> newPorosity(m_dPorosity_dt.disjointBoxLayout(), 1, m_dPorosity_dt.ghostVect());
  // Now fill them with data, making sure we fill ghost cells
  fillScalars(oldPorosity, m_time-m_dt, m_porosity, true, true);
  fillScalars(newPorosity, m_time, m_porosity, true, true);

  for (DataIterator dit = m_dPorosity_dt.dataIterator(); dit.ok(); ++dit)
  {
    m_dPorosity_dt[dit].copy(newPorosity[dit]);
    m_dPorosity_dt[dit].minus(oldPorosity[dit]);
    m_dPorosity_dt[dit].divide(m_dt);
  }

  // The 'new' variables have been calculated at the end of the previous timestep,
  // i.e. at the start of this timestep. So move them from new->old.
  // do this after we've incremented m_time for consistency with coarser levels
  copyNewToOldStates();

  //Gradually ramp up temperature forcing
//  if (m_parameters.physicalProblem == PhysicalProblems::m_poiseuilleFlow)
//  {
//    stokesDarcyForcing(*m_scalarNew[ScalarVars::m_temperature], new_time);
//    stokesDarcyForcing(*m_scalarOld[ScalarVars::m_temperature], old_time);
//  }

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

//  if (solvingFullDarcyBrinkman())
  if (doVelocityAdvection())
  {
    // This fills all the ghost cells of advectionSourceTerm
    computeAdvectionVelSourceTerm(advectionSourceTerm);

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

  if (!this->isCurrentCFLSafe(true))
  {
    doAdvectiveSrc = false;
  }

  // Another sanity check
  //  Divergence::levelDivergenceMAC(*m_scalarNew[ScalarVars::m_divUadv], m_advVel, m_dx);
  //  Real  maxDivU = ::computeNorm(*m_scalarNew[ScalarVars::m_divUadv], NULL, 1, m_dx, Interval(0,0));
  //  pout() << "  Sanity check: max(div u) = " << maxDivU << endl;

  // always* advect lambda (and update flux registers)
  // do this as soon as we have advection velocities, in case we want to
  // correct them prior to HC advection
  // * skip if we're danger of violating the CFL condition
  if (doAdvectiveSrc)
  {
    advectLambda(true);
  }

  if (m_opt.includeTracers)
  {
    this->computeRadianceIntensity();
    this->advectPassiveTracer();
    this->advectActiveTracer();
  }

//  if (m_newLevel && m_level > 0)
//  {
//    if (m_opt.skipNewLevelScalars)
//    {
//      pout() << "First time step on a new level - skipping advection and diffusion" << endl;
//      m_opt.doScalarAdvectionDiffusion = false;
//    }
//    m_newLevel = false;
//  }

  if (!(m_parameters.physicalProblem == PhysicalProblems::m_poiseuilleFlow ||
      m_parameters.physicalProblem == PhysicalProblems::m_soluteFluxTest ||
      m_parameters.physicalProblem == PhysicalProblems::m_zeroPorosityTest) && m_opt.doScalarAdvectionDiffusion)
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

    addHeatSource(srcMultiComp);

    bool doFRUpdates = true;

    exitStatus = multiCompAdvectDiffuse(HC_old, HC_new, srcMultiComp, doFRUpdates, doAdvectiveSrc);

    bool solverFailed = (exitStatus == 2 || exitStatus == 4 || exitStatus == 6);

    // Get back the answer if solver was a success
    HC_new.copyTo(Interval(0,0), *m_scalarNew[ScalarVars::m_enthalpy], Interval(0,0));
    HC_new.copyTo(Interval(1,1), *m_scalarNew[ScalarVars::m_bulkConcentration], Interval(0,0));

    updateEnthalpyVariables();

    if (solverFailed)
    {
      if (m_opt.ignoreSolveFails)
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


  // compute cell centered velocities
  // for problems where the momentum equation has time dependence

//  if (solvingFullDarcyBrinkman())
  if (isVelocityTimeDependent())
  {
    // If we're skipping advective srcs for this timestep, skip this too
    if (!doAdvectiveSrc)
    {
      DataIterator dit = m_vectorNew[VectorVars::m_UdelU]->dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
      {
        (*m_vectorNew[VectorVars::m_UdelU])[dit].setVal(0.0);
      }
    }

    bool doFRupdates = true;
//    bool compute_uDelU = !m_opt.implicitAdvectionSolve && doAdvectiveSrc;
    bool compute_uDelU = doAdvectiveSrc;
    bool doProjection = true;
    computeCCvelocity(advectionSourceTerm, m_time-m_dt, m_dt, doFRupdates, doProjection, compute_uDelU);
  }


  getExtraPlotFields();

  //supposed to return dt but the return value is never used, so don't.
  //if the return value is ever used we will know because this will break it
  return -1;
}


void AMRLevelMushyLayer::addHeatSource(LevelData<FArrayBox>& src)
{
  // This is where we add in a heat source, if required
      Real gaussian_heat_source_size = 0.0;
      Real gaussian_heat_source_width = 0.0;
      Real gaussian_heat_source_depth = 0.0;
      Real gaussian_heat_source_xpos = m_domainWidth/2;

      ParmParse ppHeatSource("heatSource");
      ppHeatSource.query("size", gaussian_heat_source_size);
      ppHeatSource.query("width", gaussian_heat_source_width);
      ppHeatSource.query("depth", gaussian_heat_source_depth);
  //    ppHeatSource.query("depth", gaussian_heat_source_depth);
      ppHeatSource.query("xpos", gaussian_heat_source_xpos);

      // enthalpy is in the first component of the source term
      int Hcomp = 0;

      if (gaussian_heat_source_size != 0.0)
      {
        for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          for (BoxIterator bit = BoxIterator(m_grids[dit]); bit.ok(); ++bit)
          {
            IntVect iv = bit();
            RealVect loc;
            ::getLocation(iv, loc, m_dx);

            src[dit](iv, Hcomp) = gaussian_heat_source_size/(gaussian_heat_source_width*sqrt(2*M_PI))
                * exp(-0.5*pow((loc[0]-gaussian_heat_source_xpos)/gaussian_heat_source_width, 2))
                * 0.5*(1 + tanh(10*(loc[1]-(m_domainHeight-gaussian_heat_source_depth) ) ));

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
  setVelZero(*m_scalarNew[VectorVars::m_fluidVel], a_limit);
}
void AMRLevelMushyLayer::setVelZero(LevelData<FArrayBox>& a_vel, Real a_limit, int a_radius)
{
  if (a_limit < 0)
  {
    a_limit = m_opt.solidPorosity;
  }

  for (DataIterator dit = a_vel.dataIterator(); dit.ok(); ++dit)
  {
    setVelZero(a_vel[dit], (*m_scalarNew[ScalarVars::m_porosity])[dit], a_limit, a_radius);
  }
}

void AMRLevelMushyLayer::setVelZero(LevelData<FluxBox>& a_vel, Real a_limit)
{
  if (a_limit < 0)
  {
    a_limit = m_opt.solidPorosity;
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


  m_advVel.exchange();

  // Do lambda advection
  m_scalarNew[ScalarVars::m_lambda]->copyTo(Interval(0,0), *m_scalarOld[ScalarVars::m_lambda], Interval(0,0));

  // Want to get the lambda flux back so we can remove it later
  LevelData<FluxBox> lambdaFlux(m_grids, 1);
  advectScalar(m_lambda, m_lambda, m_advVel, true, lambdaFlux); // advect without a diffusive source term

  setValLevel(*m_vectorNew[VectorVars::m_advVelCorr], 0.0);

}


void AMRLevelMushyLayer::updateEnthalpyVariables()
{


  CH_TIME("AMRLevelMushyLayer::updateEnthalpyVariables");

  // Apply BCs?
//  fillScalars(*m_scalarNew[ScalarVars::m_bulkConcentration], m_time, ScalarVars::m_bulkConcentration);
//  fillScalars(*m_scalarNew[ScalarVars::m_enthalpy], m_time, m_enthalpy);

  LevelData<FArrayBox> HC(m_grids, 2, IntVect::Unit);
  fillHC(HC, m_time);

//  ::updateEnthalpyVariables(*m_scalarNew[ScalarVars::m_enthalpy], *m_scalarNew[ScalarVars::m_bulkConcentration],
//                            *m_scalarNew[ScalarVars::m_temperature], *m_scalarNew[ScalarVars::m_liquidConcentration], *m_scalarNew[ScalarVars::m_solidConcentration],
//                            *m_scalarNew[ScalarVars::m_porosity],
//                            *m_scalarNew[ScalarVars::m_enthalpySolidus],*m_scalarNew[ScalarVars::m_enthalpyLiquidus],*m_scalarNew[ScalarVars::m_enthalpyEutectic],
//                            m_parameters);
  ::updateEnthalpyVariables(HC,
                            *m_scalarNew[ScalarVars::m_temperature], *m_scalarNew[ScalarVars::m_liquidConcentration], *m_scalarNew[ScalarVars::m_solidConcentration],
                            *m_scalarNew[ScalarVars::m_porosity],
                            *m_scalarNew[ScalarVars::m_enthalpySolidus],*m_scalarNew[ScalarVars::m_enthalpyLiquidus],*m_scalarNew[ScalarVars::m_enthalpyEutectic],
                            m_parameters);


  doRegularisationOps(*m_scalarNew[ScalarVars::m_liquidConcentration], m_liquidConcentration);
  doRegularisationOps(*m_scalarNew[ScalarVars::m_porosity], m_porosity);

  // A few alterations for test problems
  if (m_parameters.physicalProblem == PhysicalProblems::m_poiseuilleFlow)
  {
    initialDataPoiseuille();
  }
  else if(m_parameters.physicalProblem == PhysicalProblems::m_convectionMixedPorous ||
      m_parameters.physicalProblem == PhysicalProblems::m_zeroPorosityTest)
  {
    fillFixedPorosity(*m_scalarNew[ScalarVars::m_porosity]);
    m_scalarNew[ScalarVars::m_porosity]->copyTo(*m_scalarOld[ScalarVars::m_porosity]);
  }


  // stop doing this
//  computeLambdaPorosity();

}

void AMRLevelMushyLayer::computeLambdaPorosity()
{
  for (DataIterator dit = m_scalarNew[ScalarVars::m_lambda]->dataIterator(); dit.ok(); ++dit)
  {
    (*m_scalarNew[ScalarVars::m_lambda_porosity])[dit].copy((*m_scalarNew[ScalarVars::m_lambda])[dit]);
    (*m_scalarNew[ScalarVars::m_lambda_porosity])[dit].divide((*m_scalarNew[ScalarVars::m_porosity])[dit]);

    (*m_scalarOld[ScalarVars::m_lambda_porosity])[dit].copy((*m_scalarOld[ScalarVars::m_lambda])[dit]);
    (*m_scalarOld[ScalarVars::m_lambda_porosity])[dit].divide((*m_scalarOld[ScalarVars::m_porosity])[dit]);
  }
}


void AMRLevelMushyLayer::updateEnthalpyVariablesOld()
{
  ::updateEnthalpyVariables(*m_scalarOld[ScalarVars::m_enthalpy], *m_scalarOld[ScalarVars::m_bulkConcentration],
                            *m_scalarOld[ScalarVars::m_temperature], *m_scalarOld[ScalarVars::m_liquidConcentration], *m_scalarOld[ScalarVars::m_solidConcentration],
                            *m_scalarOld[ScalarVars::m_porosity],
                            *m_scalarOld[ScalarVars::m_enthalpySolidus],*m_scalarOld[ScalarVars::m_enthalpyLiquidus],*m_scalarOld[ScalarVars::m_enthalpyEutectic],
                            m_parameters);
}


int AMRLevelMushyLayer::getMaxLevel()
{
  return m_opt.max_possible_level;
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

        if (m_parameters.physicalProblem == PhysicalProblems::m_poiseuilleFlow)
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
    FArrayBox solidFraction((*m_scalarNew[ScalarVars::m_permeability])[dit].box(), 1);

    //First do new timestep
    solidFraction.setVal(1.0);
    solidFraction -= porosityNew[dit];
    ::calculatePermeability((*m_scalarNew[ScalarVars::m_permeability])[dit],
                            solidFraction, m_parameters, m_dx);

    // Now do old timestep
    solidFraction.setVal(1.0);
    solidFraction -= porosityOld[dit];
    ::calculatePermeability((*m_scalarOld[ScalarVars::m_permeability])[dit],
                            solidFraction, m_parameters, m_dx);
  }

  calculateCoarseFineBoundaries(ScalarVars::m_permeability);
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
      (AMRNonLinearMultiCompOp*) m_HCOpFact->AMRnewOp(m_problem_domain));

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

void AMRLevelMushyLayer::horizAverage()
{
  horizontallyAverage(*m_scalarNew[ScalarVars::m_enthalpy], *m_scalarNew[ScalarVars::m_enthalpy]);
  horizontallyAverage(*m_scalarNew[ScalarVars::m_bulkConcentration], *m_scalarNew[ScalarVars::m_bulkConcentration]);

  copyNewToOldStates();
  updateEnthalpyVariables();

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
  Box domBox = m_problem_domain.domainBox();

  IntVect smallEnd =domBox.smallEnd();
  int y_init = smallEnd[1];

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

        IntVect ivFluxBox;

        // This is a bit of a hack, but basically works.

        if (fluxDir.box().contains(ivUp))
        {

          ivFluxBox = ivUp;
        }
        else
        {
          ivFluxBox = iv;
        }

        //averaged[y_i-y_init] += fluxDir(iv, comp)*m_dx/m_opt.domainWidth;
        averaged[y_i-y_init] += fluxDir(ivFluxBox, comp)*m_dx/m_opt.domainWidth;

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
  Real width = m_opt.domainWidth;
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

  for (DataIterator dit = m_scalarNew[ScalarVars::m_bulkConcentration]->dataIterator(); dit.ok(); ++dit)
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
  restartVars.push_back(ScalarVars::m_enthalpy);
  restartVars.push_back(ScalarVars::m_bulkConcentration);
  restartVars.push_back(ScalarVars::m_pressure);
  restartVars.push_back(ScalarVars::m_lambda);

  Vector<int> vectRestartVars;
  vectRestartVars.push_back(VectorVars::m_fluidVel);

  // Make sure we copy ghost cells
  for (DataIterator dit = m_scalarNew[0]->dataIterator(); dit.ok(); ++dit)
  {

    for (int i = 0; i < restartVars.size(); i++)
    {
      int var = restartVars[i];
      if (! (var == ScalarVars::m_pressure && ignorePressure))
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

    computeScalarAdvectiveSrcHC(full_src, totalAdvectiveFlux, doFRupdates);

    if (m_opt.skipSaltUpdate)
    {
      for (dit.reset(); dit.ok(); ++dit)
      {
        full_src[dit].setVal(0.0, 1);
      }
    }

  } // end if compute advective src

  full_src.copyTo(Interval(0,0), *m_scalarNew[ScalarVars::m_enthalpySrc], Interval(0,0));
  full_src.copyTo(Interval(1,1), *m_scalarNew[ScalarVars::m_saltEqnSrcGodunov], Interval(0,0));

  if (m_opt.skipHCUpdate)
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

  if (m_opt.timeIntegrationOrder == 2)
  {
    //       MayDay::Error("multiCompAdvectDiffuse - TGA not implemented yet");
    //exitStatus = TGAUpdateScalar(a_var, a_src, converged);

    m_enthalpySalinityTGA->updateSoln(a_phi_new,
                                      a_phi_old, full_src, finerFRPtr, coarserFRPtr,
                                      coarserDataOldPtr, coarserDataNewPtr, old_time, tCoarserOld,
                                      tCoarserNew, m_dt, m_level, false); //false - don't zero phi

    baseLevBE = dynamic_cast<BaseLevelHeatSolver<LevelData<FArrayBox>, FluxBox, LevelFluxRegister> * > (&(*m_enthalpySalinityTGA));
  }
  else
  {

//    if (m_opt.noMultigrid)
//    {
//      a_phi_old.copyTo(a_phi_new);
//
//      LevelData<FArrayBox> thisSrc(m_grids, 2, IntVect::Unit);
//      full_src.copyTo(thisSrc);
//      for (DataIterator dit = thisSrc.dataIterator(); dit.ok(); ++dit)
//      {
//        thisSrc[dit].mult(m_dt);
//        thisSrc[dit].plus(a_phi_old[dit]);
//      }
//
//      for (int i=0; i < m_opt.noMultigridIter; i++)
//      {
//
//        // Need to define solvers first and foremost
//        //        defineSolvers(m_time);
//
//        // Try replacing with single level relax solver
//        RefCountedPtr<AMRNonLinearMultiCompOp> amrpop = RefCountedPtr<AMRNonLinearMultiCompOp>(
//            (AMRNonLinearMultiCompOp*) m_HCOpFact->AMRnewOp(m_problem_domain));
//
//        // Solving (1-dt Op) HC^{n+1} = HC^{n}
//        amrpop->setAlphaAndBeta(1.0, -m_dt);
//
//        LevelData<FArrayBox> a_res(m_grids, 2, IntVect::Unit);
//
//        amrpop->relax(a_phi_new, thisSrc, 1000);
//
//        amrpop->residual(a_res, a_phi_new, thisSrc, false);
//        Real maxRes = ::computeNorm(a_res, NULL, 1, m_dx, Interval(0, 1), 0);
//        pout() << "  Max residual = " << maxRes << endl;
//
//        if (maxRes < 1e-10)
//        {
//          pout() << "   Converged " << endl;
//          break;
//
//        }
//
//      }
//    }
//    else
//    {

      m_enthalpySalinityBE->updateSoln(a_phi_new,
                                       a_phi_old, full_src, finerFRPtr, coarserFRPtr,
                                       coarserDataOldPtr, coarserDataNewPtr, old_time, tCoarserOld,
                                       tCoarserNew, m_dt, m_level, false); //false - don't zero phi

      baseLevBE = dynamic_cast<BaseLevelHeatSolver<LevelData<FArrayBox>, FluxBox, LevelFluxRegister> * > (&(*m_enthalpySalinityBE));
//    }

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

    //           MAC Projection (level
    pout() << "  HC solve       (level " << m_level << "): exit status " << exitStatus << ", solver residual = " << residual << ", num MG iterations = " << num_iter << endl;

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

  IntVect ghostVect = IntVect::Unit;
  IntVect fluxGhostVect = edgeScalTotal.ghostVect();

  LevelData<FluxBox> edgeScalFluidAdv(m_grids, numComp, fluxGhostVect);
  LevelData<FluxBox> edgeScalFrameAdvection(m_grids, numComp, fluxGhostVect);

  // Need grown version of this
  LevelData<FArrayBox> TCl_old(m_grids, numComp, advect_grow);
  LevelData<FArrayBox> HC_old(m_grids, numComp, advect_grow);


  fillHC(HC_old, old_time,
         true,  // fill interior?
         (m_opt.CFinterpOrder_advection==2)); // do quadratic interpolation at CF boundaries?
  fillTCl(TCl_old, old_time,
          true,  // fill interior?
          (m_opt.CFinterpOrder_advection==2)); //  do quadratic interpolation at CF boundaries? - need this for corner cells that border CF and domain boundaries


  // Compute frame advection src term
  computeScalarAdvectiveFluxMultiComp(edgeScalFrameAdvection, m_frameAdvVel,
                                      m_patchGodHC, HC_old,
                                      old_time, m_dt);


  if (m_opt.reflux_enthalpy == m_opt.reflux_concentration && m_opt.allowMulticompAdvection)
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

    if (m_opt.reflux_enthalpy)
    {
      enthalpy_advVel = &advVel_reflux;
    }
    else
    {
      enthalpy_advVel = &advVel_noReflux;
    }

    if (m_opt.reflux_concentration)
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
    //    (*m_scalarNew[ScalarVars::m_FsVertFluid])[dit].copy(edgeScalFluidAdv[dit], 1, 0, 1);
    //    (*m_scalarNew[ScalarVars::m_FsVertFrame])[dit].copy(edgeScalFrameAdvection[dit], 1, 0, 1);
    EdgeToCell(edgeScalFluidAdv[dit], 1, (*m_scalarNew[ScalarVars::m_FsVertFluid])[dit], 0, 1);
    EdgeToCell(edgeScalFrameAdvection[dit], 1, (*m_scalarNew[ScalarVars::m_FsVertFrame])[dit], 0, 1);

    for (int dir=0; dir <SpaceDim; dir++)
    {
      EdgeToCell(edgeScalFluidAdv[dit], 1, (*m_vectorNew[VectorVars::m_FsFluid])[dit], dir, dir);
    }

    (*m_scalarNew[ScalarVars::m_FsVertFluid])[dit].mult(-1);
    (*m_scalarNew[ScalarVars::m_FsVertFrame])[dit].mult(-1);
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
//  if (m_opt.useAnalyticSource)
//  {
//    for (DataIterator dit = edgeScalTotal.dataIterator(); dit.ok(); ++dit)
//    {
//      FluxBox& thisFlux = edgeScalTotal[dit];
//
//      for (BoxIterator bit = BoxIterator(thisFlux.box()); bit.ok(); ++bit)
//      {
//        IntVect iv = bit();
//        RealVect loc;
//        getLocation(iv, loc, m_dx);
//
//        // Horizontal
//        thisFlux[0](iv, 0) =  0; // H comp
//        thisFlux[0](iv, 1) = 0.0; // C comp
//
//        // Vertical
//        //          thisFlux[1](iv, 0) =  -m_parameters.nonDimVel*loc[1]*loc[1]; // H comp
//        thisFlux[1](iv, 0) =  -m_parameters.nonDimVel*loc[1]; // H comp
//        thisFlux[1](iv, 1) = 0.0; // C comp
//      }
//
//    }
//  }


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
  IntVect advect_grow = m_numGhostAdvection*IntVect::Unit;
  LevelData<FArrayBox> diffusiveSrc(m_grids, 1, a_edgeScal.ghostVect()+IntVect::Unit); // this needs the ghost cell to cover slope boxes

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



    LevelData<FArrayBox> scalar_advection_old(m_grids, 1, advect_grow);

    fillScalars(scalar_advection_old, a_old_time, a_advectionVar,
                true, //do interior
                (m_opt.CFinterpOrder_advection==2) // quad interp - this seems to fix previous issues at insulating side walls
    );

    computeScalarAdvectiveFlux(a_edgeScal,scalar_advection_old, diffusiveSrc,
                                                        a_advVel, a_advectionVar,
                                                        a_old_time, a_dt);

}

void AMRLevelMushyLayer::computeScalarAdvectiveFlux(LevelData<FluxBox>& a_edgeScal, LevelData<FArrayBox>& a_scalar_advection_old,
                                                    LevelData<FArrayBox>& a_src,
                                                    LevelData<FluxBox>& a_advVel,
                                                    int a_advectionVar,
                                                    Real a_old_time, Real a_dt)
{



//void AMRLevelMushyLayer::computeScalarAdvectiveFlux(LevelData<FluxBox>& a_edgeScal, int a_advectionVar, int a_diffusionVar,
//                                                    LevelData<FluxBox>& a_advVel,
//                                                    Real a_old_time, Real a_dt)
//{
  // Need grown version of this
  IntVect advect_grow = m_numGhostAdvection*IntVect::Unit;
//  LevelData<FArrayBox> scalar_advection_old(m_grids, 1, advect_grow);
  LevelData<FArrayBox> vel(m_grids, SpaceDim, advect_grow);
//  LevelData<FArrayBox> diffusiveSrc(m_grids, 1, a_edgeScal.ghostVect()+IntVect::Unit); // this needs the ghost cell to cover slope boxes

  EdgeToCell(a_advVel, vel);
//  fillScalars(scalar_advection_old, a_old_time, a_advectionVar,
//              true, //do interior
//              (m_opt.CFinterpOrder_advection==2) // quad interp - this seems to fix previous issues at insulating side walls
//  );

  // Make diffusive source
//  bool doDiffusionSrc = true;
//  if (a_diffusionVar > -1 && doDiffusionSrc)
//  {
//    LevelData<FArrayBox>* crseScalarDiffusion = NULL;
//
//    if (m_level > 0)
//    {
//      // allocate crseBC info
//      AMRLevelMushyLayer* mlCrse = getCoarserLevel();
//
//      const DisjointBoxLayout& crseGrids = mlCrse->m_grids;
//      crseScalarDiffusion = new LevelData<FArrayBox>(crseGrids,1);
//
//      mlCrse->fillScalars(*crseScalarDiffusion, a_old_time, a_diffusionVar, true);
//    }
//
//    // Get something like grad^2(T) or div(chi dot grad S_l)
//    computeScalDiffusion(a_diffusionVar, diffusiveSrc, a_old_time);
//  }
//  else
//  {
//    setValLevel(diffusiveSrc, 0.0);
//  }

  // determine if we have inflow or outflow
  computeInflowOutflowAdvVel();

  // Compute advective flux
  computeScalarAdvectiveFlux(a_edgeScal, a_scalar_advection_old, a_advVel,
                             m_totalAdvVel,
                             vel,
                             a_src, *m_patchGodScalars[a_advectionVar],
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
  if (m_opt.doDiffusionSrc)
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
//  else
//  {
//    pout() << "time - prev diag time = " << m_time - m_prev_diag_output << " < diag period (" <<  m_opt.diagnostics_period << ")" << endl;
//  }




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


        // Either use the explicity discretization, or the advective fluxes computed by our algorithm
        // both seem to give very similar results (equal to around 4 significant figures)
        bool useExplicitDiscretization = false;

        if (useExplicitDiscretization)
        {
          FArrayBox& U = m_advVel[dit][0];
          FArrayBox& T = (*m_scalarNew[ScalarVars::m_temperature])[dit];
          Box TBox = T.box();

          // Check if we can get cells from the box to the right
          // note that b is face centred, so we only have to shift by 1/2
          Box bRight(b);
          bRight.enclosedCells(0);
          if (TBox.contains(bRight))
          {
            for (BoxIterator bit(b); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              IntVect ivUp = iv + BASISV(0);
              Real cellNu = U(iv)*( T(ivUp) + T(iv))/2 - ( T(ivUp) - T(iv))/m_dx ;
              thisNu += cellNu;
            }
          }
        }
        else if (numBoxCells > 0)
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
        FArrayBox& vel = (*m_vectorNew[VectorVars::m_fluidVel])[dit];

        Box b = m_grids[dit];
        b &= horizBox;

        if (b.size() >= IntVect::Unit)
        {
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

        FArrayBox& vel = (*m_vectorNew[VectorVars::m_fluidVel])[dit];
        if (b.size() >= IntVect::Unit)
        {
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
    averageVerticalFlux.copyTo(Interval(0, 0), *m_scalarNew[ScalarVars::m_averageHeatFlux], Interval(0,0));
    averageVerticalFlux.copyTo(Interval(1, 1), *m_scalarNew[ScalarVars::m_averageVerticalFlux], Interval(0,0));



    // Also want solute flux at each point in space written out
    int soluteComp = 1;
    for (DataIterator dit = totalFlux.dataIterator(); dit.ok(); ++dit)
    {
      FluxBox& flux = totalFlux[dit];
      FArrayBox& fab = flux[SpaceDim-1];
      Box b = flux.box();

      b.growDir(1, Side::Hi, -1); // need this because we also grab the flux in the cell one up

      Box b2 = (*m_scalarNew[ScalarVars::m_verticalFlux])[dit].box();
      b &= b2; // ensure we don't try and fill cells which don't exist
      for (BoxIterator bit(b); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        IntVect ivUp = iv + BASISV(1);
        (*m_scalarNew[ScalarVars::m_verticalFlux])[dit](iv) = 0.5*(fab(iv, soluteComp) + fab(ivUp, soluteComp));
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
//      LevelData<FArrayBox>& averageVertFs = *m_scalarNew[ScalarVars::m_averageVerticalFlux];
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
      //    Real domainSize = m_domainWidth*m_domainHeight;
      Real volume = 0.0;
      Real domainAverageFs = computeSum(volume, horizAvFs, refRat, lev0Dx, Interval(0,0), 0);
      //    domainAverageFs = domainAverageFs / domainSize;
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

      //    Real salt_mismatch = (ds_diff- flux_diff)/Fs_top;
      //    Real heat_mismatch = (H_diff - heat_flux_diff)/Fh_top;

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

    } // end if level 0

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

    Real sumLambda = ::computeNorm(amrLambda,nRef,m_dx, Interval(0,0), 2, m_level) ;
    Real maxLambda = ::computeMax(amrLambda, nRef, Interval(0,0), m_level); //::computeSum(amrLambda,nRef,m_dx, Interval(0,0), m_level);

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
    //    pout() << "AMRLevelMushyLayer::computeDiagnostics - compute average salinity" << endl;
    Vector<Real> averagedSalinity;
    horizontallyAverage(averagedSalinity, *m_scalarNew[ScalarVars::m_liquidConcentration]);

    // Also do average of liquid salinity over liquid regions
    Real averageSl = averageOverLiquidRegion(ScalarVars::m_liquidConcentration);


    int y_0 = round(averagedSalinity.size()*0.2*(1-1));
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
    Real mushAvBulkC = 0.0;
    Real mushAvPorosity = 0.0;
    Real mushVol = 0.0;
    int numMushyCells = 0;

//    int lo_j = m_problem_domain.domainBox().smallEnd()[1];

    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& porosity = (*m_scalarNew[ScalarVars::m_porosity])[dit];
      FArrayBox& bulkConc = (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit];

      for (BoxIterator bit = BoxIterator(m_grids[dit]); bit.ok(); ++bit)
      {
        IntVect iv = bit();

        RealVect loc;
         ::getLocation(iv, loc, m_dx);

//        bool is_sea_ice = (iv[1] - lo_j) > depth_i;
         bool is_sea_ice = loc[1] > (m_domainHeight-depth);

//        if (porosity(iv) < 1.0)
        if (is_sea_ice)
        {
          numMushyCells++;
          mushAvBulkC += bulkConc(iv);
          mushAvPorosity += porosity(iv);
        }
      }
    }
    mushAvBulkC = mushAvBulkC / numMushyCells;
    mushAvPorosity = mushAvPorosity / numMushyCells;
    mushVol = numMushyCells*m_dx*m_dx;

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

void AMRLevelMushyLayer::getTotalFlux(LevelData<FluxBox>& totalFlux)
{
  int numComp = totalFlux.nComp();
  IntVect fluxGhost = totalFlux.ghostVect();

  LevelData<FluxBox> diffusiveTSlFlux(m_grids, numComp, fluxGhost);

  // Get diffusive salt and heat fluxes
  RefCountedPtr<AMRNonLinearMultiCompOp> HCOp = RefCountedPtr<AMRNonLinearMultiCompOp>(
      (AMRNonLinearMultiCompOp*)this->m_HCOpFact->AMRnewOp(m_problem_domain));


  LevelData<FArrayBox> HC(m_grids, numComp, fluxGhost + IntVect::Unit); // Need one more ghost cell than the flux
  fillHC(HC, m_time);

  // diffusive flux
  for (DataIterator dit = m_scalarNew[ScalarVars::m_enthalpy]->dataIterator(); dit.ok(); ++dit)
  {
    HCOp->getFlux(diffusiveTSlFlux[dit], HC, diffusiveTSlFlux[dit].box(), dit(), 1.0);

    // Copy the vertical component to an farray box
    //    (*m_scalarNew[ScalarVars::m_FsVertDiffusion])[dit].copy(diffusiveTSlFlux[dit], 1, 0, 1);

    EdgeToCell(diffusiveTSlFlux[dit], 1, (*m_scalarNew[ScalarVars::m_FsVertDiffusion])[dit], 0, 1);
    (*m_scalarNew[ScalarVars::m_FsVertDiffusion])[dit].mult(-1);

    EdgeToCell(diffusiveTSlFlux[dit], 1, (*m_vectorNew[VectorVars::m_FsDiffusion])[dit], 0, 0);
    EdgeToCell(diffusiveTSlFlux[dit], 1, (*m_vectorNew[VectorVars::m_FsDiffusion])[dit], 1, 1);


  }




  // Set F = F_{fluid} + F_{frame}
  computeTotalAdvectiveFluxes(totalFlux);



  for (DataIterator dit = m_scalarNew[ScalarVars::m_bulkConcentration]->dataIterator(); dit.ok(); ++dit)
  {
    for (int idir=0; idir<SpaceDim; idir++)
    {
      Box b = totalFlux[dit][idir].box();

      totalFlux[dit][idir].minus(diffusiveTSlFlux[dit][idir], 0, 0, numComp);


      EdgeToCell(totalFlux[dit], 1, (*m_vectorNew[VectorVars::m_Fs])[dit], idir, idir);


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
      if ((*m_scalarNew[ScalarVars::m_porosity])[dit](iv) > 0.999)
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
      Real porosity = (*m_scalarNew[ScalarVars::m_porosity])[dit](iv);

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
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_chimneySpacing, m_time, averageChannelSpacing);
    m_diagnostics.addDiagnostic(DiagnosticNames::diag_chimneyWidth, m_time, averageChannelWidth);
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
  if (m_opt.fixedDt > 0 && !m_opt.useSubcycling)
  {
    m_dt = m_opt.fixedDt;
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
    grownDt = m_dt*m_opt.max_dt_growth;
  }

  Real solnDt = computeDt(m_opt.cfl);

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
        //TODO - 3D: work out how to make this work in 3d
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
    maxAdvU = ::computeNorm(*m_vectorNew[VectorVars::m_fluidVel], NULL, 1, m_dx, Interval(0,SpaceDim-1), 0);
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

Real AMRLevelMushyLayer::computeMaxUChi()
{
  LevelData<FArrayBox> U_chi(m_grids, SpaceDim, IntVect::Zero);
  fillVectorField(U_chi, m_time, m_U_porosity, true);

  Real maxUChi = computeMax(U_chi, NULL, -1, Interval(0,SpaceDim-1));
  return maxUChi;
}

Real AMRLevelMushyLayer::getMaxVelocityForCFL()
{
  Real maxAdvU = getMaxVelocity();

  // If we're doing advection with u/chi as the advection velocity,
  // then we need to use it for our cfl condition
  bool considerUChi = (m_opt.advectionMethod == m_porosityInAdvection ||
      m_opt.advectionMethod == m_porosityOutsideAdvection ||  m_opt.advectionMethod == m_noPorosity);

  if (m_opt.forceUseUChiForCFL)
  {
    considerUChi = true;
  }

  Real maxUChi = 0.0;

  if (considerUChi)
  {
    maxUChi = computeMaxUChi();
    if (maxUChi > 1e10)
    {
      maxUChi = 1e-100;
    }
  }

  if (s_verbosity >= 4)
  {
    pout() << "  Max(U) = " << maxAdvU << ", max(U/chi) = " << maxUChi << endl;
  }

  if (maxUChi < 1000*maxAdvU)
  {
    maxAdvU = max(maxAdvU, maxUChi);
  }

  return maxAdvU;

}

Real AMRLevelMushyLayer::computeDt(Real cfl)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRlevelMushyLayer::computeDt, cfl = " << cfl << endl;
  }

  if (m_grids.size() == 0)
  {
    return 0;
  }

  Real maxAdvU = getMaxVelocityForCFL();

  Real newDT = cfl * m_dx  / maxAdvU;

  if (s_verbosity >= 4)
  {
    pout() << "computeDt() new dt = cfl*dx/maxAdvU = " << cfl << " * " << m_dx << "/" << maxAdvU << endl;
  }

  AMRLevelMushyLayer* ml = this;
  while (ml->getFinerLevel())
  {
    ml = ml->getFinerLevel();
  }

  Real finest_dx = ml->m_dx;

  Real maxPorosity = ::computeMax(*m_scalarNew[ScalarVars::m_porosity], NULL, -1, Interval(0,0));

  Real maxUChi = computeMaxUChi();

  Real buoyancy_acceleration = abs(max(m_parameters.m_buoyancyTCoeff, m_parameters.m_buoyancySCoeff) * maxPorosity);
  Real darcy_acceleration = abs(m_parameters.m_darcyCoeff*maxUChi); //*maxPorosity/minPerm
  Real viscous_acceleration = abs(m_parameters.m_viscosityCoeff*maxAdvU/(m_domainHeight));
  //  Real acceleration = max(buoyancy_acceleration, darcy_acceleration);
  //    acceleration = max(acceleration, viscous_acceleration);

  Real acceleration = computeNorm(*m_vectorNew[VectorVars::m_advectionSrc], NULL, 1, m_dx, Interval(0, SpaceDim-1), 0);

  // Ignore bogus values
  if (abs(acceleration) > 1e100)
  {
    acceleration = 0.0;
  }

  // Could just compute norm of advection velocity solve source term?

  //  acceleration = max(acceleration);

  // Factor of 10 is fairly arbitrary
  Real accelCFL = cfl;
  if (m_opt.accelCFL > 0)
  {
    accelCFL = m_opt.accelCFL;
  }
  Real accelDt = sqrt(accelCFL*finest_dx/acceleration);

  if (m_opt.printAccelDt && s_verbosity >= 2)
  {
    pout() << "  Max dt computed from acceleration = " << accelDt << endl;
    pout() << "  Accleration terms: buoyancy = " << buoyancy_acceleration << ", viscous = " << viscous_acceleration << ", darcy = " << darcy_acceleration << endl;

  }

  // make sure we consider acceleration for determinig dt at time 0.
  // The initialisation may have generated a small initial velocity, with a corresponding CFL timestep,
  // which we want to override.

  if (maxAdvU == 0 || (m_time == 0 && m_opt.useInitAccelDt) || m_opt.useAccelDt)
  {
    if (s_verbosity >= 4)
    {
      pout() << "computeDt() new dt = min(new dt, accel dt) = min(" << newDT << ", " << accelDt << ")" << endl;
    }
    newDT = min(newDT, accelDt);
  }

  if (m_opt.max_dt > 0)
  {
    newDT = min(m_opt.max_dt, newDT);
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

  Real dt = computeDt(m_initial_dt_multiplier*m_opt.cfl);

  Real max_init_dt = m_opt.max_dt*m_initial_dt_multiplier;

  if (m_opt.max_init_dt > 0)
  {
    max_init_dt = m_opt.max_init_dt;
  }

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
      (*m_vectorOld[VectorVars::m_U_porosity])[dit].copy((*m_vectorOld[VectorVars::m_fluidVel])[dit]);
      (*m_vectorNew[VectorVars::m_U_porosity])[dit].copy((*m_vectorNew[VectorVars::m_fluidVel])[dit]);

      for (int dir=0; dir<SpaceDim; dir++)
      {
        (*m_vectorOld[VectorVars::m_U_porosity])[dit].divide(porosityOld[dit], 0, dir, 1);
        (*m_vectorNew[VectorVars::m_U_porosity])[dit].divide(porosityNew[dit], 0, dir, 1);
      }
    }

    //    m_vectorOld[VectorVars::m_U_porosity]->exchange();
    //    m_vectorNew[VectorVars::m_U_porosity]->exchange();
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
      if (m_parameters.isViscous())
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
        BCHolder viscousBC = m_physBCPtr->velFuncBC(idir, m_opt.viscousBCs, Interval(idir, idir) );

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


void AMRLevelMushyLayer::smoothEnthalpyBulkConc(Real a_smoothing)
{
  this->smoothScalarField(*m_scalarNew[ScalarVars::m_enthalpy], ScalarVars::m_enthalpy, a_smoothing);
  this->smoothScalarField(*m_scalarNew[ScalarVars::m_bulkConcentration], ScalarVars::m_bulkConcentration, a_smoothing);
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


  diffusionSolver.m_verbosity = m_opt.AMRMultigridVerb;
  diffusionSolver.m_eps = m_opt.viscous_solver_tol;
  diffusionSolver.m_normThresh = m_opt.viscous_solver_tol;
  diffusionSolver.m_hang = m_opt.viscous_solver_tol;

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

void AMRLevelMushyLayer::doRegularisationOps(LevelData<FluxBox>& a_scal, int a_var, int a_comp)
{
  if (a_var == m_porosity || a_var == m_permeability || a_var == m_bulkConcentration)
  {
    DataIterator dit2 = a_scal.dataIterator();

    for (dit2.reset(); dit2.ok(); ++dit2)
    {
      Box b = a_scal[dit2].box();

      for (int dir=0; dir<SpaceDim; dir++)
      {
        if (m_opt.useFortranRegularisationFace)
        {
          doRegularisationOpsNew(a_var, a_scal[dit2][dir], a_comp);
        }
        else
        {
          doRegularisationOps(a_var, a_scal[dit2][dir], a_comp);
        }
      }
    }
  }

}

void AMRLevelMushyLayer::doRegularisationOps(int a_var, FArrayBox& a_state, int a_comp)
{
  CH_TIME("AMRLevelMushyLayer::doRegularisationOpsOld");

  Box b = a_state.box();

  if (a_var == m_porosity)
  {
    // Ensure porosity is greater than 0 and less than or equal to 1
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      a_state(iv, a_comp) = max(m_opt.lowerPorosityLimit, a_state(iv));
      a_state(iv, a_comp) = min(1.0, a_state(iv));
    }

  }
  else if (a_var == m_permeability)
  {
    //Real minPermeability = m_parameters.calculatePermeability(m_lowerPorosityLimit);
    Real minPermeability = pow(m_opt.lowerPorosityLimit,3);
    // Ensure porosity is greater than 0
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();

      a_state(iv,a_comp) = max(minPermeability, a_state(iv));
    }

  }
  else if (a_var == ScalarVars::m_bulkConcentration)
  {
    // Was hoping the enthalpy-concentration update would guarantee this, but maybe not

    //Real minPermeability = m_parameters.calculatePermeability(m_lowerPorosityLimit);
    Real minVal = -m_parameters.compositionRatio;
    Real maxVal = 0;


    // Ensure porosity is greater than 0
    for (BoxIterator bit(b); bit.ok(); ++bit)
    {
      IntVect iv = bit();

      a_state(iv, a_comp) = max(minVal, a_state(iv));
      a_state(iv, a_comp) = min(maxVal, a_state(iv));
    }


  }


}

void AMRLevelMushyLayer::doRegularisationOpsNew(int a_var, FArrayBox& a_state, int a_comp)
{
  CH_TIME("AMRLevelMushyLayer::doRegularisationOpsNew");
  Box region = a_state.box();

  if (a_var == ScalarVars::m_porosity)
  {

    Real maxVal = 1.0;
    Real minVal = m_opt.lowerPorosityLimit;


    FORT_SETMINMAXVAL( CHF_FRA(a_state),
                       CHF_BOX(region),
                       CHF_CONST_REAL(minVal),
                       CHF_CONST_REAL(maxVal),
                       CHF_INT(a_comp));
  }
  else if (a_var == ScalarVars::m_permeability)
  {
    //Real minPermeability = m_parameters.calculatePermeability(m_lowerPorosityLimit);
    Real minPermeability = pow(m_opt.lowerPorosityLimit,3);

    FORT_SETMINVAL( CHF_FRA(a_state),
                    CHF_BOX(region),
                    CHF_CONST_REAL(minPermeability),
                    CHF_INT(a_comp));

  }
  else if (a_var == ScalarVars::m_bulkConcentration)
  {
    // Was hoping the enthalpy-concentration update would guarantee this, but maybe not

    Real minVal = -m_parameters.compositionRatio;
    Real maxVal = 0;

    FORT_SETMINMAXVAL( CHF_FRA(a_state),
                       CHF_BOX(region),
                       CHF_CONST_REAL(minVal),
                       CHF_CONST_REAL(maxVal),
                       CHF_INT(a_comp));

  }


}

void AMRLevelMushyLayer::doRegularisationOps(LevelData<FArrayBox>& a_scal,
                                             int a_var,
                                             int a_comp)
{
  CH_TIME("AMRLevelMushyLayer::doRegularisationOps");
  DataIterator dit2 = a_scal.dataIterator();

  for (dit2.reset(); dit2.ok(); ++dit2)
  {

    if (m_opt.useFortranRegularisation)
    {
      doRegularisationOpsNew(a_var, a_scal[dit2], a_comp);
    }
    else
    {
      doRegularisationOps(a_var, a_scal[dit2], a_comp);
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
                                       int& a_explicitBC, Real a_explicitVal)
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
  else if (a_var == ScalarVars::m_activeScalar)
  {
    for (int dir=0; dir<SpaceDim; dir++)
    {
      bcValHi[dir] = m_parameters.activeTracerInitVal;
      bcValLo[dir] = m_parameters.activeTracerInitVal;
    }

    return  m_physBCPtr->scalarTraceIBC(bcValLo, bcValHi,
                                        bcTypeLo, bcTypeHi);

  }
  else if (a_var == ScalarVars::m_passiveScalar)
  {
    for (int dir=0; dir<SpaceDim; dir++)
    {
      bcValHi[dir] = m_parameters.passiveTracerInitVal;
      bcValLo[dir] = m_parameters.passiveTracerInitVal;
    }

    return  m_physBCPtr->scalarTraceIBC(bcValLo, bcValHi,
                                        bcTypeLo, bcTypeHi);

  }
  else
  {
    pout() << "WARNING, NO IBCs for scalar var " << a_var << endl;

    MayDay::Warning("WARNING, NO IBCs for scalar var");
    return  m_physBCPtr->scalarTraceIBC(bcValLo, bcValHi,
                                        bcTypeLo, bcTypeHi);
  }
}



void AMRLevelMushyLayer::getScalarBCs(BCHolder& thisBC, int a_var, bool a_homogeneous)
{

  computeInflowOutflowAdvVel();
  // Reset this
  m_physBCPtr->setAdvVel(&m_totalAdvVel);


  if (a_var == ScalarVars::m_temperature)
  {
    thisBC = m_physBCPtr->BasicthetaFuncBC(a_homogeneous, &m_totalAdvVel);
  }
  else if (a_var == ScalarVars::m_enthalpy)
  {
    thisBC = m_physBCPtr->BasicEnthalpyFuncBC(a_homogeneous, &m_totalAdvVel);
  }
  else if(a_var == ScalarVars::m_enthalpySolidus)
  {
//    thisBC = m_physBCPtr->BasicSolidusBC(a_homogeneous, &m_totalAdvVel);
    thisBC = m_physBCPtr->noFluxBC();
  }
  else if(a_var == ScalarVars::m_enthalpyEutectic)
  {
//    thisBC = m_physBCPtr->BasicEutecticBC(a_homogeneous, &m_totalAdvVel);
    thisBC = m_physBCPtr->noFluxBC();
  }
  else if(a_var == ScalarVars::m_enthalpyLiquidus)
  {
//    thisBC = m_physBCPtr->BasicLiquidusBC(a_homogeneous, &m_totalAdvVel);
//    thisBC = m_physBCPtr->extrapFuncBC();
    thisBC = m_physBCPtr->noFluxBC();
  }
  else if (a_var == ScalarVars::m_porosity)
  {
    thisBC = m_physBCPtr->BasicPorosityFuncBC(a_homogeneous);
  }
  else if (a_var == ScalarVars::m_pressure)
  {
    thisBC = m_physBCPtr->BasicPressureFuncBC(a_homogeneous);
  }
  else if (a_var == ScalarVars::m_permeability)
  {
    thisBC = m_physBCPtr->BasicPermeabilityFuncBC(a_homogeneous);
  }
  else if (a_var == ScalarVars::m_liquidConcentration)
  {
    thisBC = m_physBCPtr->ThetaLFuncBC(a_homogeneous, &m_totalAdvVel);
  }
  else if (a_var == ScalarVars::m_bulkConcentration)
  {
    thisBC = m_physBCPtr->ThetaFuncBC(a_homogeneous, &m_totalAdvVel);
  }
  else if (a_var == ScalarVars::m_lambda)
  {
    thisBC = m_physBCPtr->lambdaBC(a_homogeneous);
  }
  else if (a_var == ScalarVars::m_lambda_porosity)
  {
    thisBC = m_physBCPtr->lambdaBC(a_homogeneous,
                                   true); // true - scale with porosity
  }
  else if (a_var == ScalarVars::m_viscosity)
  {
    thisBC = m_physBCPtr->extrapFuncBC();
  }
  else if (a_var == ScalarVars::m_activeScalar)
  {
    thisBC = m_physBCPtr->TracerBC(a_homogeneous, &m_totalAdvVel, m_parameters.activeTracerInitVal);
  }
  else if(a_var == ScalarVars::m_passiveScalar)
  {
//    thisBC = m_physBCPtr->noFluxBC();
    thisBC = m_physBCPtr->TracerBC(a_homogeneous, &m_totalAdvVel, m_parameters.passiveTracerInitVal);
  }
  else
  {
    pout() << "WARNING No BCs for " << m_scalarVarNames[a_var] << endl;
    MayDay::Warning("AMRLevelMushyLayer::getScalarBCs - No BCs for scalar variable, using extrapolation ");

    thisBC = m_physBCPtr->extrapFuncBC();
  }
}

//void AMRLevelMushyLayer::stokesDarcyForcing(LevelData<FArrayBox>& T, Real time)
//{
//
//  if (m_opt.stokesDarcyForcingTimescale > 0)
//  {
//
//    for (DataIterator dit=T.dataIterator(); dit.ok(); ++dit)
//    {
//      Box b = T[dit].box();
//      for (BoxIterator bit(b); bit.ok(); ++bit)
//      {
//        IntVect iv = bit();
//        RealVect loc;
//        getLocation(iv, loc, m_dx);
//
//        Real one = 1.0;
//        Real zero = 0.0;
//
//
//        Real sinPart = sin(M_PI*loc[0]);
//        Real thisT = min(one,  (time/m_opt.stokesDarcyForcingTimescale)-sinPart);
//
//        Real powPart = pow(loc[0]-2.5,2);
//        thisT = exp(-powPart/(time/m_opt.stokesDarcyForcingTimescale) ) ;
//        //			thisT = log(1+time/timescale)*exp(-pow(loc[0]-2.5,2)) ;
//
//        //			thisT = 0.01;
//
//        thisT = max(zero, thisT);
//        thisT = min(one, thisT);
//        T[dit](iv) = thisT;
//      }
//    }
//  }
//
//}

bool AMRLevelMushyLayer::crashed()
{

  Vector<int> vars;
  vars.append(ScalarVars::m_enthalpy);
  vars.append(ScalarVars::m_vorticity);
  vars.append(ScalarVars::m_enthalpySrc);
  vars.append(ScalarVars::m_lambda);

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

bool AMRLevelMushyLayer::isCurrentCFLSafe(bool printWarning)
{

  Real maxAdvU = getMaxVelocity();
  Real newCFL = maxAdvU * m_dt / m_dx;

  // CFL above 1 definitely unstable. CFL > 2*user limit means the velocity has increased rapidly
  if (newCFL > 1.0 || newCFL > m_opt.cfl*2)
  {
    if (printWarning)
    {
      pout() << "WARNING: new max(U) = " << maxAdvU << " means that CFL = " << newCFL << " which may be unstable" << endl;

    }

    if (m_opt.skipUnsafeCFL)
    {
      pout() << "Therefore, we will be skipping fluid advection during this timestep" << endl;
      return false;
    }
    else
    {
      return true;
    }
  }

  return true;
}

void AMRLevelMushyLayer::set_compute_diagnostics(bool compute_diags)
{
  m_opt.computeDiagnostics = true;
}




/*******/

#include "NamespaceFooter.H"



