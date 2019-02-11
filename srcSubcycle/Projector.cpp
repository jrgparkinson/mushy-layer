
#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Projector.H"
#include "Gradient.H"
#include "GradientF_F.H"
#include "Divergence.H"
#include "AMRPoissonOp.H"
#include "VCAMRPoissonOp2.H"
//#include "MyVCAMRPoissonOp2CC.H"
#include "AMRProjectionOp.h"
#include "RelaxSolver.H"
#include "CoarseAverageEdge.H"
#include "CornerCopier.H"
#include "AMRIO.H"
#include "computeSum.H"
#include "SetValLevel.H"
#include "CoarseAverage.H"
#include "VelBCHolder.H"
#include "FineInterp.H"
#include "EdgeToCell.H"
#include "RelaxSolver.H"
#include "computeNorm.H"
#include "BoxIterator.H"

#ifdef CH_USE_HDF5
#include "CH_HDF5.H"
#endif

/// define static variables here
bool Projector::m_doSyncProjection = true;
bool Projector::m_applySyncCorrection = true;
bool Projector::s_applyVDCorrection = true;
bool Projector::m_doQuadInterp = true;
Real Projector::m_etaLambda = 0.9;

#if defined(CH_USE_DOUBLE)
Real Projector::s_solver_tolerance = 1.0e-15; // pmc, 18 oct 2007:  was -10
#elif defined(CH_USE_FLOAT)
Real Projector::s_solver_tolerance = 1.0e-5; // pmc, 18 oct 2007:  was -10
#endif

Real Projector::s_solver_hang = 1e-15;

int  Projector::s_num_smooth_up = 4;
int  Projector::s_num_smooth_down = 2;
int  Projector::s_num_precond_smooth = 4;
int  Projector::s_numMG = 1;
bool  Projector::s_relax_bottom_solver = false;
bool Projector::s_constantLambdaScaling = false;
Real Projector::s_lambda_timestep = 0.0;
bool Projector::pp_init = false;
int  Projector::s_verbosity = 2;
int Projector::s_bottomSolveMaxIter = 20;
int Projector::s_multigrid_relaxation = 1; // 1 for gsrb, 4 for jacobi

/// first define quick-n-easy access functions

// ------------------------------------------------------
int Projector::getLevel() const
{
  return m_level;
}

// ------------------------------------------------------
const ProblemDomain& Projector::dProblem() const
{
  return m_domain;
}

// ------------------------------------------------------
const DisjointBoxLayout& Projector::getBoxes() const
{
  return m_Pi.getBoxes();
}

// ------------------------------------------------------
int Projector::nRefCrse() const
{
  return m_nRefCrse;
}

// ------------------------------------------------------
Real Projector::dx() const
{
  return m_dx;
}

// ------------------------------------------------------
Real Projector::etaLambda() const
{
  return m_etaLambda;
}

void Projector::etaLambda(Real a_eta)
{
  m_etaLambda = a_eta;
}

// ------------------------------------------------------
Projector* Projector::fineProjPtr() const
{
  return m_fineProjPtr;
}

// ------------------------------------------------------
Projector* Projector::crseProjPtr() const
{
  return m_crseProjPtr;
}

// ------------------------------------------------------
LevelData<FArrayBox>& Projector::phi()
{
  return m_phi;
}

void Projector::getPhi(LevelData<FArrayBox>& a_phi)
{
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
  {
    a_phi[dit].copy(m_phi[dit]);
  }
}
// ------------------------------------------------------
const LevelData<FArrayBox>& Projector::phi() const
{
  return m_phi;
}

// ------------------------------------------------------
LevelData<FArrayBox>& Projector::eSync()
{
  return m_eSync;
}

// ------------------------------------------------------
const LevelData<FArrayBox>& Projector::eSync() const
{
  return m_eSync;
}

// ------------------------------------------------------
LevelData<FArrayBox>& Projector::MACrhs()
{
  return m_MACrhs;
}

// ------------------------------------------------------
LevelData<FArrayBox>& Projector::CCrhs()
{
  return m_CCrhs;
}

LevelData<FluxBox>& Projector::MACcorrection()
{
  return m_MACcorrection;
}

LevelData<FArrayBox>& Projector::CCcorrection()
{
  return m_CCcorrection;
}



LevelData<FArrayBox>& Projector::Pi()
{
  return m_Pi;
}

void Projector::unscaledPi(LevelData<FArrayBox>& unscaledPi, Real a_dt)
{
  m_Pi.copyTo(unscaledPi);
  for (DataIterator dit = m_Pi.dataIterator(); dit.ok(); ++dit)
  {
    unscaledPi[dit].mult(a_dt);
  }

}

// ------------------------------------------------------
LevelData<FArrayBox>& Projector::eLambda()
{
  return m_eLambda;
}

// ------------------------------------------------------
const LevelData<FArrayBox>& Projector::eLambda() const
{
  return m_eLambda;
}

// ------------------------------------------------------
LevelData<FluxBox>& Projector::grad_eLambda()
{
  return m_grad_eLambda;
}

// ------------------------------------------------------
const LevelData<FluxBox>& Projector::grad_eLambda() const
{
  return m_grad_eLambda;
}

// ------------------------------------------------------
bool Projector::doSyncProjection() const
{
  return m_doSyncProjection;
}

// ------------------------------------------------------
bool Projector::doMacSync() const
{
  return (m_etaLambda != 0.0);
}

// ------------------------------------------------------
bool Projector::isInitialized() const
{
  return m_isInitialized;
}

// ------------------------------------------------------
bool Projector::doQuadInterp() const
{
  return m_doQuadInterp;
}

// ------------------------------------------------------
QuadCFInterp& Projector::quadCFInterpolator()
{
  return m_cfInterp;
}

/// returns predefined intvectset which shows coverage
const LayoutData<IntVectSet>& Projector::getGridIVS()
{
  return m_gradIVS;
}

// ------------------------------------------------------
bool Projector::isFinestLevel() const
{
  return m_finest_level;
}

// ------------------------------------------------------
void Projector::isFinestLevel(bool a_finest_level)
{
  m_finest_level = a_finest_level;
}

// ------------------------------------------------------
void Projector::verbosity(int a_verbosity)
{
  s_verbosity = a_verbosity;
}
//
//void CCProjector::setFineProjPtr(CCProjector* fineProj)
//{
//  m_fineProjPtr = fineProj;
//}

// ------------------------------------------------------
int Projector::verbosity() const
{
  return s_verbosity;
}

// ------------------------------------------------------
void Projector::setPhysBC(PhysBCUtil& a_bc)
{
  //  if (m_physBCPtr != NULL)
  //    {
  //      delete m_physBCPtr;
  //      m_physBCPtr = NULL;
  //    }
  //  m_physBCPtr = a_bc.newPhysBCUtil();
  m_physBCPtr = (&a_bc);
}

// ------------------------------------------------------
PhysBCUtil* Projector::getPhysBCPtr() const
{
  return m_physBCPtr;
}

// ------------------------------------------------------
// default constructor
Projector::Projector()
{
  m_level = -1;
  m_nRefCrse = -1;
  m_isInitialized = false;
  // have to initialize this to _something_!
  m_finest_level = false;
  m_physBCPtr = NULL;
  m_bottomSolverLevel = NULL;
  m_limitSolverCoarsening = false;
  m_crseProjPtr = NULL;
  m_dx = -1;
  m_fineProjPtr = NULL;
  m_BiCGBottomSolverLevel = NULL;
  m_sumVDrhs = 0.0;
  m_usePiAdvectionBCs = true;
  m_phiScale=1.0;
  m_MACbcScale=-0.5;
  m_scale_lambda_with_porosity = false;
  m_scale_lambda_err_with_porosity = false;
  m_scaleSyncCorrection=true;
  m_scalePressureWithPorosity=true;
}

// ------------------------------------------------------
// destructor
Projector::~Projector()
{
  //We don't own this anymore
  //  if (m_physBCPtr != NULL)
  //    {
  //      delete m_physBCPtr;
  //      m_physBCPtr = NULL;
  //    }

  if(m_bottomSolverLevel != NULL) {

    delete m_bottomSolverLevel;
    m_bottomSolverLevel = NULL;
  }
  // everything else should be automatic here
}

// ------------------------------------------------------
Projector::Projector(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_crseGridsPtr,
                         const Box& a_domain,
                         const Real a_dx,
                         Projector* a_finerProjPtr,
                         Projector* a_crseProjPtr,
                         int a_nRefCrse,
                         int a_level,
                         PhysBCUtil& a_physBC)
{
  m_physBCPtr = NULL;
  ProblemDomain physdomain(a_domain);
  define(a_grids, a_crseGridsPtr, physdomain, a_dx, a_finerProjPtr,
         a_crseProjPtr, a_nRefCrse, a_level, a_physBC);
}

// ------------------------------------------------------
Projector::Projector(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_crseGridsPtr,
                         const ProblemDomain& a_domain,
                         const Real a_dx,
                         Projector* a_finerProjPtr,
                         Projector* a_crseProjPtr,
                         int a_nRefCrse,
                         int a_level,
                         PhysBCUtil& a_physBC)
{
  m_physBCPtr = NULL;
  define(a_grids, a_crseGridsPtr, a_domain, a_dx, a_finerProjPtr,
         a_crseProjPtr, a_nRefCrse, a_level, a_physBC);
}

void Projector::copyPressure(Projector& a_proj)
{
  // Things to copy:
  // m_Pi
  // m_phi

  // Don't do anything if this projection operator isn't defined
  if (m_level == -1)
  {
    return;
  }

  LevelData<FArrayBox>& a_Pi = a_proj.Pi();

  a_Pi.copyTo(Interval(0,0), m_Pi, Interval(0,0));
  a_proj.phi().copyTo(Interval(0,0), m_phi, Interval(0,0));

}

void Projector::rescalePressure(Real a_oldDt, Real a_newDt)
{
  //Pressure is scaled by dt, such that a_Pi = pressure/dt

  Real multiplier = a_oldDt/a_newDt;

  for (DataIterator dit = m_Pi.dataIterator(); dit.ok(); ++dit)
  {
    m_Pi[dit].mult(multiplier);
    m_phi[dit].mult(multiplier);
  }
}

// ------------------------------------------------------
void Projector::define(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_crseGridsPtr,
                         const Box& a_domain,
                         const Real a_dx,
                         Projector* a_finerProjPtr,
                         Projector* a_crseProjPtr,
                         int a_nRefCrse,
                         int a_level,
                         PhysBCUtil& a_physBC,
                         bool a_usePiAdvectionBCs)
{
  ProblemDomain physdomain(a_domain);
  define(a_grids, a_crseGridsPtr, physdomain, a_dx, a_finerProjPtr,
         a_crseProjPtr, a_nRefCrse, a_level, a_physBC, a_usePiAdvectionBCs);
}

// ------------------------------------------------------
void Projector::define(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_crseGridsPtr,
                         const ProblemDomain& a_domain,
                         const Real a_dx,
                         Projector* a_finerProjPtr,
                         Projector* a_crseProjPtr,
                         int a_nRefCrse,
                         int a_level,
                         PhysBCUtil& a_physBC,
                         bool a_usePiAdvectionBCs)
{
  // set physical boundary condition object
  setPhysBC(a_physBC);

  if (!pp_init)
  {
    variableSetUp();
  }

  CH_assert (a_level > -1);
  m_level = a_level;
  m_dx = a_dx;
  m_domain = a_domain;
  m_usePiAdvectionBCs = a_usePiAdvectionBCs;

  //m_fineProjPtr = a_finerProjPtr;
  setFineProj(a_finerProjPtr);

  m_crseProjPtr = a_crseProjPtr;
  if (m_crseProjPtr != NULL)
  {
    m_nRefCrse = a_nRefCrse;
    m_crseProjPtr->setFineProj(this);
    m_crseProjPtr->isFinestLevel(false);
    CH_assert (a_crseGridsPtr != NULL);

    m_cfInterp.define(a_grids, a_crseGridsPtr, a_dx, a_nRefCrse,
                      1,a_domain);
  }

  int nGrow = 1;
  Gradient::createGridIVS(m_gradIVS, a_grids, nGrow);

  IntVect ghostVect(IntVect::Unit);

  // define persistent storage
  m_Pi.define( a_grids,    1,2*ghostVect); // need two ghost cells here for CC projection
  m_phi.define(a_grids,    1,ghostVect);
  m_eLambda.define(a_grids,1,ghostVect);
  m_eSync.define(  a_grids,1,ghostVect);
  m_MACrhs.define( a_grids,1); // don't think these need ghost vectors
  m_CCrhs.define(  a_grids,1);
  m_MACcorrection.define(a_grids, 1);
  m_CCcorrection.define(a_grids, SpaceDim);

  m_grad_eLambda.define(a_grids,1);

  // set grad_e_lambda to be 0
  setValLevel(m_grad_eLambda, 0.0);
  // also initialize all other fields to a bogus value
  Real bogus = 1.0e300;
  setValLevel(m_Pi, bogus);
  setValLevel(m_phi, 0.0);
  // these need to start as 0 for boundary conditions
  // for first intermediate solve if more than one level
  // of refinement
  setValLevel(m_eLambda, 0.0);
  setValLevel(m_eSync, 0.0);

  defineSolverMGlevel(a_grids, a_crseGridsPtr);

  // surely is initialized by now?
  m_isInitialized = true;

  if (s_verbosity > 4)
  {
    pout() << "CCProjector::define finished" << endl;
  }
}

// -------------------------------------------------------
void Projector::init(const Projector& a_oldProj)
{
  // this function performs a partial initialization
  // of this object using the oldProj...

  // at the current state of this algorithm,
  // this doesn't need to do anything (everything will
  // need to be re-initialized anyway)
}

void Projector::variableSetUp()
{
  // first set up parm parse object
  ParmParse ppProjection("projection");

  int tempBool;

  tempBool = (int) m_doSyncProjection;
  ppProjection.query("doSyncProjection",tempBool);
  m_doSyncProjection = (tempBool == 1);

  tempBool = (int) m_applySyncCorrection;
  ppProjection.query("applySyncCorrection", tempBool);
  m_applySyncCorrection = (tempBool == 1);

  tempBool = (int) s_applyVDCorrection;
  ppProjection.query("applyFreestreamCorrection", tempBool);
  s_applyVDCorrection = (tempBool == 1);

  tempBool = (int) m_doQuadInterp;
  ppProjection.query("doQuadInterp", tempBool);
  m_doQuadInterp = (tempBool == 1);

  m_scalePressureWithPorosity = false;
  ppProjection.query("scalePressureWithPorosity", m_scalePressureWithPorosity);

  m_scaleSyncCorrection=false;
  ppProjection.query("scaleSyncCorrection", m_scaleSyncCorrection);

  ppProjection.query("phiScale", m_phiScale);
  ppProjection.query("MACbcScale", m_MACbcScale);

  ppProjection.query("eta", m_etaLambda);
  ppProjection.query("scale_lambda_with_porosity", m_scale_lambda_with_porosity);
  ppProjection.query("scale_lambda_err", m_scale_lambda_err_with_porosity);

  ppProjection.query("solverTol", s_solver_tolerance);
  ppProjection.query("numSmoothUp", s_num_smooth_up);
  ppProjection.query("numSmoothDown", s_num_smooth_down);
  ppProjection.query("numPrecondSmooth", s_num_precond_smooth);
  ppProjection.query("numMG", s_numMG);
  ppProjection.query("relax_bottom_solver", s_relax_bottom_solver);
  ppProjection.query("bottomSolveMaxIter", s_bottomSolveMaxIter);
  ppProjection.query("solverHang", s_solver_hang);
  ppProjection.query("mg_relaxation", s_multigrid_relaxation);

  tempBool = (int) s_constantLambdaScaling;
  ppProjection.query("constantLambdaScaling", tempBool);
  s_constantLambdaScaling = (tempBool == 1);

  // now dump this all back out
  if (s_verbosity > 2)
  {
    pout () << "Projection inputs: " << endl;

    pout() << "  doSyncProjection = "  << m_doSyncProjection << endl;

    pout() << "  applySyncCorrection = "  << m_applySyncCorrection << endl;

    pout() << "  applyFreestreamCorrection = "  << s_applyVDCorrection
        << endl;

    pout() << "  doQuadInterp = "  << m_doQuadInterp << endl;

    pout() << "  eta = "  <<  m_etaLambda << endl;

    pout() << "  constantLambdaScaling = "  << s_constantLambdaScaling
        << endl;
  }

  // set flag so we don't have to call this function again.
  pp_init = true;

}

void Projector::setCrseProj(Projector* a_crseProjPtr, int a_nRefCrse)
{
  if (m_crseProjPtr != NULL)
  {
    // if a coarse level is defined, should never undefine it
    CH_assert (a_crseProjPtr != NULL);
    // nRefCrse should never change either.
    CH_assert (m_nRefCrse == a_nRefCrse);
    // check to see if coarse-grid info will need to be redone
    if (!(m_crseProjPtr->getBoxes()==a_crseProjPtr->getBoxes()))
    {
      // new coarse boxes -- need to redefine levelSolver
      const DisjointBoxLayout& crseGrids = a_crseProjPtr->getBoxes();
      defineSolverMGlevel(getBoxes(), &crseGrids);

    } // end if old and new grids are different
    else
    {
      // if old and new coarse grids are the same, no need to
      // do this
    }
  } // end if no coarser level
  else
  {
    // if old crse level didn't exist and new one does, then need
    // to initialize levelsolver
    if (a_crseProjPtr != NULL)
    {
      m_nRefCrse = a_nRefCrse;
      const DisjointBoxLayout& crseGrids = a_crseProjPtr->getBoxes();
      defineSolverMGlevel(getBoxes(), &crseGrids);
    } // end if there's a crse level

  } // end if new level

  // don't forget to actually _set_ the pointer!
  m_crseProjPtr = a_crseProjPtr;
}

// -----------------------------------------------------------
void Projector::setFineProj(Projector* a_fineProjPtr)
{
  // this is much simpler than resetting crseProj, since
  // there is no need to rebuild levelsolvers or anything like
  // that
  m_fineProjPtr = a_fineProjPtr;
}

void Projector::limitSolverCoarsening(bool a_limitSolverCoarsening)
{
  m_limitSolverCoarsening = a_limitSolverCoarsening;
}

#ifdef CH_USE_HDF5
// --------------------------------------------------------------
void Projector::writeCheckpointHeader(HDF5Handle& a_handle) const
{

  // nothing needs to be done here -- just
  // included for completeness...

}

// --------------------------------------------------------------
void Projector::writeCheckpointLevel(HDF5Handle& a_handle) const
{

  char level_str[20];
  sprintf (level_str, "%d", m_level);
  const std::string label = std::string("level_") + level_str
      + std::string("/projection");

  a_handle.setGroup(label);

  HDF5HeaderData header;

  header.m_real["dx"] = m_dx;
  header.m_box["domain"] = m_domain.domainBox();
  header.m_int["n_ref_crse"] = m_nRefCrse;
  if (m_finest_level)
  {
    header.m_int["finest_level"] = 1;
  }
  else
  {
    header.m_int["finest_level"] = 0;
  }

  // as in AMREuler, don't write static data (comes from ParmParse)

  header.writeToFile(a_handle);

  // now dump out this level's data
  write (a_handle, m_Pi.boxLayout());
  // also want to dump out ghost-cell vales
  write (a_handle, m_Pi, "pi", m_Pi.ghostVect());
  write (a_handle, m_phi, "phi", m_phi.ghostVect());
  write (a_handle, m_eSync, "e_sync", m_eSync.ghostVect());
  write (a_handle, m_eLambda, "e_lambda", m_eLambda.ghostVect());
  write (a_handle, m_grad_eLambda, "grad_e_lambda", m_grad_eLambda.ghostVect());
}

// --------------------------------------------------------------
void Projector::readCheckpointHeader(HDF5Handle& a_handle)
{
  //  this doesn't actually do anything...
}

// --------------------------------------------------------------
void Projector::readCheckpointLevel(HDF5Handle& a_handle)
{
  char level_str[20];
  sprintf (level_str, "%d", m_level);
  const std::string label = std::string("level_") + level_str
      + std::string("/projection");

  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout () << "hdf5 header data: " << endl;
    pout () << header << endl;
  }

  // read dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("CCProjector::readCheckpointLevel: file does not contain dx");
  }
  if (m_dx != header.m_real["dx"])
  {
    MayDay::Error("CCProjector::readCheckpointLevel: dx in checkfile does not match that in solver");
  }

  // read domain
  if (header.m_box.find("domain") == header.m_box.end())
  {
    MayDay::Error("CCProjector::readCheckpointLevel: file does not contain domain");
    if (m_domain.domainBox() != header.m_box["domain"])
    {
      MayDay::Error("CCProjector::readCheckpointLevel: domain in checkfile does not match that in solver!");
    }
  }

  // read nRefCrse
  if (header.m_int.find("n_ref_crse") == header.m_int.end())
  {
    MayDay::Error("CCProjector::readCheckpointLevel: file does not contain n_ref_crse");
    if (m_nRefCrse != header.m_int["n_ref_crse"])
    {
      MayDay::Error("CCProjector::readCheckpointLevel: n_ref_crse in checkfile does not match that in solver!");
    }
  }

  // read nRefCrse
  if (header.m_int.find("finest_level") == header.m_int.end())
  {
    MayDay::Error("CCProjector::readCheckpointLevel: file does not contain finest_level");
  }
  m_finest_level = (header.m_int["finest_level"] == 1);

  // read grids
  Vector<Box> grids;
  const int grid_status = read(a_handle, grids);
  if (grid_status != 0)
  {
    MayDay::Error("CCProjector::readCheckpointLevel: file does not contain a Vector<Box>");
  }
  // can't really check for equality, but at least check for same size...
  if (grids.size() != m_Pi.getBoxes().size())
  {
    MayDay::Error("CCProjector::readCheckpointLevel: grids in checkfile different size than in solver!");
  }

  const DisjointBoxLayout& levelGrids = m_Pi.getBoxes();
  const int piData_status = read<FArrayBox> (a_handle, m_Pi, "pi", levelGrids);

  if (piData_status != 0)
  {
    MayDay::Error("CCProjector::readCheckpointLevel: file does not contain pi data");
  }

  const int phiData_status = read<FArrayBox> (a_handle, m_phi, "phi", levelGrids);

  if (phiData_status != 0)
  {
    MayDay::Error("CCProjector::readCheckpointLevel: file does not contain phi data");
  }

  const int eSyncData_status = read<FArrayBox> (a_handle, m_eSync, "e_sync", levelGrids);

  if (eSyncData_status != 0)
  {
    MayDay::Error("CCProjector::readCheckpointLevel: file does not contain e_sync data");
  }

  const int eLambdaData_status = read<FArrayBox> (a_handle, m_eLambda, "e_lambda", levelGrids);

  if (eLambdaData_status != 0)
  {
    MayDay::Error("CCProjector::readCheckpointLevel: file does not contain eLambda data");
  }

  const int grad_eLambdaData_status = read<FluxBox> (a_handle, m_grad_eLambda, "grad_e_lambda", levelGrids);

  if (grad_eLambdaData_status != 0)
  {
    MayDay::Error("CCProjector::readCheckpointLevel: file does not contain grad_e_lambda data");
  }

  // now need to set boundary conditions on all these things
  postRestart();
}

#endif

// --------------------------------------------------------------
// note that this function assumes that all BC's on phi have been set
void Projector::gradPhi(LevelData<FArrayBox>& a_gradPhi, int a_dir) const
{
  // first compute edge-centering of gradPhi
  DataIterator dit = a_gradPhi.dataIterator();
  dit.reset();

  int edgeDir = -1;
  const Box& edgeBox = a_gradPhi[dit].box();
  for (int dir=0; dir<SpaceDim; dir++)
  {
    if (edgeBox.type(dir) == IndexType::NODE)
    {
      CH_assert (edgeDir == -1);
      edgeDir = dir;
    }
  }
  CH_assert(edgeDir != -1);

  const DisjointBoxLayout& grids = getBoxes();

  for (dit.reset(); dit.ok(); ++dit)
  {
    // call gradient subroutine directly, since this is one-directional

    // I don't think this function is called by anybody ever
    // I think this is correct, but just in case, flag it
    // (dfmartin@lbl.gov)
    CH_assert (false);

    Gradient::singleBoxMacGrad(a_gradPhi[dit],
                               m_phi[dit],
                               0, 0, 1, grids[dit], m_dx,
                               a_dir, edgeDir, m_gradIVS[dit]);
  } // end loop over grids
}

// --------------------------------------------------------------
// note that this function assumes that all BC's on phi have been set
void Projector::gradPhi(LevelData<FluxBox>& a_gradPhi) const
{
  Gradient::levelGradientMAC(a_gradPhi, m_phi, m_dx);
}

// --------------------------------------------------------------
// note that this function assumes that all BC's on Pi have been set
void Projector::gradPi(LevelData<FArrayBox>& a_gradPi, int a_dir) const
{
  DataIterator dit = a_gradPi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    // call FORTRAN subroutine directly...
    FORT_GRADCC(CHF_FRA1(a_gradPi[dit],0),
                CHF_CONST_FRA1(m_Pi[dit],0),
                CHF_BOX(getBoxes()[dit]),
                CHF_CONST_REAL(m_dx),
                CHF_INT(a_dir));
  }
}

// --------------------------------------------------------------
// note that this function assumes that all BC's on Pi have been set
void Projector::gradPi(LevelData<FArrayBox>& a_gradPi) const
{
  Gradient::levelGradientCC(a_gradPi, m_Pi, m_dx);
}


void Projector::gradPiBCs(LevelData<FArrayBox>& a_gradPi, bool extrapBCs, bool a_usePhi)
{
  //  BCHolder gradPbcHolder = m_physBCPtr->gradPiFuncBC();
  BCHolder bcHolder;

  if (extrapBCs)
  {
    bcHolder  = m_physBCPtr->extrapFuncBC();
  }
  else
  {
    bcHolder  = m_physBCPtr->BasicPressureFuncBC(false);
  }
  //  BCHolder bcHolder = m_physBCPtr->gradPiFuncBC();
  DataIterator dit = a_gradPi.dataIterator();
  const DisjointBoxLayout& levelGrids = getBoxes();

  LevelData<FArrayBox> pressureTemp(m_Pi.disjointBoxLayout(), 1, m_Pi.ghostVect());
  if ( a_usePhi)
  {
    pressureTemp[dit].copy(m_phi[dit]);
  }
  else
  {
    for (dit.reset(); dit.ok(); ++dit)
    {
      pressureTemp[dit].copy(m_Pi[dit]);
    }
  }

  for (dit.reset(); dit.ok(); ++dit)
  {
    const Box& b = levelGrids[dit];

    // Make sure we copy across ghost cells
    bcHolder.operator()(pressureTemp[dit],
                        b,
                        m_domain,
                        m_dx,
                        false); // not homogeneous
  }

  pressureTemp.exchange();

  Gradient::levelGradientCC(a_gradPi, pressureTemp, m_dx);

}

// -------------------------------------------------------------
// note that this function assumes all BC's on eSync have been set
// also only returns valid grad(eSync) in valid regions of grids
void Projector::grad_eSync(LevelData<FArrayBox>& a_grad_eSync) const
{
  if (m_fineProjPtr != NULL)
  {
    const LevelData<FArrayBox>& fine_eSync = m_fineProjPtr->eSync();
    int nRefFine = m_fineProjPtr->nRefCrse();

    Gradient::compGradientCC(a_grad_eSync, m_eSync,  &fine_eSync,
                             m_dx, nRefFine, m_domain);
  }
  else
  {
    // no finer level -- can do level-operator gradient
    Gradient::levelGradientCC(a_grad_eSync, m_eSync, m_dx);
  }
}

// --------------------------------------------------------------
// note that this function assumes all BC's on eSync have been set
// also only returns valid grad(eSync) in valid regions of grids
void Projector::grad_eSync(LevelData<FArrayBox>& a_grad_eSync, int a_dir) const
{
  // this is not yet implemented (not sure if i'll really need it!)
  MayDay::Warning("directional grad_eSync not yet implemented!");
}


void Projector::setSubcycledMACBCs()
{
  m_usePiAdvectionBCs = true;
  m_phiScale=1.0;
  m_MACbcScale=-0.5;
}

void Projector::setNonSubcycledMACBCs()
{
  m_usePiAdvectionBCs = false;
  m_phiScale=1.0;
  m_MACbcScale=1.0;
}

// This is used for projecting a time independent advection velocity
int Projector::levelMacProject(LevelData<FluxBox>& a_uEdge,
                                  Real a_dt,
                                  const RefCountedPtr<LevelData<FArrayBox> > a_crsePressurePtr,
                                  const RefCountedPtr<LevelData<FArrayBox> > a_pressureScalePtr,
                                  const RefCountedPtr<LevelData<FArrayBox> > a_crsePressureScalePtr,
                                  const RefCountedPtr<LevelData<FluxBox> > a_pressureScaleEdgePtr,
                                  const RefCountedPtr<LevelData<FluxBox> > a_crsePressureScaleEdgePtr,
                                  bool alreadyHasPhi,
                                  Real correctScale)
{
  if (s_verbosity >= 5)
  {
    pout() << "CCProjector::levelMacProject (level " << m_level << ")"    << endl;
  }


//  if (s_multigrid_relaxation == 1)
//  {
//    s_multigrid_relaxation = 4;
//  }
//  else
//  {
//    s_multigrid_relaxation = 1;
//  }
//  pout() << "  Multigrid relaxation = " << s_multigrid_relaxation << endl;

  // MAC rhs should have no ghost values
  //  LevelData<FArrayBox> MacRHS(getBoxes(),1);

  Divergence::levelDivergenceMAC(MACrhs(), a_uEdge, m_dx);

  // report sum(rhs)
  if (s_verbosity >= 2)
  {
    DisjointBoxLayout* finerGridsPtr = NULL;
    int nRefFine = -1;
    Real sumRHS = computeSum(MACrhs(), finerGridsPtr,
                             nRefFine, m_dx, MACrhs().interval());
    pout() << "  MAC projection (level " << m_level << ") -- sum(RHS) = " << sumRHS << endl;


    // Compute sum over just the interior
    DisjointBoxLayout interiorGrids(MACrhs().disjointBoxLayout());
//    interiorGrids.grow(-2);

    LevelData<FArrayBox> rhsInterior(interiorGrids, 1);
    MACrhs().copyTo(rhsInterior);

    Box domBox = m_domain.domainBox();

    for (DataIterator dit = rhsInterior.dataIterator(); dit.ok(); ++dit)
    {
      SideIterator sit;
      for (sit.reset(); sit.ok(); ++sit)
      {
        Side::LoHiSide side = sit();

        for (int dir=0; dir < SpaceDim; dir++)
        {
          int boxSize = 5;
          Box zeroBox = adjCellBox(domBox, dir, side, boxSize);
          zeroBox.shift(dir, -boxSize);

          Box stateBox = rhsInterior[dit].box();

          zeroBox &= stateBox;
          rhsInterior[dit].setVal(0.0, zeroBox, 0);
        }
      }
    }

    Real sumRHSInterior = ::computeSum(rhsInterior, finerGridsPtr,
                                       nRefFine, m_dx, MACrhs().interval());
    pout() << "  MAC projection (level " << m_level << ") -- sum(interior RHS) = " << sumRHSInterior << endl;

  }

  DataIterator dit = m_phi.dataIterator();

  LevelData<FArrayBox> oldPhi(m_phi.disjointBoxLayout(), m_phi.nComp(), m_phi.ghostVect());

  // now solve for phi
  for (dit.reset(); dit.ok(); ++dit)
  {
    oldPhi[dit].copy(m_phi[dit]);
    MACrhs()[dit].mult(correctScale);
  }

  MACrhs().exchange();

  int exitStatus = solveMGlevel(m_phi, a_crsePressurePtr, MACrhs(), a_pressureScalePtr, a_crsePressureScalePtr,
               a_pressureScaleEdgePtr, a_crsePressureScaleEdgePtr);

  // now correct

  // first re-apply physical and copy BC's
  Interval phiComps(0,0);
  m_phi.exchange(phiComps);

  //    BCHolder bcHolder = m_physBCPtr->gradMacPressureFuncBC();
  BCHolder bcHolder = m_physBCPtr->BasicPressureFuncBC(false);

  const DisjointBoxLayout& levelGrids = getBoxes();

  for (dit.reset(); dit.ok(); ++dit)
  {
    bcHolder.operator()(m_phi[dit],
                        levelGrids[dit],
                        m_domain,
                        m_dx,
                        false); // not homogeneous
  }

  //#define CORNERCOPIER
  //#ifdef CORNERCOPIER
  CornerCopier cornerExchangeCopier(m_phi.getBoxes(),
                                    m_phi.getBoxes(),
                                    m_domain,
                                    m_phi.ghostVect(),
                                    true);

  m_phi.exchange(m_phi.interval(), cornerExchangeCopier);
  //#endif

  LevelData<FluxBox> gradPhi(levelGrids,1);


  // compute gradient
  Gradient::levelGradientMAC(gradPhi,m_phi, a_crsePressurePtr, m_dx,
                             m_gradIVS,
                             m_cfInterp);

  // now rescale and subtract from uEdge

  // assumes that uEdge and phi can use same dataIterator
  for (dit.reset(); dit.ok(); ++dit)
  {
    FluxBox& thisGrad = gradPhi[dit];
    FluxBox& thisEdgeVel = a_uEdge[dit];

    for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& thisGradDir = thisGrad[dir];
      FArrayBox& thisVelDir = thisEdgeVel[dir];

      if (m_porosityEdgePtr != NULL)
      {
        thisGradDir.mult((*m_porosityEdgePtr)[dit][dir],thisGradDir.box(), 0, 0);

      }
//      thisGradDir.mult(correctScale);

      thisVelDir -= thisGradDir;
    }
  }

  if(alreadyHasPhi)
  {
    for (dit.reset(); dit.ok(); ++dit)
    {
      m_phi[dit].plus(oldPhi[dit]);
    }
  }

  checkDivergence(a_uEdge);

  return exitStatus;

}

void Projector::checkDivergence(LevelData<FluxBox>& a_uEdge)
{
  LevelData<FArrayBox> DivU(getBoxes(),1);
  Divergence::levelDivergenceMAC(DivU, a_uEdge, m_dx);

  Real maxDivU = ::computeNorm(DivU, NULL, 1, m_dx, Interval(0,0), 0);

  int minNumSmooth = 2;
  int maxNumSmooth = 20;
  Real tol = 1e-11;

  ParmParse pp("projection");
  pp.query("maxNumSmooth", maxNumSmooth);
  pp.query("minNumSmooth", minNumSmooth);

  if (maxDivU < tol)
  {
    // We might be able to reduce the number of smoothings
    if (s_num_smooth_down > minNumSmooth || s_num_smooth_up > minNumSmooth)
    {
      s_num_smooth_down = max(s_num_smooth_down-2, minNumSmooth);
      s_num_smooth_up = max(s_num_smooth_up-2, minNumSmooth);

      pout() << "  max(div U)  = " << maxDivU << " < " << tol << ", decreasing number of smoothing steps to " <<
          s_num_smooth_down << " (down) " << s_num_smooth_up << " (up) "<< endl;
    }

  }
  else
  {
    // Increasing the number of smoothings may help
    if (s_num_smooth_down < maxNumSmooth || s_num_smooth_up < maxNumSmooth)
    {
      s_num_smooth_down = min(s_num_smooth_down+2, maxNumSmooth);
      s_num_smooth_up = min(s_num_smooth_up+2, maxNumSmooth);
      pout() << "  max(div U)  = " << maxDivU << " > " << tol << ", increasing number of smoothing steps to " <<
          s_num_smooth_down << " (down) " << s_num_smooth_up << " (up) " << endl;
    }
  }

}

Real Projector::getPhiScale(Real a_dt)
{
  return getScale(m_phiScale, a_dt);
}

//Real CCProjector::getScale(const char* param, Real a_dt)
Real Projector::getScale(Real a_scale, Real a_dt)
{


  Real newScale = a_scale;

  if (newScale < 0)
  {
    newScale = -newScale*a_dt;
  }

  return newScale;
}

// This is used for projecting a time dependent advection velocity
// --------------------------------------------------------------
int Projector::levelMacProject(LevelData<FluxBox>& a_uEdge,
                                  Real a_oldTime, Real a_dt,
                                  const RefCountedPtr<LevelData<FArrayBox> > a_porosityPtr,
                                  const RefCountedPtr<LevelData<FArrayBox> > a_crsePorosityPtr,
                                  const RefCountedPtr<LevelData<FluxBox> > a_porosityEdgePtr,
                                  const RefCountedPtr<LevelData<FluxBox> > a_crsePorosityEdgePtr)
{

  if (s_verbosity > 3)
  {
    pout() << "CCProjector::levelMacProject" << endl;
  }

  Divergence::levelDivergenceMAC(MACrhs(), a_uEdge, m_dx);

  LevelData<FArrayBox>* pressureBCPtr = NULL;

  // NB this used to be 1
  Real phiScale = getScale(m_phiScale, a_dt); //getPhiScale(a_dt);
  Real CFscale = getScale(m_MACbcScale, a_dt);

  bool usePhiBCs = !m_usePiAdvectionBCs;

  if (m_crseProjPtr != NULL)
  {
    // coarse-fine BC is 0.5*dt*(coarse Pi)
    const DisjointBoxLayout crseGrids = m_crseProjPtr->getBoxes();
    pressureBCPtr = new LevelData<FArrayBox>(crseGrids,1);
    LevelData<FArrayBox>& crseBC = *pressureBCPtr;

    Interval comps(0,0);
    if (usePhiBCs)
    {
      const LevelData<FArrayBox>& crsePhi = m_crseProjPtr->phi();
      crsePhi.copyTo(comps,crseBC,comps);
    }
    else
    {
      const LevelData<FArrayBox>& crsePi = m_crseProjPtr->Pi();
      crsePi.copyTo(comps,crseBC,comps);
    }

    DataIterator dit = crseBC.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      crseBC[dit] *= CFscale;

    } // end loop over grids

  } // end if crse porosity

  // report sum(rhs)
  if (s_verbosity >= 2)
  {
    DisjointBoxLayout* finerGridsPtr = NULL;
    int nRefFine = -1;
    Real sumRHS = computeSum(MACrhs(), finerGridsPtr,
                             nRefFine, m_dx, MACrhs().interval());
    pout() << "    MAC projection (level " << m_level << ") -- sum(RHS) = " << sumRHS << endl;
  }

  for (DataIterator dit = MACrhs().dataIterator(); dit.ok(); ++dit)
  {
    MACrhs()[dit].divide(phiScale);
  }

  MACrhs().exchange();


  // now solve for phi
  int exitStatus = solveMGlevel(m_phi, pressureBCPtr, MACrhs(),
               a_porosityPtr, a_crsePorosityPtr,
               a_porosityEdgePtr, a_crsePorosityEdgePtr);

  // now correct

  // first re-apply physical and copy BC's
  Interval phiComps(0,0);
  m_phi.exchange(phiComps);
  //  BCHolder bcHolder = m_physBCPtr->gradMacPressureFuncBC();
  BCHolder bcHolder = m_physBCPtr->BasicPressureFuncBC(false);
  const DisjointBoxLayout& levelGrids = getBoxes();

  DataIterator dit = m_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    bcHolder.operator()(m_phi[dit],
                        levelGrids[dit],
                        m_domain,
                        m_dx,
                        false); // not homogeneous
  }

#define CORNERCOPIER
#ifdef CORNERCOPIER
  CornerCopier cornerExchangeCopier(m_phi.getBoxes(),
                                    m_phi.getBoxes(),
                                    m_domain,
                                    m_phi.ghostVect(),
                                    true);

  m_phi.exchange(m_phi.interval(), cornerExchangeCopier);
#endif

  applyMacCorrection(a_uEdge, pressureBCPtr, phiScale);

  if (pressureBCPtr != NULL && phiScale != a_dt)
  {
    delete pressureBCPtr;
    pressureBCPtr = NULL;
  }
  // that should be it!


  // Apply filter. Will need to make this AMR safe
  //  u = u-alpha*grad(div(u))
  Real alpha = 0.0;
  ParmParse pp("projection");
  pp.query("filter", alpha);

  if (abs(alpha) > 0.0)
  {
    LevelData<FluxBox> correction(a_uEdge.disjointBoxLayout(), 1);
    LevelData<FArrayBox> divU(MACrhs().disjointBoxLayout(), 1, IntVect::Unit);

    for (dit.reset(); dit.ok(); ++dit)
    {
      divU[dit].setVal(0.0);
    }
    Divergence::levelDivergenceMAC(divU, a_uEdge, m_dx);
    divU.exchange();

    Real maxDivU = ::computeNorm(divU, NULL, 1, m_dx, Interval(0,0), 0);
    Real maxVel = ::norm(a_uEdge, Interval(0,0), 0);

    // Don't correct if div(u) is very small, or of order U
    if (maxDivU > 1e-10 && maxVel/maxDivU > 100)
    {
      pout() << "  CCProjector - smoothing div(u) " << endl;

      Gradient::levelGradientMAC(correction, divU, m_dx);

      correction.exchange();

      // Correct
      for (dit.reset(); dit.ok(); ++dit)
      {
        FluxBox& thisCorr = correction[dit];
        FluxBox& thisEdgeVel = a_uEdge[dit];

        for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisCorrDir = thisCorr[dir];
          FArrayBox& thisVelDir = thisEdgeVel[dir];

          //        if (m_porosityEdgePtr != NULL)
          //        {
          //          thisGradDir.mult((*m_porosityEdgePtr)[dit][dir],thisGradDir.box(), 0, 0);
          //        }

          thisCorrDir.mult(alpha);

          thisVelDir -= thisCorrDir;
        }
      }

    }
  }

  checkDivergence(a_uEdge);

  return exitStatus;
}

// --------------------------------------------------------------
void
Projector::applyMacCorrection(LevelData<FluxBox>& a_uEdge,
                                     LevelData<FArrayBox>* crseBCDataPtr,
                                     Real phiScale)
{
  // this function computes a_uEdge - chi*grad(phi)
  // (assumes all BC's already set)

  const DisjointBoxLayout& levelGrids = a_uEdge.getBoxes();
  LevelData<FluxBox> gradPhi(levelGrids,1);

  // compute gradient
  Gradient::levelGradientMAC(gradPhi,m_phi, crseBCDataPtr, m_dx,
                             m_gradIVS,
                             m_cfInterp);

  // now rescale and subtract from uEdge
  DataIterator dit = a_uEdge.dataIterator();

  // assumes that uEdge and phi can use same dataIterator
  for (dit.reset(); dit.ok(); ++dit)
  {
    FluxBox& thisGrad = gradPhi[dit];
    FluxBox& thisEdgeVel = a_uEdge[dit];

    for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& thisGradDir = thisGrad[dir];
      FArrayBox& thisVelDir = thisEdgeVel[dir];

      if (m_porosityEdgePtr != NULL)
      {
        thisGradDir.mult((*m_porosityEdgePtr)[dit][dir],thisGradDir.box(), 0, 0);
      }

      thisGradDir.mult(phiScale);

      thisVelDir -= thisGradDir;
    }
  }

  gradPhi.copyTo(MACcorrection());

  //  if (crseBCDataPtr!=NULL)
  //  {
  //    delete crseBCDataPtr;
  //    crseBCDataPtr = NULL;
  //  }
}

void Projector::setPressureScalePtr(RefCountedPtr<LevelData<FArrayBox> > a_pressureScalePtr)
{
  m_porosityPtr = a_pressureScalePtr;
}

void Projector::setPressureScaleEdgePtr(RefCountedPtr<LevelData<FluxBox> > a_pressureScaleEdgePtr)
{
  m_porosityEdgePtr = a_pressureScaleEdgePtr;
}

void Projector::getCFBC(LevelData<FArrayBox>& velBC, LevelData<FArrayBox>* a_crseVelPtr,
                          LevelData<FArrayBox>* a_crsePorosityPtr, Real a_dt)
{

  bool addSubtractGradP = true;
  ParmParse pp("main");
  pp.query("addSubtractGradP", addSubtractGradP);

  if (doQuadInterp())
  {
    CH_assert(a_crseVelPtr != NULL);
    CH_assert(m_crseProjPtr != NULL);

    // will need to add dt*gradPi from crse vel.
    const DisjointBoxLayout& crseLevelGrids = m_crseProjPtr->getBoxes();
    velBC.define(crseLevelGrids, SpaceDim);
    //              velBCPtr = new LevelData<FArrayBox>(crseLevelGrids, SpaceDim);
    //              LevelData<FArrayBox>& crseData = *velBCPtr;
    LevelData<FArrayBox>& crseData = velBC;
    m_crseProjPtr->gradPi(crseData);

    // now multiply by dt and add crseVel
    DataIterator ditCrse = crseData.dataIterator();
    LevelData<FArrayBox>& crseVel = *a_crseVelPtr;

    // inherent assumption that crseVel and crseGradPi share
    // the same layout here
    for (ditCrse.reset(); ditCrse.ok(); ++ditCrse)
    {
      FArrayBox& thisCrseVel = crseVel[ditCrse];
      FArrayBox& thisCrseGradPi = crseData[ditCrse];

      if (addSubtractGradP)
      {
        // note that we multiply by dt on _this_ level
        thisCrseGradPi *= a_dt;
      }
      else
      {
        // if we're not doing the whole add/subtract grad(P) thing, then
        // we need to use dt on the crse level so we remove the entire grad(P) from the crse vel
        thisCrseGradPi *= 0;
      }

      if (a_crsePorosityPtr!=NULL)
      {

        for (int dir=0; dir<SpaceDim; dir++)
        {
          thisCrseGradPi.mult((*a_crsePorosityPtr)[ditCrse], 0, dir);
        }

      }

      // now add crseVel
      thisCrseGradPi.plus(thisCrseVel);
    }
  }
  else
  {
    // alternative is not implemented at this point
    MayDay::Error("non-quadratic interpolation BC's not implemented");
  }

}

void Projector::scaleRHS(LevelData<FArrayBox>& a_rhs, Real a_scale)
{
  DataIterator dit = a_rhs.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    a_rhs[dit] *= a_scale;
  }

}

// ---------------------------------------------------------------
void Projector::LevelProject(LevelData<FArrayBox>& a_velocity,
                               LevelData<FArrayBox>* a_crseVelPtr,
                               const Real a_newTime, const Real a_dt,
                               const RefCountedPtr<LevelData<FArrayBox> > a_porosityPtr,
                               const RefCountedPtr<LevelData<FArrayBox> > a_crsePorosityPtr,
                               const RefCountedPtr<LevelData<FluxBox> > a_porosityEdgePtr,
                               const RefCountedPtr<LevelData<FluxBox> > a_crsePorosityEdgePtr,
                               const bool a_isViscous)
{
  if (s_verbosity >= 5)
  {
    pout() << "CCProjector::LevelProject "            << endl;
  }

  CH_TIME("CCProjector::levelProject");

  LevelData<FArrayBox>* velBCPtr=NULL;
  LevelData<FArrayBox> velBC;
  ParmParse pp("projection");

  // just to be safe.  proper place for this may be outside this function
  Interval velComps(0,SpaceDim-1);
  a_velocity.exchange(velComps);

  // Get CF BC if we need it
  QuadCFInterp velCFInterp;

  if (m_level > 0)
  {
    getCFBC(velBC, a_crseVelPtr, a_crsePorosityPtr, a_dt);

    velBCPtr = &velBC;

    if (doQuadInterp() && m_crseProjPtr != NULL)
    {
      // define two-component CF-interp object for velocities
      velCFInterp.define(a_velocity.getBoxes(),
                         &(velBCPtr->getBoxes()),
                         m_dx, m_nRefCrse, SpaceDim, m_domain);
    }
  }

  if (s_verbosity >= 5)
  {
    pout() << "CCProjector::LevelProject - calculateDivergence "            << endl;
  }

  // now compute RHS for projection
  // BCs should already be set
  Divergence::levelDivergenceCC(CCrhs(), a_velocity, velBCPtr,
                                m_dx, doQuadInterp(), velCFInterp);

  // for proper scaling of pi, divide this by dt
  // so solving lap(pi) = div(u)/dt
  Real dtScale = 1.0/a_dt;
  bool applyScaling = true;
  pp.query("scaleCCRHS", applyScaling);
  if (!applyScaling)
  {
    dtScale = 1.0;
  }
  scaleRHS(CCrhs(), dtScale);

  // set up coarse BC's for solve, then solve
  const LevelData<FArrayBox>* crsePiPtr = NULL;
  if (m_crseProjPtr != NULL) crsePiPtr = &(m_crseProjPtr->Pi());

  // report sum(rhs)
  if (s_verbosity >= 2)
  {
    DisjointBoxLayout* finerGridsPtr = NULL;
    int nRefFine = -1;
    Real sumRHS = computeSum(CCrhs(), finerGridsPtr,
                             nRefFine, m_dx, CCrhs().interval());
    pout() << "    Level projection -- sum(RHS) = " << sumRHS << endl;
  }

  if (s_verbosity >= 5)
  {
    pout() << "CCProjector::LevelProject - solve for pressure correction "            << endl;
  }


  // now solve for pi
  solveMGlevel(m_Pi, crsePiPtr, CCrhs(),
               a_porosityPtr, a_crsePorosityPtr,
               a_porosityEdgePtr, a_crsePorosityEdgePtr,
               true); // cell centred

  // apply appropriate physical BC's
  BCHolder bcHolder = m_physBCPtr->gradPiFuncBC(); // this is what CC project used to use
  //  BCHolder bcHolder = m_physBCPtr->BasicPressureFuncBC(false); // this is what MAC project uses

  const DisjointBoxLayout& levelGrids = getBoxes();

  if (s_verbosity >= 5)
  {
    pout() << "CCProjector::LevelProject - apply BCs "            << endl;
  }

  // loop over grids...
  DataIterator dit = CCrhs().dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    bcHolder.operator()(m_Pi[dit],
                        levelGrids[dit],
                        m_domain,
                        m_dx,
                        false); // not homogeneous
  }

  m_Pi.exchange();

  if (m_crseProjPtr != NULL)
  {
    // reapply coarse-fine BC's here if necessary
    m_cfInterp.coarseFineInterp(m_Pi, *crsePiPtr);
    // also clean up after ourselves!

    // Don't need to do this anymore
    //    delete velBCPtr;
    //    velBCPtr = NULL;
  }

  // dt scale was 1/dt, so correct scale is dt
  Real correctScale = 1/dtScale; //1/dtScale;
  applyCCcorrection(a_velocity, &m_Pi, a_porosityPtr, correctScale);

  // check resulting velocity field
  // to do this, need to reset physical BC's

  // Smoothing
  Real alpha = 0.0;

  pp.query("CCfilter", alpha);

  if (abs(alpha) > 0.0)
  {
    LevelData<FArrayBox> correction(a_velocity.disjointBoxLayout(), SpaceDim);
    LevelData<FArrayBox> divU(CCrhs().disjointBoxLayout(), 1, IntVect::Unit);

    for (dit.reset(); dit.ok(); ++dit)
    {
      divU[dit].setVal(0.0);
    }
    Divergence::levelDivergenceCC(divU, a_velocity, velBCPtr,
                                  m_dx, doQuadInterp(), velCFInterp);

    divU.exchange();

    pout() << "  CCProjector - smoothing div(u) " << endl;

    Gradient::levelGradientCC(correction, divU, m_dx);

    correction.exchange();

    // Correct
    for (dit.reset(); dit.ok(); ++dit)
    {
      FArrayBox& thisCorr = correction[dit];
      FArrayBox& thisVel = a_velocity[dit];

      thisCorr.mult(alpha*m_dx*m_dx);

      thisVel -= thisCorr;

    }

  }

}

void Projector::AdditionalLevelProject(LevelData<FArrayBox>& a_velocity,
                               LevelData<FArrayBox>* a_crseVelPtr,
                               const Real a_newTime, const Real a_dt,
                               const bool a_isViscous)
{
  if (s_verbosity >= 5)
  {
    pout() << "CCProjector::LevelProject "            << endl;
  }

  CH_TIME("CCProjector::levelProject");

  LevelData<FArrayBox>* velBCPtr=NULL;
  LevelData<FArrayBox> velBC;
  ParmParse pp("projection");

  // just to be safe.  proper place for this may be outside this function
  Interval velComps(0,SpaceDim-1);
  a_velocity.exchange(velComps);

  LevelData<FArrayBox> rhs(a_velocity.disjointBoxLayout(), 1);
  LevelData<FArrayBox> pressure(a_velocity.disjointBoxLayout(), 1, IntVect::Unit);

  // Get CF BC if we need it
  QuadCFInterp velCFInterp;

  if (m_level > 0)
  {
//    getCFBC(velBC, a_crseVelPtr, a_crsePorosityPtr, a_dt);
//
//    velBCPtr = &velBC;
//
//    if (doQuadInterp() && m_crseProjPtr != NULL)
//    {
//      // define two-component CF-interp object for velocities
//      velCFInterp.define(a_velocity.getBoxes(),
//                         &(velBCPtr->getBoxes()),
//                         m_dx, m_nRefCrse, SpaceDim, m_domain);
//    }
  }

  if (s_verbosity >= 5)
  {
    pout() << "CCProjector::LevelProject - calculateDivergence "            << endl;
  }

  // now compute RHS for projection
  // BCs should already be set
  Divergence::levelDivergenceCC(rhs, a_velocity, velBCPtr,
                                m_dx, doQuadInterp(), velCFInterp);

  // for proper scaling of pi, divide this by dt
  // so solving lap(pi) = div(u)/dt
  Real dtScale = 1.0/a_dt; //1.0/a_dt;
  bool applyScaling = true;
  pp.query("scaleCCRHS", applyScaling);
  if (!applyScaling)
  {
    dtScale = 1.0;
  }
  scaleRHS(rhs, dtScale);

  // this won't work with AMR yet
  CH_assert(m_level==0);

  // set up coarse BC's for solve, then solve
  const LevelData<FArrayBox>* crsePiPtr = NULL;
  if (m_crseProjPtr != NULL) crsePiPtr = &(m_crseProjPtr->Pi());

  // report sum(rhs)
  if (s_verbosity >= 2)
  {
    DisjointBoxLayout* finerGridsPtr = NULL;
    int nRefFine = -1;
    Real sumRHS = computeSum(CCrhs(), finerGridsPtr,
                             nRefFine, m_dx, CCrhs().interval());
    pout() << "    Level projection -- sum(RHS) = " << sumRHS << endl;
  }

  if (s_verbosity >= 5)
  {
    pout() << "CCProjector::LevelProject - solve for pressure correction "            << endl;
  }


  // now solve for pi
  solveMGlevel(pressure, crsePiPtr, rhs, true); // cell centred

  // apply appropriate physical BC's
  BCHolder bcHolder = m_physBCPtr->gradPiFuncBC(); // this is what CC project used to use
  //  BCHolder bcHolder = m_physBCPtr->BasicPressureFuncBC(false); // this is what MAC project uses

  const DisjointBoxLayout& levelGrids = getBoxes();

  if (s_verbosity >= 5)
  {
    pout() << "CCProjector::LevelProject - apply BCs "            << endl;
  }

  // loop over grids...
  DataIterator dit = CCrhs().dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    bcHolder.operator()(pressure[dit],
                        levelGrids[dit],
                        m_domain,
                        m_dx,
                        false); // not homogeneous
  }

  pressure.exchange();

  if (m_crseProjPtr != NULL)
  {
    // reapply coarse-fine BC's here if necessary
    m_cfInterp.coarseFineInterp(pressure, *crsePiPtr);
    // also clean up after ourselves!

    // Don't need to do this anymore
    //    delete velBCPtr;
    //    velBCPtr = NULL;
  }

  // dt scale was 1/dt, so correct scale is dt
  Real correctScale = 1/dtScale; //1/dtScale;
  applyCCcorrection(a_velocity, &pressure, NULL, correctScale);



}

// ---------------------------------------------------------------
void Projector::applyCCcorrection(LevelData<FArrayBox>& a_velocity,
                                  LevelData<FArrayBox>* a_pressure,
                                  LevelData<FArrayBox>* a_pressureScalePtr,
                                    const Real scale) //const
{
  if (s_verbosity >= 3)
  {
    pout() << "CCProjector::applyCCcorrection with scale " << scale             << endl;
  }

  const DisjointBoxLayout& grids = getBoxes();
  LevelData<FArrayBox> gradPi(grids, SpaceDim);

  // assumes that all relevant BC's have already been set
  // trying higher order gradient

  bool HO_corr = false;
  ParmParse ppProjection("projection");
  ppProjection.query("HO_CC_Corr", HO_corr);
  if (HO_corr)
  {
  Gradient::levelGradientCC_HO(gradPi, *a_pressure, m_dx);
  }
  else
  {
    // Default
      Gradient::levelGradientCC(gradPi, *a_pressure, m_dx);
  }
  //  this->gradPi(gradPi);

  // vel = vel - scale*pressureScale*gradPi
  DataIterator dit = gradPi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    FArrayBox& thisGradPi = gradPi[dit];
    FArrayBox& thisVel = a_velocity[dit];

    if (a_pressureScalePtr!=NULL)
    {
      for (int dir=0; dir<SpaceDim; dir++)
      {
        thisGradPi.mult((*a_pressureScalePtr)[dit],0,dir);
      }
    }

    thisGradPi *= scale;
    thisVel -= thisGradPi;
  }

  gradPi.copyTo(m_CCcorrection);

}

// ---------------------------------------------------------------
void Projector::doSyncOperations(Vector<LevelData<FArrayBox>* >& a_velocity,
                                   Vector<LevelData<FArrayBox>* >& a_lambda,
                                   Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosityFace,
                                   Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity,
                                   const Real a_newTime, const Real a_dtLevel)
{
  CH_TIME("CCProjector::doSyncOperations");

  if (s_verbosity > 3)
  {
    pout() << "CCProjector::doSyncOperations" << endl;
  }



  AMRMultiGrid<LevelData<FArrayBox> >* dspSolver = new
      AMRMultiGrid<LevelData<FArrayBox> >;

  defineMultiGrid(*dspSolver, a_velocity, a_porosityFace, false);

  // now call sync projection
  doSyncProjection(a_velocity, a_porosity, a_newTime, a_dtLevel, *dspSolver);

  // Need to delete dspSolver->m_bottomSolver from AMRMultiGrid.
  delete dspSolver;

  bool vdBCs = true;
  ParmParse ppProj("projection");
  ppProj.query("vd_sync_bcs", vdBCs);

  AMRMultiGrid<LevelData<FArrayBox> >* cvdcSolver = new
      AMRMultiGrid<LevelData<FArrayBox> >;
  defineMultiGrid(*cvdcSolver, a_velocity, a_porosityFace,
                  vdBCs); // surely this should be true!??!

  // now do freestream preservation solve
  computeVDCorrection(a_lambda, a_porosityFace, a_newTime, a_dtLevel, *cvdcSolver);

  Projector* proj = this;

  while(proj)
  {
    proj = proj->fineProjPtr();

  }


  // Need to delete cvdcSolver->m_bottomSolver from AMRMultiGrid.
  delete cvdcSolver;

}


// ---------------------------------------------------------------
void Projector::doPostRegridOps(Vector<LevelData<FArrayBox>* >& a_lambda,
                                  Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosity,
                                  const Real a_dt, const Real a_time, const Real a_etaScale)
{
  //


  AMRMultiGrid<LevelData<FArrayBox> >* bigSolverPtr = new
      AMRMultiGrid<LevelData<FArrayBox> >;

  defineMultiGrid(*bigSolverPtr, a_lambda, a_porosity, true); // true - freestream solve

  // for inviscid flow, only do this
  // now do freestream preservation solve
  m_etaLambda *= a_etaScale;
  computeVDCorrection(a_lambda, a_porosity, a_time, a_dt, *bigSolverPtr);
  m_etaLambda /= a_etaScale;

  // Need to delete bigSolverPtr->m_bottomSolver from AMRMultiGrid.
  delete bigSolverPtr;


}

// ---------------------------------------------------------------
void Projector::doSyncProjection(Vector<LevelData<FArrayBox>* >& a_velocity,
                                   Vector<RefCountedPtr< LevelData<FArrayBox> > >& a_porosity,
                                   const Real a_newTime, const Real a_dtSync,
                                   AMRMultiGrid<LevelData<FArrayBox> >& a_solver
)
{
  if (doSyncProjection())
  {
    // set up temp storage
    int vectorSize = a_velocity.size();

    Vector<LevelData<FArrayBox>* > syncRHS(vectorSize,NULL);
    Vector<LevelData<FArrayBox>* > syncCorr(vectorSize,NULL);

    // initially
    Projector* levelProjPtr = this;
    //    while (!levelProjPtr->isFinestLevel())
    //    {
    //      levelProjPtr = levelProjPtr->fineProjPtr();
    //    }

    // Need to be consistent about this
    // Can't have syncRHS with more/less elements
    // than the finest level that we pass to the solver
    int finestLevel = vectorSize - 1;

    // reset projection ptr to this level
    levelProjPtr = this;

    // this is a bogus leveldata pointer for the composite divergences
    LevelData<FArrayBox>* crseVelPtr=NULL;
    LevelData<FArrayBox>* fineVelPtr = NULL;

    // loop over levels to allocate storage and compute RHS
    for (int lev = m_level; lev<= finestLevel; lev++)
    {
      const DisjointBoxLayout& levelGrids = a_velocity[lev]->getBoxes();

      syncRHS[lev]= new LevelData<FArrayBox>(levelGrids,1);
      syncCorr[lev] = &(levelProjPtr->eSync());

      Real dx = levelProjPtr->dx();
      const ProblemDomain& levelDomain = levelProjPtr->dProblem();
      int nRefCrse = levelProjPtr->nRefCrse();
      int nRefFine = -1;

      if (lev > 0) crseVelPtr = a_velocity[lev-1];
      if (lev < finestLevel)
      {
        fineVelPtr = a_velocity[lev+1];
        nRefFine = levelProjPtr->fineProjPtr()->nRefCrse();
      }
      else
      {
        fineVelPtr = NULL;
      }

      // do exchange here
      a_velocity[lev]->exchange(a_velocity[lev]->interval());

      // now compute divergence
      Divergence::compDivergenceCC(*syncRHS[lev],*a_velocity[lev],
                                   crseVelPtr, fineVelPtr, dx,
                                   nRefCrse, nRefFine, levelDomain,
                                   doQuadInterp());

      // increment projection level
      levelProjPtr = levelProjPtr->fineProjPtr();

    } // end loop over levels for computation of RHS

    // take care of C/F Boundary condition for lBase
    if (m_level > 0)
    {
      syncCorr[m_level-1] = &(m_crseProjPtr->eSync());
      // need to rescale crse-level BC
      LevelData<FArrayBox>& crseSyncCorr = *(syncCorr[m_level-1]);
      DataIterator ditCrse = crseSyncCorr.dataIterator();
      for (ditCrse.reset(); ditCrse.ok(); ++ditCrse)
      {
        crseSyncCorr[ditCrse] *= a_dtSync;
      }
    }

    // debugging tool -- want sum(rhs) = 0
    Vector<int> nRefFineVect(syncRHS.size());
    levelProjPtr = this;
    for (int lev=m_level; lev<=finestLevel; lev++)
    {
      if (lev < finestLevel)
      {
        nRefFineVect[lev] = levelProjPtr->fineProjPtr()->nRefCrse();
      }
      else
      {
        nRefFineVect[lev] = 1;
      }
      levelProjPtr = levelProjPtr->fineProjPtr();
    }

    Interval sumComps(0,0);
    Real sumRHS;
    sumRHS =  computeSum(syncRHS, nRefFineVect, m_dx, sumComps, m_level);
    if (s_verbosity >= 3)
    {
      pout() << "  Sum(RHS) for sync solve = "
          << setiosflags(ios::scientific) << sumRHS << endl;

    }

    // now solve
    a_solver.solve(syncCorr, syncRHS, vectorSize-1, m_level,
                   true); // initialize syncCorr to zero

    // apply physical boundary conditions before taking gradients
    levelProjPtr = this;

    BCHolder bcHolder = m_physBCPtr->gradESyncFuncBC();

    for (int lev=m_level; lev<=finestLevel; lev++)
    {

      // loop over grids
      LevelData<FArrayBox>& level_eSync = *syncCorr[lev];
      const ProblemDomain& levelDomain = levelProjPtr->dProblem();
      Real levelDx = levelProjPtr->dx();
      const DisjointBoxLayout& levelGrids = levelProjPtr->getBoxes();
      DataIterator ditLev = level_eSync.dataIterator();

      for (ditLev.reset(); ditLev.ok(); ++ditLev)
      {
        bcHolder.operator()(level_eSync[ditLev],
                            levelGrids[ditLev],
                            levelDomain,
                            levelDx,
                            false); // not homogeneous
      }
      // exchange here should be unnecessary
      //level_eSync.exchange(corrComps);
      levelProjPtr = levelProjPtr->fineProjPtr();
    }

    // now apply sync correction
    Real scale = 1.0;
    LevelData<FArrayBox>* crseCorrPtr = NULL;
    if (m_level > 0) crseCorrPtr = syncCorr[m_level-1];

    applySyncCorrection(a_velocity, a_porosity, scale, crseCorrPtr);

    // cleaning up after we're done --
    // now need to rescale eSync to look like other pressures
    // start with level-1 because we rescaled it for the BC
    // also need to delete syncRHS

    Real corrScale = 1.0/a_dtSync;
    for (int lev=m_level-1; lev<=finestLevel; lev++)
    {
      if (lev >=0)
      {
        LevelData<FArrayBox>& levelCorr = *syncCorr[lev];
        DataIterator dit = levelCorr.dataIterator();
        for (dit.reset(); dit.ok(); ++dit)
        {
          levelCorr[dit] *= corrScale;
        }
      } // end if level > 0
    } // end loop over levels for rescaling eSync

    // loop over levels to delete syncRHS
    for (int lev= m_level; lev<=finestLevel; lev++)
    {
      delete syncRHS[lev];
      syncRHS[lev] = NULL;
    }
  } // end if do sync in the first place
}

// ---------------------------------------------------------------
void Projector::applySyncCorrection(Vector<LevelData<FArrayBox>* >& a_velocities,
                                      Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity,
                                      const Real a_scale,
                                      LevelData<FArrayBox>* crseCorrPtr)
{
  if (s_verbosity >= 3)
  {
    pout() << "CCProjector::applySyncCorrection " << endl;
  }

  if (m_applySyncCorrection)
  {
    const DisjointBoxLayout& levelGrids = getBoxes();

    LevelData<FArrayBox> levelGradCorr(levelGrids, SpaceDim);
    LevelData<FArrayBox>* fineCorrPtr = NULL;
    int nRefFine = -1;


    if (!isFinestLevel())
    {
      fineCorrPtr = &(m_fineProjPtr->eSync());
      nRefFine = m_fineProjPtr->nRefCrse();
    }

    Gradient::compGradientCC(levelGradCorr, m_eSync,
                             crseCorrPtr, fineCorrPtr,
                             m_dx, nRefFine, m_domain,
                             m_gradIVS,
                             m_cfInterp);

    // now loop over grids and correct
    LevelData<FArrayBox>& levelVel = *a_velocities[m_level];
    LevelData<FArrayBox>& porosity = *a_porosity[m_level];

    DataIterator dit = levelVel.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      levelGradCorr[dit] *= a_scale;

      if (m_scaleSyncCorrection)
      {
        for (int dir=0; dir<SpaceDim; dir++)
        {
          levelGradCorr[dit].mult(porosity[dit], 0, dir, 1);
        }
      }

      levelVel[dit] -= levelGradCorr[dit];
    }

    // now recursively call finer level
    if (!isFinestLevel())
    {
      m_fineProjPtr->applySyncCorrection(a_velocities, a_porosity, a_scale,
                                         &m_eSync);

      // average down finer level to this level
      const DisjointBoxLayout& fineGrids = m_fineProjPtr->getBoxes();
      int nRefFine = m_fineProjPtr->nRefCrse();
      CoarseAverage avgDownThingy(fineGrids, SpaceDim, nRefFine);

      avgDownThingy.averageToCoarse(*a_velocities[m_level],
                                    *a_velocities[m_level+1]);
    }
  }
}


void Projector::computeVDCorrection(Vector<LevelData<FArrayBox>* >& a_lambda,
                                      Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosity,
                                      const Real a_newTime, const Real a_dtSync)
{

  AMRMultiGrid<LevelData<FArrayBox> >* cvdcSolver = new
      AMRMultiGrid<LevelData<FArrayBox> >;
  defineMultiGrid(*cvdcSolver, a_lambda, a_porosity,
                  true);

  // now do freestream preservation solve
  computeVDCorrection(a_lambda, a_porosity, a_newTime, a_dtSync, *cvdcSolver);


  // Need to delete cvdcSolver->m_bottomSolver from AMRMultiGrid.
  delete cvdcSolver;

}

// ---------------------------------------------------------------
void Projector::computeVDCorrection(Vector<LevelData<FArrayBox>* >& a_lambda,
                                      Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosity,
                                      const Real a_newTime, const Real a_dtSync,
                                      AMRMultiGrid<LevelData<FArrayBox> >& a_solver
)
{
  if (s_verbosity > 3)
  {
    pout() << "CCProjector::computeVDCorrection, eta = " << m_etaLambda << endl;
  }

  // only do this if eta != 0
  if (doMacSync())
  {

    if (m_level == 0) s_lambda_timestep = a_dtSync;
    // also need to handle initial case where level 0
    // hasn't been hit yet. in this case, back out
    // what the level 0 timestep is.
    if (s_lambda_timestep ==0.0)
    {
      int timestepMult = 1;
      Projector* levelProjPtr = this;
      for (int lev=m_level; lev>0; lev--)
      {
        timestepMult *= levelProjPtr->nRefCrse();
        levelProjPtr = levelProjPtr->crseProjPtr();
      }
      s_lambda_timestep = a_dtSync*timestepMult;
    }

    // set up temp storage
    int vectorSize = a_lambda.size();

    Vector<LevelData<FArrayBox>* > VDCorr(vectorSize);

    Projector* levelProjPtr = this;
    while (!levelProjPtr->isFinestLevel())
    {
      levelProjPtr = levelProjPtr->fineProjPtr();
    }

    int finestLevel = levelProjPtr->getLevel();

    // reset projection ptr to this level
    levelProjPtr = this;

    CH_assert (a_dtSync != 0);
    Real lambdaMult;
    if (s_constantLambdaScaling)
    {
      lambdaMult = m_etaLambda/s_lambda_timestep;
    }
    else
    {
      lambdaMult = m_etaLambda/a_dtSync;
    }

    Vector<LevelData<FArrayBox>* > newLambda(a_lambda.size());

    Vector<int> nRefFineVect(a_lambda.size());
    levelProjPtr = this;
    for (int lev=m_level; lev<=finestLevel; lev++)
    {
      if (lev < finestLevel)
      {
        nRefFineVect[lev] = levelProjPtr->fineProjPtr()->nRefCrse();
      }
      else
      {
        nRefFineVect[lev] = 1;
      }
      levelProjPtr = levelProjPtr->fineProjPtr();
    }

    // Compute lambda-1
    levelProjPtr = this;
    for (int lev = m_level; lev<vectorSize; lev++)
    {

      LevelData<FArrayBox>& oldLambda = *(a_lambda[lev]);

      DataIterator dit = oldLambda.dataIterator();

      // RHS = (lambda-1)
      for (dit.reset(); dit.ok(); ++dit)
      {
        if (m_scale_lambda_with_porosity)
        {
          // Compute 1/porosity
          // Note that inversePorosity has one component in each direction
          // so that EdgeToCell works. We only end up using the x-direction,
          // but this shouldn't make too much difference

          FArrayBox inversePorosity(oldLambda[dit].box(), SpaceDim);
          //          inversePorosity.copy((*a_porosity[lev])[dit]);
          EdgeToCell((*a_porosity[lev])[dit], inversePorosity);
          inversePorosity.invert(1);

          oldLambda[dit] -= inversePorosity;
        }
        else
        {
          oldLambda[dit] -= 1.0;
        }


      }
    }


    //    Real maxLambda = ::computeNorm(a_lambda, nRefFineVect, m_dx, Interval(0,0), 0, m_level);

    levelProjPtr = this;
    for (int lev = m_level; lev<vectorSize; lev++)
    {
      // stuff eLambda into VDCorr array
      VDCorr[lev] = &(levelProjPtr->eLambda());

      // use lambda array as RHS for solve
      //      LevelData<FArrayBox>& thisLambda = *(a_lambda[lev]);

      //      // Make new array for RHS so we can do (lambda-1)^3
      newLambda[lev] = new LevelData<FArrayBox>(a_lambda[lev]->disjointBoxLayout(), 1, a_lambda[lev]->ghostVect());
      //      a_lambda[lev]->copyTo(*newLambda[lev]);
      LevelData<FArrayBox>& thisLambda = *(newLambda[lev]);
      LevelData<FArrayBox>& oldLambda = *(a_lambda[lev]);

      DataIterator dit = thisLambda.dataIterator();


      // RHS = (lambda-1)*eta/dtSync
      for (dit.reset(); dit.ok(); ++dit)
      {
        thisLambda[dit].copy(oldLambda[dit]);

        thisLambda[dit] *= lambdaMult;

        if (m_scale_lambda_err_with_porosity)
        {
          FArrayBox ccPorosity(thisLambda[dit].box(), SpaceDim);
          EdgeToCell((*a_porosity[lev])[dit], ccPorosity);
          thisLambda[dit].divide(ccPorosity);
        }

        levelProjPtr->eLambda()[dit].setVal(0.0);
      }

      // advance projection pointer to next finer level
      levelProjPtr = levelProjPtr->fineProjPtr();
    } // end loop over levels

    // if necessary, define coarse BC here
    if (m_level > 0)
    {
      VDCorr[m_level-1] = &(m_crseProjPtr->eLambda());
      // rescale crse-level BC
      // try setting coarse-level BC to 0
      setValLevel(m_crseProjPtr->eLambda(), 0.0);
    }

    // debugging tool -- want sum(rhs) = 0

    Interval sumComps(0,0);
    m_sumVDrhs =  computeSum(a_lambda, nRefFineVect, m_dx, sumComps, m_level);
    m_sumVDrhs = m_sumVDrhs/m_etaLambda;
    pout() << "  Sum(RHS) for VD solve on level " << m_level << " (divided by eta) = "
        << setiosflags(ios::scientific) << m_sumVDrhs << " (current time " <<a_newTime << ")" <<  endl;
    // now solve

    a_solver.solve(VDCorr, newLambda,
                   vectorSize-1,  // max
                   m_level, // base
                   true); // initialize VDCorr to zero

    // apply all bc's here

    // Unsure why we were using extrap BCs here before. Using the freestream BCs instead.
    // If we don't use freestream BCs, we don't conserve fluxes with neumann BCs...
    //            BCHolder bcHolder = m_physBCPtr->gradELambdaFuncBC();
    BCHolder bcHolder = m_physBCPtr->FreestreamCorrFuncBC();

    Interval corrComps(0,0);
    levelProjPtr = this;
    for (int lev=m_level; lev<=finestLevel; lev++)
    {
      LevelData<FArrayBox>& levelCorr = *(VDCorr[lev]);
      const ProblemDomain& levelDomain = levelProjPtr->dProblem();
      Real levelDx = levelProjPtr->dx();
      const DisjointBoxLayout& levelGrids = levelProjPtr->getBoxes();
      DataIterator ditLev = levelCorr.dataIterator();
      for (ditLev.reset(); ditLev.ok(); ++ditLev)
      {
        bcHolder.operator()(levelCorr[ditLev],
                            levelGrids[ditLev],
                            levelDomain,
                            levelDx,
                            false); // not homogeneous
      }
      levelCorr.exchange(corrComps);
      levelProjPtr = levelProjPtr->fineProjPtr();
    }

    // compute correction field
    computeGrad_eLambda(a_porosity);

    // Just need to add 1.0 now
    //    lambdaMult = 1.0/lambdaMult;
    // now return lambda to previous state
    for (int lev = m_level; lev<=finestLevel; lev++)
    {
      LevelData<FArrayBox>& levelLambda = *(a_lambda[lev]);
      DataIterator dit = levelLambda.dataIterator();

      for (dit.reset(); dit.ok(); ++dit)
      {
        //        levelLambda[dit] *= lambdaMult;
        //        levelLambda[dit] += 1.0;

        if (m_scale_lambda_with_porosity)
        {
          // Compute 1/porosity
          FArrayBox inversePorosity(levelLambda[dit].box(), SpaceDim);
          //          inversePorosity.copy((*a_porosity[lev])[dit]);
          EdgeToCell((*a_porosity[lev])[dit], inversePorosity);
          inversePorosity.invert(1);

          levelLambda[dit] += inversePorosity;
        }
        else
        {
          levelLambda[dit] += 1.0;
        }
      }
    }


    for (int lev = m_level; lev<vectorSize; lev++)
    {
      if (newLambda[lev] != NULL)
      {
        delete newLambda[lev];
        newLambda[lev] = NULL;
      }
    }

  } // end if eta > 0
}

// ---------------------------------------------------------------
void Projector::computeGrad_eLambda(Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosity           )
{
  // do this recursively starting at the finest level

  if (!isFinestLevel())
  {
    m_fineProjPtr->computeGrad_eLambda(a_porosity);
  }

  // now compute composite gradient on this level
  LevelData<FArrayBox>* crseBCPtr = NULL;

  if (m_crseProjPtr != NULL)
  {
    crseBCPtr = &(m_crseProjPtr->eLambda());
  }

  // recall that composite MAC gradient is the same as the level
  // gradient, since finer level is not considered to be part
  // of this level (take care of covered regions with avgDown)

  LevelData<FluxBox> temp_grad_eLambda(m_grad_eLambda.disjointBoxLayout(), m_grad_eLambda.nComp(), m_grad_eLambda.ghostVect());
  for (DataIterator dit = temp_grad_eLambda.dataIterator(); dit.ok(); ++dit)
  {
    temp_grad_eLambda[dit].setVal(0.0);
  }

  Gradient::levelGradientMAC(temp_grad_eLambda, m_eLambda,
                             crseBCPtr, m_dx,
                             m_gradIVS,
                             m_cfInterp);

  //  LevelData<FArrayBox> temp(m_grad_eLambda.disjointBoxLayout(), 1);
  //  Divergence::levelDivergenceMAC(temp, temp_grad_eLambda,m_dx);

  if (!s_applyVDCorrection)
  {
    setValLevel(temp_grad_eLambda, 0.0);
  }

  for (DataIterator dit = temp_grad_eLambda.dataIterator(); dit.ok(); ++dit)
  {
    FluxBox& porosity = (*a_porosity[m_level])[dit];
    if (m_scaleSyncCorrection || m_scale_lambda_err_with_porosity) // || scale_lambda_err_with_porosity ?
    {
      temp_grad_eLambda[dit].mult(porosity, porosity.box(), 0, 0, 1);
    }

    ParmParse pp("projection");
    bool divide_lambda_correction_porosity = false;
    pp.query("divide_lambda_correction_porosity", divide_lambda_correction_porosity);
    if (divide_lambda_correction_porosity)
    {
      temp_grad_eLambda[dit].divide(porosity, porosity.box(), 0, 0, 1);
    }


    m_grad_eLambda[dit].copy(temp_grad_eLambda[dit]);

  }


  // now average down finer level grad_eLambda if necessary
  if (!isFinestLevel())
  {
    const DisjointBoxLayout& fineGrids = m_fineProjPtr->getBoxes();
    int nRefFine = m_fineProjPtr->nRefCrse();
    int nComp = 1;

    CoarseAverageEdge avgDownObject(fineGrids, nComp, nRefFine);

    const LevelData<FluxBox> & fineEdgeGrad = m_fineProjPtr->grad_eLambda();

    avgDownObject.averageToCoarse(m_grad_eLambda, fineEdgeGrad);
  }
}

void Projector::applyFreestreamCorrection(LevelData<FluxBox>& a_advVel, Real scale)
{
  // Compute u = u + grad(e_lambda),
  // where grad(e_lambda) was calculated at the end of previous
  // timestep during synchronisation

  // scale is positive (1.0) for adding the correction,
  // and negative (-1.0) if we want to remove the correction

  for (DataIterator dit = a_advVel.dataIterator(); dit.ok(); ++dit)
  {

    for (int dir = 0; dir < SpaceDim; dir++)
    {
      a_advVel[dit][dir].plus(m_grad_eLambda[dit][dir], scale);
    }

  }
}

void Projector::removeFreestreamCorrection(LevelData<FluxBox>& a_advVel)
{
  applyFreestreamCorrection(a_advVel, -1);
}


// ---------------------------------------------------------------
void Projector::initialLevelProject(LevelData<FArrayBox>& a_velocity,
                                      LevelData<FArrayBox>* a_crseVelPtr,
                                      const Real a_newTime, const Real a_dt,
                                      const RefCountedPtr<LevelData<FArrayBox> > a_porosityPtr,
                                      const RefCountedPtr<LevelData<FArrayBox> > a_crsePorosityPtr,
                                      const RefCountedPtr<LevelData<FluxBox> > a_porosityEdgePtr,
                                      const RefCountedPtr<LevelData<FluxBox> > a_crsePorosityEdgePtr)
{
  if (s_verbosity >= 2)
  {
    pout() << "  CCProjector::initialLevelProject with dt = " << a_dt << endl;
  }
  LevelData<FArrayBox>* crseDataPtr=NULL;
  LevelData<FArrayBox> levelProjRHS(getBoxes(),1);

  // do coarse-fine BC
  if (m_level > 0)
  {
    if (doQuadInterp())
    {
      CH_assert(a_crseVelPtr != NULL);
      CH_assert(m_crseProjPtr != NULL);

      // will need to add dt*gradPi to crse vel
      const DisjointBoxLayout& crseLevelGrids = m_crseProjPtr->getBoxes();
      crseDataPtr = new LevelData<FArrayBox>(crseLevelGrids, SpaceDim);
      LevelData<FArrayBox>& crseData = *crseDataPtr;
      m_crseProjPtr->gradPi(crseData);

      // now multiply by dt and add crseVel
      DataIterator ditCrse = crseData.dataIterator();
      LevelData<FArrayBox>& crseVel = *a_crseVelPtr;

      // inherent assumption that crseVel & crseGradPi share same layout
      for (ditCrse.reset(); ditCrse.ok(); ++ditCrse)
      {
        FArrayBox& thisCrseVel = crseVel[ditCrse];
        FArrayBox& thisCrseGradPi = crseData[ditCrse];

        if (a_crsePorosityPtr!=NULL)
        {
          for (int dir=0; dir<SpaceDim; dir++)
          {
            thisCrseGradPi.mult((*a_crsePorosityPtr)[ditCrse],thisCrseGradPi.box(), 0,dir);
          }
        }

        // note that we multiply by this level's dt
        thisCrseGradPi *= a_dt;
        // now add crseVel
        thisCrseGradPi.plus(thisCrseVel);
      }

    }
    else
    {
      // alternative is not implemented at this point
      MayDay::Error("non-quadratic interpolation BC's not implemented");
    }
  } // end if coarser level exists

  // just to be safe...
  Interval velComps(0,SpaceDim-1);
  a_velocity.exchange(velComps);

  // now compute RHS for projection
  Divergence::levelDivergenceCC(levelProjRHS, a_velocity, crseDataPtr,
                                m_dx, doQuadInterp(), m_cfInterp);

  // for proper scaling of pi, divide this by dt
  DataIterator dit = levelProjRHS.dataIterator();
  Real dtScale = 1.0/a_dt;

  for (dit.reset(); dit.ok(); ++dit)
  {
    levelProjRHS[dit] *= dtScale;
  }

  // set up coarse BC's for solve, then solve
  const LevelData<FArrayBox>* crsePiPtr = NULL;
  if (m_crseProjPtr != NULL) crsePiPtr = &(m_crseProjPtr->Pi());

  // now solve for Pi
  solveMGlevel(m_Pi, crsePiPtr, levelProjRHS, a_porosityPtr, a_crsePorosityPtr,
               a_porosityEdgePtr, a_crsePorosityEdgePtr,
               true); // cell centred

  if (m_crseProjPtr != NULL)
  {
    // reapply coarse-fine BC's here if necessary
    m_cfInterp.coarseFineInterp(m_Pi, *crsePiPtr);
  }

  // apply appropriate physical BC's
  BCHolder bcHolder = m_physBCPtr->gradPiFuncBC();
  const DisjointBoxLayout& levelGrids = getBoxes();

  // loop over grids...
  for (dit.reset(); dit.ok(); ++dit)
  {
    bcHolder.operator()(m_Pi[dit],
                        levelGrids[dit],
                        m_domain,
                        m_dx,
                        false); // not homogeneous
  }

  Interval piComps(0,0);
  m_Pi.exchange(piComps);

  dtScale = a_dt;
  applyCCcorrection(a_velocity, &m_Pi, a_porosityPtr ,dtScale);

  // clean up storage
  if (crseDataPtr!=NULL)
  {
    delete crseDataPtr;
    crseDataPtr=NULL;
  }
}

// ---------------------------------------------------------------
void Projector::initialSyncProjection(Vector<LevelData<FArrayBox>* >& a_vel,
                                        Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity,
                                        const Real a_newTime, const Real a_dtSync,
                                        AMRMultiGrid<LevelData<FArrayBox> >& a_solver
)
{
  // actually, in the current algorithm, i _think_ this is the same
  // as the normal sync projection
  doSyncProjection(a_vel, a_porosity, a_newTime, a_dtSync, a_solver);
}

// ---------------------------------------------------------------
void Projector::doInitialSyncOperations(Vector<LevelData<FArrayBox>* >& a_vel,
                                          Vector<LevelData<FArrayBox>* >& a_lambda,
                                          Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosityFace,
                                          Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity,
                                          const Real a_newTime, const Real a_dtSync)
{


  // now can define multilevel solver
  AMRMultiGrid<LevelData<FArrayBox> >* ispSolverPtr = new
      AMRMultiGrid<LevelData<FArrayBox> >;

  defineMultiGrid(*ispSolverPtr, a_vel, a_porosityFace, false);

  // now call sync projection
  initialSyncProjection(a_vel, a_porosity, a_newTime, a_dtSync, *ispSolverPtr);

  // Need to delete ispSolverPtr->m_bottomSolver from AMRMultiGrid.
  delete ispSolverPtr;
  AMRMultiGrid<LevelData<FArrayBox> >* cvdcSolverPtr = new
      AMRMultiGrid<LevelData<FArrayBox> >;
  defineMultiGrid(*cvdcSolverPtr, a_lambda, a_porosityFace, true);

  // now do freestream preservation solve
  computeVDCorrection(a_lambda, a_porosityFace, a_newTime, a_dtSync, *cvdcSolverPtr);

  delete cvdcSolverPtr; //Kris R.

}

// ---------------------------------------------------------------
void Projector::initialVelocityProject(Vector<LevelData<FArrayBox>* >& a_vel,
                                         Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosityFace,
                                         Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity,
                                         bool a_homogeneousCFBC)
{
  CH_assert(m_level==0);

  // first need to define solver
  AMRMultiGrid<LevelData<FArrayBox> > bigSolver;
  defineMultiGrid(bigSolver,a_vel, a_porosityFace, false);

  initialVelocityProject(a_vel, a_porosity, bigSolver, a_homogeneousCFBC);


}

// ---------------------------------------------------------------
void Projector::initialVelocityProject(Vector<LevelData<FArrayBox>* >& a_vel,
                                         Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity,
                                         AMRMultiGrid<LevelData<FArrayBox> >& a_bigSolver,
                                         bool a_homogeneousCFBC)
{
  // set up temp storage
  int vectorSize = a_vel.size();

  Vector<LevelData<FArrayBox>* > projRHS(vectorSize,NULL);
  Vector<LevelData<FArrayBox>* > projCorr(vectorSize,NULL);

  Projector* levelProjPtr = this;

  int finestLevel = vectorSize - 1;

  // this is a bogus levelData pointer for the composite divergence
  LevelData<FArrayBox>* crseVelPtr = NULL;
  LevelData<FArrayBox>* fineVelPtr = NULL;
  Vector<int> nRefFineVect(vectorSize, -1);

  // loop over levels to allocate storage and compute RHS
  for (int lev=m_level; lev <= finestLevel; lev++)
  {
    const DisjointBoxLayout& levelGrids = a_vel[lev]->getBoxes();

    projRHS[lev] = new LevelData<FArrayBox>(levelGrids,1);
    projCorr[lev] = &(levelProjPtr->eSync());

    Real dx = levelProjPtr->dx();
    const ProblemDomain& levelDomain = levelProjPtr->dProblem();
    int nRefCrse = levelProjPtr->nRefCrse();
    int nRefFine = -1;

    if (lev > 0)
    {
      crseVelPtr = a_vel[lev-1];
    }
    if (lev < finestLevel)
    {
      fineVelPtr = a_vel[lev+1];
      nRefFine = levelProjPtr->fineProjPtr()->nRefCrse();
    }
    else
    {
      fineVelPtr = NULL;
    }
    nRefFineVect[lev] = nRefFine;

    // just in case...
    Interval velComps(0,SpaceDim-1);
    a_vel[lev]->exchange(velComps);

    // now compute divergence
    //    fineVelPtr = NULL; // this was just testing
    Divergence::compDivergenceCC(*projRHS[lev],*a_vel[lev],crseVelPtr,
                                 fineVelPtr,dx,nRefCrse,nRefFine,
                                 levelDomain,doQuadInterp());


    // increment levelProjPtr
    levelProjPtr= levelProjPtr->fineProjPtr();

  } // end loop over levels for computation of RHS

  // take care of C/F Boundary condition for lBase

  LevelData<FArrayBox>* crseBCPtr = NULL;
  if (m_level > 0)
  {
    CH_assert(m_crseProjPtr != NULL);
    if (a_homogeneousCFBC)
    {
      // need to define coarse BC and set to 0
      const DisjointBoxLayout& crseGrids = m_crseProjPtr->getBoxes();

      crseBCPtr = new LevelData<FArrayBox>(crseGrids,1);

      setValLevel(*crseBCPtr, 0.0);

      projCorr[m_level-1] = crseBCPtr;
    }
    else
    {
      // not sure if this should ever be called
      MayDay::Warning("inhomogeneous BC's on initialVelProject");
      projCorr[m_level-1] = &m_crseProjPtr->eSync();
    }

  } // end if level > 0

  Vector<DisjointBoxLayout> allGrids(vectorSize);
  for (int lev=m_level; lev<=finestLevel; lev++)
  {
    allGrids[lev] = projCorr[lev]->getBoxes();
  }

  Interval sumComps(0,0);
  Real sumRHS;
  sumRHS =  computeSum(projRHS, nRefFineVect, m_dx, sumComps, m_level);
  pout() << "Sum(RHS) for IVP solve = "
      << setiosflags(ios::scientific) << sumRHS << endl;

  // now solve!
  a_bigSolver.solve(projCorr, projRHS, finestLevel, m_level,
                    true); // initialize projCorr to zero

  if (s_verbosity >= 2)
  {
    pout() << "CCProjector::initialVelocityProject - finished solve " << endl;
  }

  // apply physical boundary conditions before
  // taking gradients
  Interval corrComps(0,0);
  levelProjPtr = this;
  BCHolder bcHolder = m_physBCPtr->gradESyncFuncBC();

  for (int lev=m_level; lev<=finestLevel; lev++)
  {

    // loop over grids
    LevelData<FArrayBox>& levelCorr = *projCorr[lev];
    const ProblemDomain& levelDomain = levelProjPtr->dProblem();
    Real levelDx = levelProjPtr->dx();
    const DisjointBoxLayout& levelGrids = allGrids[lev];
    // levelProjPtr->getBoxes();
    DataIterator ditLev = levelCorr.dataIterator();

    for (ditLev.reset(); ditLev.ok(); ++ditLev)
    {
      bcHolder.operator()(levelCorr[ditLev],
                          levelGrids[ditLev],
                          levelDomain,
                          levelDx,
                          false); // not homogeneous
    }

    levelCorr.exchange(corrComps);

    levelProjPtr = levelProjPtr->fineProjPtr();
  }

  // now apply correction
  Real scale = 1.0;

  applySyncCorrection(a_vel, a_porosity, scale, crseBCPtr);

  // delete temp storage -- RHS and crse BC
  for (int lev=m_level; lev<=finestLevel; lev++)
  {
    delete projRHS[lev];
    projRHS[lev] = NULL;
  }

  if (crseBCPtr != NULL)
  {
    delete crseBCPtr;
    crseBCPtr = NULL;
  }
}

void Projector::scalePwithPorosity(bool a_scalePwithPorosity)
{
  m_scalePressureWithPorosity = a_scalePwithPorosity;
}

// ---------------------------------------------------------------
void Projector::defineMultiGrid(AMRMultiGrid<LevelData<FArrayBox> >& a_solver,
                                  const Vector<LevelData<FArrayBox>* >& a_vel,
                                  Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosity,
                                  bool a_freestreamSolve)
{
  // assume that this level is lBase...
  int vectorSize = a_vel.size();

  // determine finest existing level
  Projector* levelProj = this;

  int finestLevel = vectorSize - 1;

  // put level grids togther into a vector
  Vector<DisjointBoxLayout> allGrids(vectorSize);
  Vector<ProblemDomain> allDomains(vectorSize);
  Vector<Real> allDx(vectorSize);
  Vector<int> allRefRatios(vectorSize);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(vectorSize);

  levelProj = this; // Now levelProj is at level m_level.

  //Since we are always going to need the grids and coefficients on this level, do that now
  allGrids[m_level] = a_vel[m_level]->getBoxes();
  aCoef[m_level] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(allGrids[m_level], 1));
  setValLevel(*aCoef[m_level], 1.0);

  // Now deal with finer levels
  for (int lev = m_level+1; lev<=finestLevel; lev++)
  {
    levelProj = levelProj->fineProjPtr(); //levelProj now points to the same level as index lev
    allGrids[lev] = a_vel[lev]->getBoxes();
    allRefRatios[lev-1] = levelProj->nRefCrse(); //ref ratio to next finer level belongs to that level for CCProjector

    aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(allGrids[lev], 1));
    setValLevel(*aCoef[lev], 1.0);
  }
  // Now levelProj is at level finestLevel.

  //Kris R.
  levelProj = this; //now levelProj is at m_level

  //now deal with the case where there are coarser levels than this one
  for(int lev = m_level-1; lev >=0; lev--) {

    allRefRatios[lev] = levelProj->nRefCrse(); //levelProj points to the level above lev

    levelProj = levelProj->crseProjPtr(); //lecvelProj now points to the same level as lev
    allGrids[lev] = levelProj->m_Pi.getBoxes();

    aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(allGrids[lev], 1));
    setValLevel(*aCoef[lev], 1.0);
  }

  //levelProj is now at level 0. This is crucial for making the LHS operator consistent
  CH_assert(levelProj->m_level == 0);

  ProblemDomain baseDomain = levelProj->dProblem(); //Kris R.
  Real baseDx = levelProj->dx(); //Kris R.

  // Old version -solves Lap(p) = div(u)
  AMRPoissonOpFactory localPoissonOpFactory;
  // physBCPtr = m_physBCPtr->FreestreamCorrBC();
  // physBCPtr = m_physBCPtr->SyncProjBC();
  localPoissonOpFactory.define(baseDomain,
                               allGrids,
                               allRefRatios,
                               baseDx,
                               ((a_freestreamSolve) ?
                                   m_physBCPtr->FreestreamCorrFuncBC() :
                                   m_physBCPtr->SyncProjFuncBC() ));

  // New version = solves div(porosity grad p) = div(u)
  Real alpha=0.0; Real beta=-1.0;
//  VCAMRPoissonOp2Factory opFact;
  AMRProjectionOpFactory opFact;
  opFact.define(baseDomain,
                allGrids,
                allRefRatios,
                baseDx,
                ((a_freestreamSolve) ?
                    m_physBCPtr->FreestreamCorrFuncBC() :
                    m_physBCPtr->SyncProjFuncBC() ),
                    alpha, aCoef, beta, a_porosity);


  opFact.m_relaxMode = s_multigrid_relaxation;


  makeBottomSolvers();

  if (m_scaleSyncCorrection)
  {
    a_solver.define(baseDomain,
                    opFact,
                    m_bottomSolverLevel,
                    finestLevel+1);
  }
  else
  {
    a_solver.define(baseDomain,
                    localPoissonOpFactory,
                    m_bottomSolverLevel,
                    finestLevel+1);
  }


  // We also want to use m_limitCoarsening:
  // if true, only do multigrid coarsening down to next coarser
  // AMR level (only coarsen by a_nRefCrse).
  // If false, coarsen as far as possible. Only relevant when lBase > 0.
  a_solver.m_verbosity = max(s_verbosity-1, 0);
  a_solver.m_eps = s_solver_tolerance;
  a_solver.m_pre = s_num_smooth_down; // smoothings before avging
  a_solver.m_post = s_num_smooth_up; // smoothings after avging
  a_solver.m_numMG = s_numMG;
  //a_solver.m_iterMax = 20;

}


// -------------------------------------------------------------
void Projector::postRestart()
{
  // set C/F BC's -- this depends on the fact that the coarser level
  // has already been read in!
  if (m_level > 0)
  {
    // need to do Pi
    CH_assert (m_cfInterp.isDefined());
    LevelData<FArrayBox>& crsePi = m_crseProjPtr->Pi();
    m_cfInterp.coarseFineInterp(m_Pi, crsePi);
  }

  Interval piComps = m_Pi.interval();
  m_Pi.exchange(piComps);

  // also need to do physical BC's
  BCHolder bcHolder = m_physBCPtr->gradPiFuncBC();
  const DisjointBoxLayout& levelGrids = getBoxes();

  // loop over grids
  DataIterator dit = m_Pi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    bcHolder.operator()(m_Pi[dit],
                        levelGrids[dit],
                        m_domain,
                        m_dx,
                        false); // not homogeneous
  }
}


void Projector::defineSolverMGLevelVCOp(
    const RefCountedPtr<LevelData<FArrayBox> > a_porosityPtr,
    const RefCountedPtr<LevelData<FArrayBox> > a_crsePorosityPtr,
    const RefCountedPtr<LevelData<FluxBox> > a_porosityEdgePtr,
    const RefCountedPtr<LevelData<FluxBox> > a_crsePorosityEdgePtr,
    bool cellCentred,
    Real beta)
{
  if (s_verbosity >= 5)
  {
    pout() << "CCProjector::defineSolverMGLevelVCOp "            << endl;
  }

  m_porosityPtr = a_porosityPtr;
  m_crsePorosityPtr = a_crsePorosityPtr;
  m_porosityEdgePtr = a_porosityEdgePtr;
  m_crsePorosityEdgePtr = a_crsePorosityEdgePtr;


  int numSolverLevels = 1;
  if (a_crsePorosityEdgePtr != NULL)
  {
    numSolverLevels= 2;
  }

  m_aCoef.resize(numSolverLevels);
  m_bCoef.resize(numSolverLevels);
  m_bCoefCC.resize(numSolverLevels);

  ProblemDomain baseDomain(m_domain); // on this level

  Real alpha=0;

  Real dxCrse=m_dx;

  if (a_crsePorosityEdgePtr != NULL)
  {
    m_aCoef[0] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_allGrids[0], 1, IntVect::Unit));
    m_aCoef[1] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_allGrids[1], 1, IntVect::Unit));
    //    m_bCoefCC[1] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_allGrids[1], 1, 2*IntVect::Unit));
    setValLevel(*m_aCoef[0], 1.0);
    setValLevel(*m_aCoef[1], 1.0);

    //    m_bCoef[0] = a_crsePorosityEdgePtr;
    //    m_bCoef[1] = a_porosityEdgePtr;
    m_bCoef[0] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_allGrids[0], 1, IntVect::Unit));
    m_bCoef[1] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_allGrids[1], 1, IntVect::Unit));
    m_bCoefCC[0] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_allGrids[0], 1, IntVect::Unit));
    m_bCoefCC[1] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_allGrids[1], 1, IntVect::Unit));

    //    a_porosityEdgePtr->copyTo(m_bCoef[0]);
    for (DataIterator dit=m_bCoef[0]->dataIterator(); dit.ok(); ++dit)
    {
      (*m_bCoef[0])[dit].copy((*a_crsePorosityEdgePtr)[dit]);
      //      (*m_bCoef[0])[dit].setVal(1.0);
      //      (*m_bCoefCC[0])[dit].setVal(1.0);

      //      (*m_bCoefCC[0])[dit].mult(-beta);

      for (int dir=0; dir<SpaceDim; dir++)
      {
        (*m_bCoef[0])[dit][dir].mult(-beta);
      }
    }

    for (DataIterator dit=m_bCoef[1]->dataIterator(); dit.ok(); ++dit)
    {
      (*m_bCoef[1])[dit].copy((*a_porosityEdgePtr)[dit]);
      //      (*m_bCoef[1])[dit].setVal(1.0);
      for (int dir=0; dir<SpaceDim; dir++)
      {
        (*m_bCoef[1])[dit][dir].mult(-beta);
      }
    }

    //    m_bCoefCC[0] = a_crsePorosityPtr;
    //    m_bCoefCC[1] = a_porosityPtr;

    baseDomain.coarsen(m_nRefCrse);
    dxCrse = m_nRefCrse * m_dx;
  }
  else
  {
    m_aCoef[0] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_allGrids[0], 1, IntVect::Unit));
    setValLevel(*m_aCoef[0], 1.0);
    //    m_bCoef[0] = a_porosityEdgePtr;
    m_bCoef[0] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_allGrids[0], 1, IntVect::Unit));
    for (DataIterator dit=m_bCoef[0]->dataIterator(); dit.ok(); ++dit)
    {
      (*m_bCoef[0])[dit].copy((*a_porosityEdgePtr)[dit]);
      //      (*m_bCoef[0])[dit].setVal(1.0);
      for (int dir=0; dir<SpaceDim; dir++)
      {
        (*m_bCoef[0])[dit][dir].mult(-beta);
      }
    }
    //    m_bCoefCC[0] = a_porosityPtr;
  }


  Vector<int> refRatios(1, m_nRefCrse);

  makeBottomSolvers();

  if (s_verbosity > 3)
  {
    pout() << "CCProjector::defineSolverMGlevel - define multigrids" << endl;
  }

//  VCAMRPoissonOp2Factory faceOpFact;
  int average_type = CoarseAverage::arithmetic;
  ParmParse pp("projection");
  pp.query("average_type", average_type);
  AMRProjectionOpFactory faceOpFact;

  faceOpFact.define(baseDomain,
                    m_allGrids,
                    refRatios,
                    dxCrse,
                    m_physBCPtr->LevelPressureFuncBC(),
                    alpha, m_aCoef, beta, m_bCoef,
                    average_type);
  faceOpFact.m_relaxMode = s_multigrid_relaxation;

  if (s_relax_bottom_solver)
  {
    m_solverMGlevel.define(baseDomain, // on either this level or coarser level
                           faceOpFact,
                           m_bottomSolverLevel,
                           numSolverLevels);
  }
  else
  {
    m_solverMGlevel.define(baseDomain, // on either this level or coarser level
                           faceOpFact,
                           m_BiCGBottomSolverLevel,
                           numSolverLevels);
  }

  setSolverParameters();


  //  m_numMG, m_normThresh;
}
void Projector::makeBottomSolvers()
{
  if (s_verbosity > 3)
  {
    pout() << "CCProjector::makeBottomSolvers()" << endl;
  }


  //Delete previous bottom solver -- Kris R.
  if (m_BiCGBottomSolverLevel != NULL)
  {
    delete m_BiCGBottomSolverLevel;
    m_BiCGBottomSolverLevel = NULL;
  }

  if (m_bottomSolverLevel!= NULL)
  {
    delete m_bottomSolverLevel;
    m_bottomSolverLevel = NULL;
  }

  BiCGStabSolver<LevelData<FArrayBox> >* newBottomPtr = new BiCGStabSolver<LevelData<FArrayBox> >;
  RelaxSolver<LevelData<FArrayBox> >* newRelaxBottomPtr = new RelaxSolver<LevelData<FArrayBox> >;
  newRelaxBottomPtr->m_verbosity = max(s_verbosity-2, 0);
  newBottomPtr->m_verbosity = max(s_verbosity-2, 0);

  newBottomPtr->m_imax = s_bottomSolveMaxIter;
  newRelaxBottomPtr->m_imax = s_bottomSolveMaxIter;

  m_BiCGBottomSolverLevel = newBottomPtr;
  m_bottomSolverLevel = newRelaxBottomPtr;
}
void Projector::setSolverParameters()
{
  m_solverMGlevel.m_verbosity = s_verbosity;
  m_solverMGlevel.m_eps = s_solver_tolerance;
  m_solverMGlevel.m_pre = s_num_smooth_down; // smoothings before avging
  m_solverMGlevel.m_post = s_num_smooth_up; // smoothings after avging
  m_solverMGlevel.m_bottom = s_num_precond_smooth; // smoothing before bottom solve
  m_solverMGlevel.m_numMG = s_numMG;
  m_solverMGlevel.m_hang = s_solver_hang;

//  if (s_solver_tolerance > 0)
//  {
//    m_solverMGlevel.m_eps = s_solver_tolerance;
//  }
}

// -------------------------------------------------------------
void Projector::defineSolverMGlevel(const DisjointBoxLayout& a_grids,
                                      const DisjointBoxLayout* a_crseGridsPtr)
{
  if (s_verbosity > 3)
  {
    pout() << "CCProjector::defineSolverMGlevel" << endl;
  }

  int numSolverLevels = 1;
  if (a_crseGridsPtr != NULL)
  {
    numSolverLevels= 2;
  }

  AMRPoissonOpFactory localPoissonOpFactory;

  m_allGrids = Vector<DisjointBoxLayout>(numSolverLevels);

  if (a_crseGridsPtr != NULL)
  {
    m_allGrids[0] = *a_crseGridsPtr;
    m_allGrids[1] = a_grids;
  }
  else
  {
    m_allGrids[0] = a_grids;
  }


  Vector<int> refRatios(1, m_nRefCrse);

  ProblemDomain baseDomain(m_domain); // on this level

  // Note that AMRPoissonOp and VCAMR.. require a different sign for beta
//  Real alpha=0.0;
//  Real beta = 1.0;

  if (a_crseGridsPtr != NULL)
  { // coarser level exists:  define solver on two levels


    baseDomain.coarsen(m_nRefCrse);


    m_allGrids[0] = *a_crseGridsPtr;
    m_allGrids[1] = a_grids;

    // this returns zero for me:
    // Real dxCrse = m_crseProjPtr->dx();
    Real dxCrse = m_nRefCrse * m_dx;



    if (s_verbosity > 3)
    {
      pout() << "CCProjector::defineSolverMGlevel - define Poisson Op Factory" << endl;
    }
    localPoissonOpFactory.define(baseDomain,
                                 m_allGrids,
                                 refRatios,
                                 dxCrse,
                                 m_physBCPtr->LevelPressureFuncBC());
//                                 alpha, beta);

  }
  else
  { // no coarser level:  define solver on only one level

    localPoissonOpFactory.define(m_domain,
                                 a_grids,
                                 m_dx,
                                 m_physBCPtr->LevelPressureFuncBC());
//                                 alpha, beta);

  }

  // You really need to delete this when you're done with a_solver.

  makeBottomSolvers();

  if (s_verbosity > 3)
  {
    pout() << "CCProjector::defineSolverMGlevel - define multigrids" << endl;
  }

  if (s_relax_bottom_solver)
  {
    m_solverMGlevel.define(baseDomain, // on either this level or coarser level
                           localPoissonOpFactory,
                           m_bottomSolverLevel,
                           numSolverLevels);
  }
  else
  {
    m_solverMGlevel.define(baseDomain, // on either this level or coarser level
                           localPoissonOpFactory,
                           m_BiCGBottomSolverLevel,
                           numSolverLevels);
  }


  setSolverParameters();


}

void Projector::solveMGlevel(LevelData<FArrayBox>&   a_phi,
                               const LevelData<FArrayBox>*   a_phiCoarsePtr,
                               const LevelData<FArrayBox>&   a_rhs,
                               bool cellCentred)
{
RefCountedPtr<LevelData<FArrayBox> > a_pressureScalePtr,  a_crsePressureScalePtr;
RefCountedPtr<LevelData<FluxBox> > a_pressureScaleEdgePtr,  a_crsePressureScaleEdgePtr;

solveMGlevel(a_phi, a_phiCoarsePtr, a_rhs, a_pressureScalePtr,  a_crsePressureScalePtr, a_pressureScaleEdgePtr,  a_crsePressureScaleEdgePtr,
             cellCentred);
}

// -------------------------------------------------------------
int Projector::solveMGlevel(LevelData<FArrayBox>&   a_phi,
                               const LevelData<FArrayBox>*   a_phiCoarsePtr,
                               const LevelData<FArrayBox>&   a_rhs,
                               const RefCountedPtr<LevelData<FArrayBox> > a_pressureScalePtr,
                               const RefCountedPtr<LevelData<FArrayBox> > a_crsePressureScalePtr,
                               const RefCountedPtr<LevelData<FluxBox> > a_pressureScaleEdgePtr,
                               const RefCountedPtr<LevelData<FluxBox> > a_crsePressureScaleEdgePtr,
                               bool cellCentred)
{
  if (s_verbosity >= 5)
  {
    pout() << "CCProjector::solveMGlevel "            << endl;
  }

  // Turn off initialising phi by default - it's dodgy with regridding
  ParmParse pp("projection");
  bool initPhi = false;
  pp.query("initPhi", initPhi);

  // Initialize a_phi to zero.
  if (!initPhi)
  {
    setValLevel(a_phi, 0.0);
  }

  Vector< LevelData<FArrayBox>* > phiVect;
  Vector< LevelData<FArrayBox>* > rhsVect;
  LevelData<FArrayBox>& rhsRef =
      const_cast< LevelData<FArrayBox>& >(a_rhs);

  int maxLevel;
  if (a_phiCoarsePtr != NULL)
  {
    // try initializing a_phi to the crseBC interpolated to here
    if (initPhi)
    {
      FineInterp interp(a_phi.disjointBoxLayout(), 1, m_nRefCrse, m_domain);
      interp.interpToFine(a_phi, *a_phiCoarsePtr);
    }

    maxLevel = 1;
    LevelData<FArrayBox>* phiCoarsePtrRef =
        const_cast< LevelData<FArrayBox>* >(a_phiCoarsePtr);
    phiVect.push_back(phiCoarsePtrRef);
    rhsVect.push_back(NULL); // I don't think this will be used


  }
  else
  {
    maxLevel = 0;
  }
  phiVect.push_back(&a_phi);
  rhsVect.push_back(&rhsRef);

  Real beta=1;

  if (a_pressureScaleEdgePtr!=NULL)
  {
    defineSolverMGLevelVCOp(a_pressureScalePtr, a_crsePressureScalePtr,
                            a_pressureScaleEdgePtr,  a_crsePressureScaleEdgePtr,
                            cellCentred, beta);
  }

  if (s_verbosity >= 5)
  {
    pout() << "CCProjector::solveMGlevel - do actual solve "            << endl;
  }


  // l_max = maxLevel, l_base = maxLevel
  m_solverMGlevel.solve(phiVect, rhsVect, maxLevel, maxLevel,
                        false); // don't initialize to zero



//  if (s_verbosity >= 2
//      || m_solverMGlevel.m_exitStatus != 0)
//  {
//    if (cellCentred)
//    {
//      pout() << "  CC Projection - solveMGlevel - solved, exitStatus =  " <<  m_solverMGlevel.m_exitStatus           << endl;
//    }
//    else
//    {
//      pout() << "  MAC Projection - solveMGlevel - solved, exitStatus =  " <<  m_solverMGlevel.m_exitStatus           << endl;
//    }
//  }

  //	int exitStatus = m_solverMGlevel.m_exitStatus;

  return m_solverMGlevel.m_exitStatus;
}
