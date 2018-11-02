#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//#ifndef _PhysBCUtil_H_
//
//#endif

#ifndef _CCProjectorComp_H_
#define _CCProjectorComp_H_
//#define PHYSBCUTIL_H



#include <cstdlib>
#include <cmath>
#include <cassert>

#include <iostream>
using std::istream;
using std::cout;
using std::cerr;

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "AMRMultiGrid.H"
#include "QuadCFInterp.H"
#include "PhysBCUtil.H"
#include "RelaxSolver.H" 
#include "FluxBox.H"
#include "ParmParse.H"
#include "AMRPoissonOp.H"
#include "VCAMRPoissonOp2.H"
#include "BiCGStabSolver.H"
#include "PiecewiseLinearFillPatchFace.H"

#include "UsingNamespace.H"

/// this class manages the various forms of the CC projection
class CCProjectorComp
/** this class performs the various MAC composite projections required
by the MushyLayer code, all on a composite AMR hierarchy */
{
public:
  ///
  CCProjectorComp();

  /// full constructor
  CCProjectorComp(const Vector<DisjointBoxLayout>& a_amrGrids,
              const Vector<ProblemDomain>& a_amrDomains,
              const Vector<Real>&  a_amrDx,
              const Vector<int>& a_refinementRatios,
			  const int a_finest_level,
              PhysBCUtil& a_physBC,
			  Vector<RefCountedPtr<LevelData<FArrayBox> > > a_theta,
			  Vector<RefCountedPtr<LevelData<FArrayBox> > > a_permeability,
			  int a_numGhost);

  /// destructor
  ~CCProjectorComp();

  void undefine();

  /// define function
  void define(const Vector<DisjointBoxLayout>& a_amrGrids,
          const Vector<ProblemDomain>& a_amrDomains,
          const Vector<Real>&  a_amrDx,
          const Vector<int>& a_refinementRatios,
		  const int a_finest_level,
          PhysBCUtil& a_physBC,
		  Vector<RefCountedPtr<LevelData<FArrayBox> > > a_theta,
		  Vector<RefCountedPtr<LevelData<FArrayBox> > > a_permeability,
		  int a_numGhost);

  /// initialize new projection with data from old projection
  void init(const CCProjectorComp& a_oldProj);

  /// define static parts
  void variableSetUp();

  ///
  void setCrseProj(CCProjectorComp* a_crseProj, int nRefCrse);

  ///
  void setFineProj(CCProjectorComp* a_fineProj);

  /// set solver parameter
  void limitSolverCoarsening(bool a_limitSolverCoarsening);

  /// write checkpoint header
  void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// write this class to a checkpoint file for later restart
  void  writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// read the checkpoint header
  void readCheckpointHeader(HDF5Handle& a_handle);

  /// read this class from a checkpoint file
  void readCheckpointLevel(HDF5Handle& a_handle);

  /// do level MAC projection, correct uEdge.
  /*** a_oldTime is
       old-time on AMRLevel which is calling the projection.
       assumes BC's already set on uEdge values */
  void levelMacProject(LevelData<FluxBox>& uEdge,
                       Real a_oldTime, Real a_dt);

  /// do level projection and correct (cell-centered) velocities;
  /*** assumes physical and exchange (copy) boundary
       conditions already set.
  */
  void LevelProject(LevelData<FArrayBox>& a_velocity,
                    LevelData<FArrayBox>* a_crseVelPtr,
                    const Real a_newTime, const Real a_dt);

  /** defines AMRSolver and calls sync projection
      and freestream preservation solve. assumes all physical
      and copy BC's already set (does C/F BC's however) */
  void doSyncOperations(Vector<LevelData<FArrayBox>* >& a_velocity,
//                        Vector<LevelData<FArrayBox>* >& a_lambda,
                        const Real a_newTime, const Real a_dtSync);

  /// Initialization functions
  /** performs initial level projection -- phys and copy BC's
      should be preset */
  void initialLevelProject(LevelData<FArrayBox>& a_velocity,
                           LevelData<FArrayBox>* a_crseVelPtr,
                           const Real a_oldTime, const Real a_newTime);

  /** defines AMRSolver and calls initial sync projection and
      freestream preservation solves */
  void doInitialSyncOperations(Vector<LevelData<FArrayBox>* >& a_vel,
                               Vector<LevelData<FArrayBox>* >& a_lambda,
                               const Real a_newTime, const Real a_dtSync);

  void  doLambdaCorrection(Vector<LevelData<FluxBox>* >& a_vel,
			Vector<LevelData<FArrayBox> * >& a_lambda,
			Real a_newTime, Real a_dtSync);

  /// performs multilevel projection on velocity, modifies velocity in place
  /**  dtSync is used as a scaling factor for coarse-level BC's,
       not used in actual projection.  if homogeneousCFBC is false, C/F bc
       for l_base is taken from coarse-level eSync, o/w use l_base-1
       correction field = 0.  if applyCorrection == false, then correction
       is computed, but not applied (can be useful as a diagnostic).
       physical and copy BC's should already be set before this is called.
  */

  void initialVelocityProject(Vector<LevelData<FArrayBox>* >& a_velocity,
                              bool a_homogeneousCFBC = true);

  /// same as other initialVelocityProject, but uses a pre-defined AMRSolver
  void initialVelocityProject(Vector<LevelData<FArrayBox>* >& a_velocity,
                              AMRMultiGrid<LevelData<FArrayBox> >& a_solver,
                              bool a_homogeneousCFBC = true);

  ////  performs the necessary re-initialization after regridding
  /*** does composite projection on new velocity field
       and also recomputes e_lambda and grad_eLambda */
  void doPostRegridOps(Vector<LevelData<FArrayBox>* >& a_velocity,
                       Vector<LevelData<FArrayBox>* >& a_lambda,
                       const Real a_dt, const Real a_time);

  /// access functions

  ///
  int getLevel() const;

  ///
  const ProblemDomain& dProblem() const;

  ///
  const DisjointBoxLayout& getBoxes() const;

  ///
  int nRefCrse() const;

  ///
  Real dx() const;

  /// returns coefficient for volume-discrepancy solve.
  Real etaLambda() const;

  ///
  CCProjectorComp* fineProjPtr() const;

  ///
  CCProjectorComp* crseProjPtr() const;

  ///
  void getPressure(Vector<RefCountedPtr<LevelData<FArrayBox> > >& pi);

  void getGradPressure(Vector<RefCountedPtr<LevelData<FArrayBox> > >& gradP);
  void getGradPressure(Vector<LevelData<FArrayBox> * >& gradP);
  void getGradPressureEdge(Vector<LevelData<FluxBox> * >& gradPedge);

  /// returns edge-centered grad(phi) in direction dir
  void gradPhi(LevelData<FArrayBox>& a_gradPhi, int a_dir) const;

  /// returns all components of grad(phi) (gradPhi should be correct size)
  void gradPhi(LevelData<FluxBox>& a_gradPhi) const;

  void gradPhiComp(Vector<RefCountedPtr<LevelData<FluxBox> > >& a_gradPhi) const;

  /// returns level-projection pressure
  Vector<RefCountedPtr<LevelData<FArrayBox> > >& Pi();

  /// returns grad(pi) in direction dir
  void gradPi(LevelData<FArrayBox>& a_gradPi, int a_dir) const;

  /// returns grad(pi) in all directions into SpaceDim-dimensioned gradPi
  void gradPi(LevelData<FArrayBox>& a_gradPi) const;

  /// returns synchronization correction
  LevelData<FArrayBox>& eSync();

  /// returns synchronization correction (const version)
  const LevelData<FArrayBox>& eSync() const;

  /// returns cell-centered grad(eSync) (composite gradient)
  void grad_eSync(LevelData<FArrayBox>& a_grad_eSync, int a_dir) const;

  /// returns cell-centered G^{comp}(eSync)
  void grad_eSync(LevelData<FArrayBox>& a_grad_eSync) const;

  /// returns volume-discrepancy correction
  LevelData<FArrayBox>& eLambda();

  /// returns volume-discrepancy correction (const version)
  const LevelData<FArrayBox>& eLambda() const;

  /// returns edge-centered grad(eLambda) in direction dir
  void grad_eLambda(LevelData<FArrayBox>& a_grad_eLambda, int a_dir) const;

  /// returns edge-centered grad(eLambda) in all directions
  /** non-constant reference returned because of need to rescale. */
  LevelData<FluxBox>&  grad_eLambda();

  ///
  const LevelData<FluxBox>& grad_eLambda() const;

  /// do sync projection?
  bool doSyncProjection() const;

  /// do volume discrepancy correction?
  bool doMacSync() const;

  /// has this object been completely initialized?
  bool isInitialized() const;

  /// use quadratic interpolation (instead of extrap) for c/f bc
  bool doQuadInterp() const;

  /// returns predefined quadratic coarse-fine interpolation object
  QuadCFInterp& quadCFInterpolator();

  /// returns predefined intvectset which shows coverage
  const LayoutData<IntVectSet>& getGridIVS();

  /// is this the finest level?
  bool isFinestLevel() const;

  /// set whether this is the finest level
  void isFinestLevel(bool a_finest_level);

  /// set verbosity
  void verbosity(int a_verbosity);

  /// returns verbosity
  int verbosity() const;

  /// sets physBCs
  void setPhysBC(PhysBCUtil& a_bc);

  ///
  PhysBCUtil* getPhysBCPtr() const;

  void projectVelocity(Vector<LevelData<FArrayBox> *> a_U,
		  Vector<RefCountedPtr<LevelData<FArrayBox> > > a_permeability,
		  Vector<LevelData<FluxBox>* > a_Uedge,
		  const Vector<LevelData<FArrayBox> *> a_Ustar,
		  const int order);

  void getDivUStar(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_div);

  void averageCoarseToFineFace(Vector<LevelData<FluxBox>* > a_phi);

protected:


  /// returns uEdge - G(phi_mac); C/F BC is CFscale*Pi(crse)
  void applyMacCorrection(LevelData<FluxBox>& a_uEdge, Real CFscale, int lev);

  /// vel := vel - scale*G_CC(pi)
  void correctCCVelocities(LevelData<FArrayBox>& a_velocity,
                           const Real scale) const;

  /// do sync projection with this level as lBase; corrects velocity in place
  void doSyncProjection(Vector<LevelData<FArrayBox>* >& a_velocity,
                        const Real a_newTime, const Real a_dtSync,
                        AMRMultiGrid<LevelData<FArrayBox> >& a_solver
                        );

  /// perform initial sync projection
  void initialSyncProjection(Vector<LevelData<FArrayBox>* >& a_velocity,
                             const Real a_newTime, const Real a_dtSync,
                             AMRMultiGrid<LevelData<FArrayBox> >& a_solver
                             );

  /// do volume discrepancy solve, with this level as lBase
  void computeVDCorrection(Vector<LevelData<FArrayBox> * >& a_lambda,
                           const Real a_newTime,
                           const Real a_dtSync,
                           AMRMultiGrid<LevelData<FArrayBox> >& a_solver,
						   Vector<LevelData<FluxBox>* >& a_vel
                           );

  ///  apply sync correction vel = vel - scale*G^{comp} eSync
  /** (recursively applies correction to all finer levels) */
  void applySyncCorrection(Vector<LevelData<FArrayBox>* >& a_velocity,
                           const Real scale,
                           LevelData<FArrayBox>* crseCorr);

  /// compute composite gradient of eLambda -- does this recursively
  void computeGrad_eLambda();

  /// rescales composite gradient of eLambda, recursively all finer levels
  /** rescales grad(eLambda) by (a_dx_lbase)/m_dx for this level and
      recursively calls this function for finer levels.  This is done
      for stability reasons.
  */
  void rescaleGrad_eLambda(Real a_dx_lbase);

  /// defines Multilevel AMRSolver used for composite projections
  /** uses grids in vel, lBase is this level.
      a_freestream solver is true if solver will be used for
      freestream preservation solve. */
  void defineMultiGrid(AMRMultiGrid<LevelData<FArrayBox> >& a_solver,
                       LinearSolver<LevelData<FArrayBox> >* a_bottomSolver, 
                       bool a_freestreamSolver);

  /// defines m_solverMGlevel
  void defineSolverMGlevel(const DisjointBoxLayout& a_grids,
		  	  	  	  	   const Vector<RefCountedPtr<LevelData<FluxBox> > > a_permeability,
						   const Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef,
                           const DisjointBoxLayout* a_crseGridsPtr,
						   const int lev);
  void defineSolverMG();

  /// solves using m_solverMGlevel
  void solveMGlevel(LevelData<FArrayBox>&         a_phi,
                    const LevelData<FArrayBox>*   a_phiCoarsePtr,
                    const LevelData<FArrayBox>&   a_rhs);

  /// takes care of setting ghost cell values after restarting
  void postRestart();



private:

  ///
  int m_level;

  ///
  Vector<Real> m_amrDx;

  ///
  Vector<DisjointBoxLayout> m_amrGrids;

  ///
  Vector<ProblemDomain> m_amrDomains;

  /// Quadratic C/F interpolation object
  QuadCFInterp m_cfInterp;
  Vector<QuadCFInterp*> m_amrCFInterps;

  // Linear C/F interpolation object
  Vector<PiecewiseLinearFillPatchFace*> m_amrFillPatchFace;

  /// cache intvectsets for gradient here (if relevant)
  LayoutData<IntVectSet> m_gradIVS;
  Vector<LayoutData<IntVectSet>* > m_amrGradIVS;

  /// if coarse level exists, refinement ratio down to it
  Vector<int> m_refinement_ratios;

  /// cell centered pressure data at time n+1/2
//  LevelData<FArrayBox> m_Pi;
  Vector<LevelData<FArrayBox> * > m_pressure;
  Vector<LevelData<FArrayBox> * > m_gradPressure;
  Vector<LevelData<FluxBox>* > m_gradPressureEdge;

  ///Permeability, a collection of pointers to the permeability passed into the operator
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_permeability;

  /// To store div(U^*)
  Vector<LevelData<FArrayBox>* > m_divUstar;

  /// MAC correction field
//  LevelData<FArrayBox> m_phi;
  Vector<LevelData<FArrayBox> * > m_phi;

  /// synchronization correction
  Vector<LevelData<FArrayBox>* > m_eSync;
  /// volume-discrepancy correction
  Vector<LevelData<FArrayBox>* > m_eLambda;

  /// edge-centered gradient of eLambda
  Vector<LevelData<FluxBox>* > m_grad_eLambda;

  /// For BCs
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_theta;

  ///
  AMRMultiGrid<LevelData<FArrayBox> > m_solverMGlevel;
  AMRMultiGrid<LevelData<FArrayBox> > m_solverMG;

  BiCGStabSolver<LevelData<FArrayBox> >* m_bottomSolver;

//  AMRPoissonOpFactory m_solverOpFact;
  VCAMRPoissonOp2Factory m_solverOpFact;

  ///
  static bool m_doSyncProjection;

  ///
  static bool m_applySyncCorrection;

  ///
  bool m_isInitialized;

  ///
  static bool m_doQuadInterp;

  /// solver parameter -- passed on through when defining solvers
  bool m_limitSolverCoarsening;

  ///
  static Real m_etaLambda;

  /// multigrid solver tolerance
  static Real s_solver_tolerance;

  /// multigrid solver parameter
  static int s_num_smooth_up;

  /// multigrid solver parameter
  static int s_num_smooth_down;

  /// number of smoothings in preconditioner for projection
  static int s_num_precond_smooth;

  /// use consistent scaling for lamdba solve?
  static bool s_constantLambdaScaling;

  /// timestep used for all subcycled lamdba solves
  static Real s_lambda_timestep;

  /// have the static variables been initialized (ParmParse)
  static bool pp_init;

  /// verbosity of output
  static int s_verbosity;

  /// if false, VD correction is set to zero once it's computed
  static bool s_applyVDCorrection;

  /// is this the finest extant level?
  bool m_finest_level;
  int finest_level;
  int num_ghost;

  ///bottom solver for level solves
  RelaxSolver<LevelData<FArrayBox> > * m_relaxBottomSolver;
  BiCGStabSolver<LevelData<FArrayBox> > * m_bottomSolverLevel;
//  RelaxSolver<LevelData<FArrayBox> > * m_bottomSolverLevel;

  ///
  PhysBCUtil* m_physBCPtr;
};

#endif
