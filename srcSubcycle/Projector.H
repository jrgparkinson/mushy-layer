#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Projector_H_
#define _Projector_H_

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
#include"BiCGStabSolver.H"

#include "UsingNamespace.H"

/// this class manages the various forms of the projection
class Projector
/** this class performs the various cell-centered projections required
by the IAMR algorithm. These include the single-level operators levelMac
and levelProject, as well as the multilevel sync projection and volume
discrepancy solves. */
{
public:
  /// default constructor
  Projector();

  /// full constructor
  Projector(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_crseGridsPtr,
              const Box& a_domain,
              const Real  a_dx,
              Projector* a_finerProj,
              Projector* a_crseProj,
              int a_nRefCrse,
              int a_level,
              PhysBCUtil& a_physBC);

  /// full constructor
  Projector(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_crseGridsPtr,
              const ProblemDomain& a_domain,
              const Real  a_dx,
              Projector* a_finerProj,
              Projector* a_crseProj,
              int a_nRefCrse,
              int a_level,
              PhysBCUtil& a_physBC);

  /// Copy pressure from another projection object
  void copyPressure(Projector& a_proj);

  /// May need to rescale pressures after changing the timestep?
  void rescalePressure(Real a_oldDt, Real a_newDt);

  /// destructor
  ~Projector();

  /// define function
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_crseGridsPtr,
              const Box& a_domain,
              const Real  a_dx,
              Projector* a_finerProj,
              Projector* a_crseProj,
              int a_nRefCrse,
              int a_level,
              PhysBCUtil& a_physBC,
              bool a_usePiAdvectionBCs = true);

  /// define function
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_crseGridsPtr,
              const ProblemDomain& a_domain,
              const Real  a_dx,
              Projector* a_finerProj,
              Projector* a_crseProj,
              int a_nRefCrse,
              int a_level,
              PhysBCUtil& a_physBC,
              bool a_usePiAdvectionBCs = true);

  /// initialize new projection with data from old projection
  void init(const Projector& a_oldProj);

  /// define static parts
  void variableSetUp();

  /// set coarse projection operator
  void setCrseProj(Projector* a_crseProj, int nRefCrse);

  /// set fine projection operator
  void setFineProj(Projector* a_fineProj);

  /// set whether to scale pressure with porosity
  void scalePwithPorosity(bool a_scalePwithPorosity);

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
  int levelMacProject(LevelData<FluxBox>& uEdge,
                       Real a_oldTime, Real a_dt,
						const RefCountedPtr<LevelData<FArrayBox> > a_porosityPtr,
						const RefCountedPtr<LevelData<FArrayBox> > a_crsePorosityPtr,
						const RefCountedPtr<LevelData<FluxBox> > a_porosityEdgePtr,
						const RefCountedPtr<LevelData<FluxBox> > a_crsePorosityEdgePtr);

  /// Scale for \f$ \phi \f$
  /**
   * Important to get this right for coarse-fine interface matching
   */
  Real getPhiScale(Real a_dt);

  /// Check if the divergence of the velocity is sufficiently small
  /**
   * If not, increase the number of smoothing steps
   */
  void checkDivergence(LevelData<FluxBox>& a_uEdge);

  /// Return the scale for \f$ \phi \f$
  Real getScale(Real a_scale, Real a_dt);

  /// Set up for nonsubcycled solves
  void setNonSubcycledMACBCs();

  /// Set up for subcycled solves
  void setSubcycledMACBCs();

  /// Projection function that let's user specify pressure coarse fine boundary condition
  int levelMacProject(LevelData<FluxBox>& a_uEdge,
                       Real a_dt,
                       const RefCountedPtr<LevelData<FArrayBox> > a_crsePressurePtr,
                       const RefCountedPtr<LevelData<FArrayBox> > a_porosityPtr,
                       const RefCountedPtr<LevelData<FArrayBox> > a_crsePorosityPtr,
                       const RefCountedPtr<LevelData<FluxBox> > a_porosityEdgePtr,
                       const RefCountedPtr<LevelData<FluxBox> > a_crsePorosityEdgePtr,
                       bool alreadyhasPhi = false);


  /// Refactored this so we can use for an implicit MAC solve
  void getCFBC(LevelData<FArrayBox>& crseVelPtr, LevelData<FArrayBox>* a_crseVelPtr,
               LevelData<FArrayBox>* a_crsePorosityPtr, Real a_dt);

  /// Refactored to use in implicit MAC solve
  void scaleRHS(LevelData<FArrayBox>& a_rhs, Real a_scale);

  /// do level projection and correct (cell-centered) velocities;
  /*** assumes physical and exchange (copy) boundary
       conditions already set.
  */
  void LevelProject(LevelData<FArrayBox>& a_velocity,
                    LevelData<FArrayBox>* a_crseVelPtr,
                    const Real a_newTime, const Real a_dt,
					const RefCountedPtr<LevelData<FArrayBox> > a_porosityPtr,
					const RefCountedPtr<LevelData<FArrayBox> > a_crsePorosityPtr,
					const RefCountedPtr<LevelData<FluxBox> > a_porosityEdgePtr,
					const RefCountedPtr<LevelData<FluxBox> > a_crsePorosityEdgePtr,
					const bool a_isViscous);

  /// For doing a further (constant coefficient) projection
  void AdditionalLevelProject(LevelData<FArrayBox>& a_velocity,
                              LevelData<FArrayBox>* a_crseVelPtr,
                              const Real a_newTime, const Real a_dt,
                              const bool a_isViscous);

  /// do synchronisation operations
  /** defines AMRSolver and calls sync projection
      and freestream preservation solve. assumes all physical
      and copy BC's already set (does C/F BC's however) */
  void doSyncOperations(Vector<LevelData<FArrayBox>* >& a_velocity,
                        Vector<LevelData<FArrayBox>* >& a_lambda,
                        Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosityFace,
                        Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity,
                        const Real a_newTime, const Real a_dtSync);

  /// Initialization functions
  /** performs initial level projection -- phys and copy BC's
      should be preset */
  void initialLevelProject(LevelData<FArrayBox>& a_velocity,
                           LevelData<FArrayBox>* a_crseVelPtr,
                           const Real a_oldTime, const Real a_newTime,
                           const RefCountedPtr<LevelData<FArrayBox> > a_porosityPtr,
                           const RefCountedPtr<LevelData<FArrayBox> > a_crsePorosityPtr,
                           const RefCountedPtr<LevelData<FluxBox> > a_porosityEdgePtr,
                           const RefCountedPtr<LevelData<FluxBox> > a_crsePorosityEdgePtr);

  /** defines AMRSolver and calls initial sync projection and
      freestream preservation solves */
  void doInitialSyncOperations(Vector<LevelData<FArrayBox>* >& a_vel,
                               Vector<LevelData<FArrayBox>* >& a_lambda,
                               Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosityFace,
                               Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity,
                               const Real a_newTime, const Real a_dtSync);

  /// performs multilevel projection on velocity, modifies velocity in place
  /**  dtSync is used as a scaling factor for coarse-level BC's,
       not used in actual projection.  if homogeneousCFBC is false, C/F bc
       for l_base is taken from coarse-level eSync, o/w use l_base-1
       correction field = 0.  if applyCorrection == false, then correction
       is computed, but not applied (can be useful as a diagnostic).
       physical and copy BC's should already be set before this is called.
  */

  void initialVelocityProject(Vector<LevelData<FArrayBox>* >& a_velocity,
                              Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosityFace,
                              Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity,
                              bool a_homogeneousCFBC = true);

  /// same as other initialVelocityProject, but uses a pre-defined AMRSolver
  void initialVelocityProject(Vector<LevelData<FArrayBox>* >& a_velocity,
                              Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity,
                              AMRMultiGrid<LevelData<FArrayBox> >& a_solver,
                              bool a_homogeneousCFBC = true);

  ///  performs the necessary re-initialization after regridding
  /*** does composite projection on new velocity field
       and also recomputes e_lambda and grad_eLambda */
  void doPostRegridOps(Vector<LevelData<FArrayBox>* >& a_lambda,
                       Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosity,
                       const Real a_dt, const Real a_time, const Real a_etaScale=1.0);

  /// access functions

  /// get the level this Projector is defined for
  int getLevel() const;

  /// get the problem domain
  const ProblemDomain& dProblem() const;

  /// get the DisjointBoxLayout for this projection operator
  const DisjointBoxLayout& getBoxes() const;

  /// get the refinement ratio to the coarser level
  int nRefCrse() const;

  /// get the cell spacing
  Real dx() const;

  /// returns coefficient for volume-discrepancy solve.
  Real etaLambda() const;

  /// We may want to change the value of \f$ \eta \f$ during a simulation
  void etaLambda(Real a_eta);

  /// pointer to fine projection operator
  Projector* fineProjPtr() const;

  /// pointer to coarse projection operator
  Projector* crseProjPtr() const;

  /// returns MAC correction
  LevelData<FArrayBox>& phi();

  /// Get \f$ \phi \f$
  void getPhi(LevelData<FArrayBox>& a_phi);

  /// const version of accessor
  const LevelData<FArrayBox>& phi() const;

  /// returns edge-centered grad(phi) in direction dir
  void gradPhi(LevelData<FArrayBox>& a_gradPhi, int a_dir) const;

  /// returns all components of grad(phi) (gradPhi should be correct size)
  void gradPhi(LevelData<FluxBox>& a_gradPhi) const;

  /// returns level-projection pressure
  LevelData<FArrayBox>& Pi();

  /// RHS for projection, i.e. \f$ \nabla \cdot \mathbf{U}^* \f$ (cell centred)
  LevelData<FArrayBox>& CCrhs();

  /// RHS for projection solve, i.e. \f$ \nabla \cdot \mathbf{U}^* \f$ (face centred)
  LevelData<FArrayBox>& MACrhs();

  /// Correction to make velocity incompressible (cell centred)
  LevelData<FArrayBox>& CCcorrection();

  /// Correction to make velocity incompressible (face centred)
  LevelData<FluxBox>& MACcorrection();

  /// Pressure, not scaled by the timestep
  void unscaledPi(LevelData<FArrayBox>& unscaledPi, Real a_dt);

  /// returns grad(pi) in direction dir
  void gradPi(LevelData<FArrayBox>& a_gradPi, int a_dir) const;

  /// returns grad(pi) in all directions into SpaceDim-dimensioned gradPi
  void gradPi(LevelData<FArrayBox>& a_gradPi) const;

  /// Applys BCs
  void gradPiBCs(LevelData<FArrayBox>& a_gradPi, bool extrapBCs=false, bool a_usePhi=false);

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

  /// correct face centred velocities with the already calculated grad(e_lambda)
    void applyFreestreamCorrection(LevelData<FluxBox>& a_advVel, Real scale=1.0);

    /// Remove freestream correction from face centred velocities
    void removeFreestreamCorrection(LevelData<FluxBox>& a_advVel);

  /// returns edge-centered grad(eLambda) in all directions
  /** non-constant reference returned because of need to rescale. */
  LevelData<FluxBox>&  grad_eLambda();

  /// returns edge-centered grad(eLambda) in all directions
  const LevelData<FluxBox>& grad_eLambda() const;

  /// should we do sync projection?
  bool doSyncProjection() const;

  /// should we do volume discrepancy correction?
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

  /// Set fine proj ptr
//  void setFineProjPtr(Projector* fineProj);

  /// returns verbosity
  int verbosity() const;

  /// sets physBCs
  void setPhysBC(PhysBCUtil& a_bc);

  /// pointer to the boundary condition object
  PhysBCUtil* getPhysBCPtr() const;

  /// Sum of RHS for VD correction
  Real m_sumVDrhs;

  /// Version which defines it's own amr multigrid
   void computeVDCorrection(Vector<LevelData<FArrayBox>* >& a_lambda,
                              Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosity,
                              const Real a_newTime,
                              const Real a_dtSync
                             );


protected:
  /// returns uEdge - G(phi_mac); C/F BC is CFscale*Pi(crse)
  /**
   * \f[ \mathbf{u}^{edge} - \nabla \phi^{mac} \f]
   */
  void applyMacCorrection(LevelData<FluxBox>& a_uEdge, LevelData<FArrayBox>* crseBCDataPtr, Real phiScale=1.0);

  /// vel := vel - scale*G_CC(pi)
  /**
   * \f[ \mathbf{u} = \mathbf{u} - \alpha \nabla \Pi \f]
   * where \f$ \alpha \f$ is the scale
   */
  void applyCCcorrection(LevelData<FArrayBox>& a_velocity,
                         LevelData<FArrayBox>* a_pressure,
                         LevelData<FArrayBox>* a_pressureScalePtr,
                           const Real scale); //const;

  /// do sync projection with this level as lBase; corrects velocity in place
  // assumes multigrid already set up
  void doSyncProjection(Vector<LevelData<FArrayBox>* >& a_velocity,
                        Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity,
                        const Real a_newTime, const Real a_dtSync,
                        AMRMultiGrid<LevelData<FArrayBox> >& a_solver
                        );

  /// perform initial sync projection
  void initialSyncProjection(Vector<LevelData<FArrayBox>* >& a_velocity,
                             Vector<RefCountedPtr< LevelData<FArrayBox> > >& a_porosity,
                             const Real a_newTime, const Real a_dtSync,
                             AMRMultiGrid<LevelData<FArrayBox> >& a_solver
                             );

  /// do volume discrepancy solve, with this level as lBase
  void computeVDCorrection(Vector<LevelData<FArrayBox>* >& a_lambda,
                           Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosity,
                           const Real a_newTime,
                           const Real a_dtSync,
                           AMRMultiGrid<LevelData<FArrayBox> >& a_solver
                           );



  /// Compute volume discrepancy with this level as lMax
//  void computeVDCorrection2(Vector<LevelData<FArrayBox>* >& a_lambda,
//                             const Real a_dtSync,
//                             AMRMultiGrid<LevelData<FArrayBox> >& a_solver
//                             );



  ///  apply sync correction vel = vel - scale*G^{comp} eSync
  /** (recursively applies correction to all finer levels) */
  void applySyncCorrection(Vector<LevelData<FArrayBox>* >& a_velocity,
                           Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity,
                           const Real scale,
                           LevelData<FArrayBox>* crseCorr);

  /// compute composite gradient of eLambda -- does this recursively
  void computeGrad_eLambda(Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosity);



  /// rescales composite gradient of eLambda, recursively all finer levels
  /** rescales grad(eLambda) by (a_dx_lbase)/m_dx for this level and
      recursively calls this function for finer levels.  This is done
      for stability reasons.
  */
//  void rescaleGrad_eLambda(Real a_dx_lbase);

  /// defines Multilevel AMRSolver used for composite projections
  /** uses grids in vel, lBase is this level.
      a_freestream solver is true if solver will be used for
      freestream preservation solve. */
  void defineMultiGrid(AMRMultiGrid<LevelData<FArrayBox> >& a_solver,
                       const Vector<LevelData<FArrayBox>* >& a_vel,
                       Vector<RefCountedPtr<LevelData<FluxBox> > >& a_porosity,
                       bool a_freestreamSolver);

  /// Define multigrid solver for freestream preservation solve
//  void defineMultiGridFreestream(AMRMultiGrid<LevelData<FArrayBox> >& a_solver,
//                                 LinearSolver<LevelData<FArrayBox> >* a_bottomSolver,
//                                 Vector<LevelData<FArrayBox>* > a_lambda);

  /// defines multigrid solver for this level
  void defineSolverMGlevel(const DisjointBoxLayout& a_grids,
                           const DisjointBoxLayout* a_crseGridsPtr);

  /// Set things like tolerance, number of smoothing steps etc.
  void setSolverParameters();

  /// Setup bototm solvers
  void makeBottomSolvers();

  /// define multigrid solver for this level
  /**
   * variable coefficient version to deal with variable porosity
   */
  void defineSolverMGLevelVCOp(const RefCountedPtr<LevelData<FArrayBox> > a_porosityPtr,
                               const RefCountedPtr<LevelData<FArrayBox> > a_crsePorosityPtr,
		const RefCountedPtr<LevelData<FluxBox> > a_porosityEdgePtr,
		const RefCountedPtr<LevelData<FluxBox> > a_crsePorosityEdgePtr,
		bool cellCentred = false,
		Real beta=-1);

  /// solve for pressure on this level
  int solveMGlevel(LevelData<FArrayBox>&         a_phi,
                    const LevelData<FArrayBox>*   a_phiCoarsePtr,
                    const LevelData<FArrayBox>&   a_rhs,
					const RefCountedPtr<LevelData<FArrayBox> > a_porosityPtr,
					const RefCountedPtr<LevelData<FArrayBox> > a_crsePressureScalePtr,
					const RefCountedPtr<LevelData<FluxBox> > a_porosityEdgePtr,
					const RefCountedPtr<LevelData<FluxBox> > a_crsePorosityEdgePtr,
					bool cellCentred = false);

  /// For when we don't have any pressure scale
  void solveMGlevel(LevelData<FArrayBox>&   a_phi,
                    const LevelData<FArrayBox>*   a_phiCoarsePtr,
                    const LevelData<FArrayBox>&   a_rhs,
                    bool cellCentred);

  /// takes care of setting ghost cell values after restarting
  void postRestart();



private:

  /// Whether to use pi or phi for the pressure BC in the MAC projection
  bool m_usePiAdvectionBCs;

  /// Whether or not to scale pressure with porosity
  /**
   * i.e, do we want to solve \f$ \nabla^2 p = \nabla \cdot \mathbf{u} \f$ or
   * \f$ \nabla \cdot \chi \nabla p = \nabla \cdot \mathbf{u} \f$?
   */
  bool m_scalePressureWithPorosity;

  /// Also scale sync corrections
  bool m_scaleSyncCorrection;

  /// level of refinement
  int m_level;

  /// cell spacing
  Real m_dx;

  /// problem domain
  ProblemDomain m_domain;

  /// Quadratic C/F interpolation object
  QuadCFInterp m_cfInterp;

  /// cache intvectsets for gradient here (if relevant)
  LayoutData<IntVectSet> m_gradIVS;

  /// ptr to finer-level Projector object (if it exists -- o/w == NULL)
  Projector* m_fineProjPtr;
  /// ptr to coarser-level Projector (if level == 0, ptr == NULL)
  Projector* m_crseProjPtr;

  /// if coarse level exists, refinement ratio down to it
  int m_nRefCrse;

  /// cell centered pressure data at time n+1/2
  LevelData<FArrayBox> m_Pi;

  /// MAC correction field
  LevelData<FArrayBox> m_phi;

  /// Divergence of advection velocity (MAC RHS)
  LevelData<FArrayBox> m_MACrhs;

  /// Divergence of CC velocity (CC RHS)
  LevelData<FArrayBox> m_CCrhs;

  LevelData<FArrayBox> m_CCcorrection;
  LevelData<FluxBox> m_MACcorrection;

  /// synchronization correction
  LevelData<FArrayBox> m_eSync;
  /// volume-discrepancy correction
  LevelData<FArrayBox> m_eLambda;

  /// edge-centered gradient of eLambda
  LevelData<FluxBox> m_grad_eLambda;

  /// For variable coefficient pressure solve
  Vector< RefCountedPtr<LevelData<FArrayBox> > > m_aCoef;
  Vector< RefCountedPtr<LevelData<FArrayBox> > > m_bCoefCC;

  /// For variable coefficient pressure solve
  Vector< RefCountedPtr<LevelData<FluxBox> > > m_bCoef;

  /// Porosity - used for pressure solve
  RefCountedPtr<LevelData<FArrayBox> > m_porosityPtr;
  RefCountedPtr<LevelData<FArrayBox> > m_crsePorosityPtr;

  /// Porosity - used for pressure solve (edge-centred)
  RefCountedPtr<LevelData<FluxBox> > m_porosityEdgePtr;

  /// Porosity on coarse level (edge centred)
  RefCountedPtr<LevelData<FluxBox> > m_crsePorosityEdgePtr;

  /// Multigrid solver on a level
  AMRMultiGrid<LevelData<FArrayBox> > m_solverMGlevel;

  /// All disjoint box layouts in the AMR hierarchy
  Vector<DisjointBoxLayout> m_allGrids;

  /// Should we do the synchronisation projection?
  static bool m_doSyncProjection;

  /// Should we apply the synchronisation correction
  static bool m_applySyncCorrection;

  /// Is object initialized
  bool m_isInitialized;

  /// Should we do quadratic interpolation at coarse fine boundaries?
  static bool m_doQuadInterp;

  /// solver parameter -- passed on through when defining solvers
  bool m_limitSolverCoarsening;

  /// \f$ \eta_\lambda \f$ is the coefficient for the freestream preservation correction
  /**
   * Should be between 0 and 1. Choosing 1 is probably unstable.
   * Choosing 0 is pointless (no correction is applied)
   */
  static Real m_etaLambda;

  /// multigrid solver tolerance
  static Real s_solver_tolerance;

  /// multigrid hang condition
  static Real s_solver_hang;

  /// multigrid solver parameter
  static int s_num_smooth_up;

  /// multigrid solver parameter
  static int s_num_smooth_down;

  /// number of smoothings in preconditioner for projection
  static int s_num_precond_smooth;

  static int s_numMG;

  /// Whether to use relax solver for bottom solve
  /**
   * Turned off by default. This works better for non 2^n grids though.
   * Generally want block factor >= max grid size/2 to accurately do MAC projection
   */
  static bool s_relax_bottom_solver;

  /// use consistent scaling for lamdba solve?
  static bool s_constantLambdaScaling;

  /// timestep used for all subcycled lamdba solves
  static Real s_lambda_timestep;

  /// have the static variables been initialized (ParmParse)
  static bool pp_init;

  /// verbosity of output
  static int s_verbosity;

  /// Max number of relax solves on the bottom level of the multigrid
  static int s_bottomSolveMaxIter;

  /// multigrid relaxation scheme
  static int s_multigrid_relaxation;

  /// if false, VD correction is set to zero once it's computed
  static bool s_applyVDCorrection;

  /// is this the finest extant level?
  bool m_finest_level;

  Real m_phiScale;
  Real m_MACbcScale;

  bool m_scale_lambda_with_porosity;
  bool m_scale_lambda_err_with_porosity;

  ///bottom solver for level solves BiCGStab works much better for cases where we don't have 2^N grid cells
//  RelaxSolver<LevelData<FArrayBox> > * m_bottomSolverLevel;
  RelaxSolver<LevelData<FArrayBox> > * m_bottomSolverLevel;
  BiCGStabSolver<LevelData<FArrayBox> > * m_BiCGBottomSolverLevel;

  ///
  PhysBCUtil* m_physBCPtr;
};

#endif


