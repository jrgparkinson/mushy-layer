
#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELMUSHYLAYER_H_
#define _AMRLEVELMUSHYLAYER_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "LevelFluxRegister.H"
#include "LevelDomainFluxRegister.h"

#include "PatchGodunov.H"
#include "PiecewiseLinearFillPatch.H"
#include "ProblemDomain.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "RealVect.H"
#include "LevelFluxRegister.H"
#include "DisjointBoxLayout.H"
#include "BCFunc.H"
#include "LevelTGA.H"
#include "AMRPoissonOp.H"
#include "BiCGStabSolver.H"
#include "RelaxSolver.H"
#include "AMRFASMultiGrid.H"
#include "Gradient.H"
#include "ParmParse.H"

#include "computeSum.H"
#include "computeNorm.H"
#include "SetValLevel.H"

#include "VCAMRPoissonOp2.H"
#include "AdvectionPhysics.H"
#include "VelBCHolder.H"

#include "CellToEdge.H"
#include "CellToEdge2.H"
#include "EdgeToCell.H"
#include "FourthOrderFineInterp.H"
#include "Divergence.H"
#include "PiecewiseLinearFillPatchFace.H"
#include "timeInterp.H"

// Files I made myself
#include "DarcyBrinkmanOp.H"
#include "MushyLayerParams.h"
#include "MushyLayerUtils.H"
#include "Diagnostics.h"
#include "PhysBCUtil.H"
#include "Projector.H"
#include "PhysIBC.H"
#include "phaseDiagram.H"
#include "AMRNonLinearMultiCompOp.H"
#include "mushyLayerOpt.h"

// Fortran files
#include "AdvectUtilF_F.H"
#include "AdvectPhysicsF_F.H"
#include "GodunovUtilitiesF_F.H"
#include "EnthalpyVariablesF_F.H"
#include "UtilF_F.H"

#include "NamespaceHeader.H"

#ifdef CH_USE_DOUBLE
#define TIME_EPS 1.0e-10
#else
#define TIME_EPS 1.0e-5
#endif

/// AMRLevel for mushy layer calculations
/**
 * Big class to manage integration of nonlinear mushy layer equations on an AMR level
 */
class AMRLevelMushyLayer : public AMRLevel
{
public:

  /// Default constructor
  AMRLevelMushyLayer()
{
    m_isDefined = false;
    m_usePrevPressureForUStar = true;

    m_domainWidth = -1.0;
    m_domainHeight = -1.0;

    m_dx = -1.0;
    m_hasFiner = false;
    m_hasCoarser = false;
    m_max_dt_growth = 1.0;
    m_numGhost = 0;
    m_numGhostAdvection = 0;
    m_computedCFL = 0;
    m_timestepFailed = false;
    m_timestepReduced = false;
    m_useLimiting = false;

    m_physBCPtr = NULL;

//    m_newLevel=false;

    m_enforceGradP = false;
    m_dtReduction = 0.0;

    m_numOutputComps = 0;
    m_newGrids_different = false;
    m_regrid_smoothing_done = false;
    m_adv_vel_centering = 0.5;
    m_adv_vel_centering_growth = 1.01;

    s_regrid_smoothing_coeff = 0.5;
    s_implicit_reflux = true;

    AMREnthalpySum_new = 0;
    AMREnthalpySum_old = 0;
    AMRSaltSum_new = 0;
    AMRSaltSum_old = 0;
    m_pressureScaleVar = -1;

    m_doAutomaticRestart = false;

    m_maxLambda=0.0;

}

  /// Full constructor. Arguments are same as in define()
  AMRLevelMushyLayer(MushyLayerOptions a_opt, MushyLayerParams a_params)
  {
    define(a_opt, a_params);
  }

  /// Defines this AMRLevelMushyLayer
  void define(MushyLayerOptions a_opt, MushyLayerParams a_params);

  /// Destructor
  virtual ~AMRLevelMushyLayer();

  /// Never called: historical
  virtual void define(AMRLevel*  a_coarserLevelPtr,
                      const Box& a_problemDomain,
                      int        a_level,
                      int        a_refRatio)
  {
    MayDay::Error("never called--historical");
  }

  /// Set default values
  void setDefaults();

  /// Define new AMRLevelMushyLayer from coarser
  virtual void define(
      AMRLevel*            a_coarserLevelPtr,
      const ProblemDomain& a_problemDomain,
      int                  a_level,
      int                  a_refRatio);

  /// Advance by one timestep
  virtual Real advance();


  /// Things to do after a timestep
  virtual void postTimeStep();

  /// Explicitly reflux scalar field
  void doExplicitReflux(int a_var);

  /// Reflux enthalpy and bulk concentration
  Real doHCreflux();

  /// Do reflux on momentum field
  void doMomentumReflux(Vector<LevelData<FArrayBox>*>& compVel);

  /// Create tags for regridding
  virtual void tagCells(IntVectSet& a_tags) ;

  /// Tags cells on the mush side of the mush-liquid boundary
  void tagMushLiquidBoundary(IntVectSet& localTags);

  /// Tag cells on the domain boundary
  void tagBoundaryLayerCells(IntVectSet& localTags);

  /// Tag cells inside grid
  void tagCenterCells(IntVectSet& localTags,  Real radius = 0);

  /// Tag cells with porosity < 1
  void tagMushyCells(IntVectSet& localTags);

  /// Add cells that specify criteria
  void tagCellsVar(IntVectSet& localTags, Real refineThresh,
                   int taggingVar, int vectorTaggingVar,
                   TaggingMethod taggingMethod, int comp=-1);

  /// Setup menagerie of data structures
  void levelSetup();

  /// Setup CF interp objects
  void defineCFInterp();

  /// Set eta on all projection operators
  /**
   * Only call from level 0
   */
  void setEta(Real a_eta);

  /// Maximum allowed eta
  Real maxAllowedEta();

  /// Set adv_vel_centering on all levels
  void setAdvVelCentering(Real a_fraction);

  /// Create tags at initialization
  virtual void tagCellsInit(IntVectSet& a_tags) ;

  /// Set up data on this level after regridding
  virtual void regrid(const Vector<Box>& a_newGrids);

  /// Initialize grids
  virtual void initialGrid(const Vector<Box>& a_newGrids);

  /// Operations to do  after initial grids have been set up
  virtual void postInitialGrid(const bool a_restart);

  /// Operations to execute after we've done regridding
  virtual void postRegrid(int a_base_level);

  /// Create structures to hold data
  void createDataStructures();

  /// Initialize data structures - resize arrays etc.
  void initDataStructures();

  /// Compute cell-centred and face-centred velocities
  void computeAllVelocities(bool doFRupdates);

  /// Initialize data
  virtual void initialData();

  /// Initial data for corner flow problem
  void initialDataCornerFlow();

  /// Initial data for  Horton-Rogers-Lapwood problem
  void initialDataHRL();

  /// Initial data for convection in a hetrogeneous porous medium
  void initialDataConvectionMixedPorous();

  /// Initial data for poiseuille flow
  void initialDataPoiseuille();

  /// Initial data for a solid ice block
  void initialDataIceBlock();

  /// Initial data for mushy layer simulations
  void initialDataMushyLayer();

  /// Mushy layer with a porous hole in the middle
  void initialDataPorousHole();

  /// Initial data for a test case with regions of zero porosity
  void initialDataZeroPorosityTest();

  /// Initial data for diffusion problem
  void initialDataDiffusion();

  /// Initial data for vortex pair problem
  void initialDataVortexPair();

  /// Add vortex to domain
  void addVortex(RealVect center, Real strength, Real radius);

  /// Initial data for Rayleigh-Benard problem
  void initialDataRayleighBenard();

  /// Initial data for solute flux test problem
  void initialDataSoluteFlux();

  /// Initial data for sidewall heating problem
  void initialDataSidewallHeating();

  /// Default initial data
  void initialDataDefault();

  /// Initial data for burgers problem
  void initialDataBurgers();

  /// Initial data for problem to test refluxing
  void initialDataRefluxTest();

  /// Whether or not we're loading advection velocities
  bool loadAdvVel() { return m_opt.load_advVel; }

  /// Have we converged to steady state?
  virtual bool convergedToSteadyState();

  /// Things to do after initialization
  virtual void postInitialize();

#ifdef CH_USE_HDF5
  /// Write checkpoint header
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// Write checkpoint data for this level
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// Read checkpoint header
  virtual void readCheckpointHeader(HDF5Handle& a_handle);

  /// Read checkpoint data for this level
  virtual void readCheckpointLevel(HDF5Handle& a_handle);

  /// Write plotfile header
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  /// Write plotfile data for this level
  virtual void writePlotLevel(HDF5Handle& a_handle) const;

  void writePlotFile(int iter);
#endif

  /// Calculate vorticity and streamfunction so we can plot them
  void computeVorticityStreamfunction();

  /// Compute vorticity
  void computeVorticity();

  /// Get \f$ \pi \f$ and \f$  \phi \f$ from projection so we can plot them
  void getExtraPlotFields();

  /// Write AMR hierarchy to file
  void writeAMRHierarchy(string filename);

  /// Set reference point for nondimensionalisation to be the initial temperature/salinity
  void setDimensionlessReferenceInitial();

  /// Set reference point for nondimensionalisation to be the eutectic temperature/salinity
  void setDimensionlessReferenceEutectic();

  /// Returns the dt computed earlier for this level
  virtual Real computeDt();

  /// Compute dt using initial data
  virtual Real computeInitialDt();

  /// Compute nusselt number, salt fluxes etc.
  void computeDiagnostics();

  /// Define solvers
  void defineSolvers(Real a_time);

  /// Compute chimney diagnostics
  void computeChimneyDiagnostics();

  /// Diagnostics e.g. Nusselt number, solute flux
  Diagnostics m_diagnostics;

  /// Whether or not to compute diagnostics
//  bool m_computeDiagnostics;

  /// Get the grids on this level
  DisjointBoxLayout grids();

  /// Shift data in domain
  void shiftData(int dir, int distance);

  /// Set \f$ \Delta x \f$ on this level
  void dx(Real newDx);

  /// Get \f$ \Delta x \f$ on this level
  Real dx();

  /// Transfer data to new disjoint box layout
  void reshapeData(DisjointBoxLayout newGrids,  ProblemDomain newDomain);

  /// Refine data
  void refine(Real ref_ratio, DisjointBoxLayout a_grids, ProblemDomain a_domain);

  /// Do smoothing after regridding
  void doPostRegridSmoothing(bool a_smoothVel=true,
                             bool a_smoothScalar=true);

  /// Set the coefficient for smoothing fields
  void setSmoothingCoeff(Real a_coeff);

  /// Add a melt pond at the top of the domain
  void addMeltPond(int depth, Real salinity, Real enthalpy, bool rescaleExistingSolution=false);

  /// Make sure zero porosity regions have no velocity (cell centred velocity)
  void setVelZero(LevelData<FArrayBox>& a_vel, Real a_limit=-1, int a_radius=0);

  /// Make sure zero porosity regions have no velocity (face centred velocity)
  void setVelZero(LevelData<FluxBox>& a_vel, Real a_limit=-1);

  /// Set velocity to zero when porosity is below specified limit (cell centred)
  void setCCVelZero(Real a_limit);

  /// Update thermodynamic variables based on phase diagram
  void updateEnthalpyVariables();

  /// Compute lambda/porosity
  void computeLambdaPorosity();

  void set_compute_diagnostics(bool compute_diags);

protected:

  /// Has a particular variable converged to steady state?
  Real convergedToSteadyState(const int a_var, bool vector=false);

  /// Compute \f$ \partial / \partial t \f$ of some field
  void compute_d_dt(const int a_var, LevelData<FArrayBox>& diff, bool vector = false);

  /// Determine if we're solving the full time dependent equation, or just darcy's law
  bool solvingFullDarcyBrinkman();

  /// Update thermodynamic variables based on phase diagram (old time)
  void updateEnthalpyVariablesOld();

  /// Returns true if there is not a finer level with a grid defined on it
  bool finestLevel();

  /// Return the finest level
  int getFinestLevel();

  /// Returns \f$ \Delta t \f$
  Real computeDt(Real cfl);

  /// Get largest velocity in domain
  Real getMaxVelocity();

  /// get largest advection velocity, considering \f$ \mathbf{U}/\chi \f$ if appropriate
  Real getMaxVelocityForCFL();

  /// Compute max \f$ \mathbf{U}/\chi \f$
  Real computeMaxUChi();

  /// Get max advection velocity
  Real getMaxAdvVel();

  /// Returns \f$ \Delta t \f$
  Real computeDt(bool growdt);



  /// Replace with data from before the last backup
  void restartTimestepFromBackup(bool ignorePressure = false);

  /// Save the current data
  void backupTimestep();

  /// Get information on the entire AMR hierarchy
  void getHierarchyAndGrids(Vector<AMRLevelMushyLayer*>&        a_hierarchy,
                            Vector<DisjointBoxLayout>&             a_grids,
                            Vector<int>&                           a_refRat,
                            ProblemDomain&                         a_lev0Dom,
                            Real&                                  a_lev0Dx);

  /// Calculate permeability \f$ \Pi \f$
  void calculatePermeability();

  /// Fills coarse-fine boundaries on variable by quadratic interpolation
  void calculateCoarseFineBoundaries(int a_var, bool vector=false);

  /// Advection diffusion for multiple components
  int multiCompAdvectDiffuse(LevelData<FArrayBox>& a_phi_old, LevelData<FArrayBox>& a_phi_new,
                             LevelData<FArrayBox>& a_src,
                             bool doFRupdates = false, bool computeAdvectiveSrc = true);

  /// Compute src term and add it to a_src
  void computeScalarAdvectiveSrc(LevelData<FArrayBox>& a_src, int a_var,  bool converged, int a_comp = 0);

  /// Compute multi component source term: \f$ (\mathbf{U} \cdot \nabla T, \mathbf{U} \cdot \nabla S_l)\f$
  void computeScalarAdvectiveSrcHC(LevelData<FArrayBox>& a_src,
                                   LevelData<FluxBox>& edgeScalTotal,
                                   bool converged);

  /// Increment Enthalpy-Bulk concentration flux registers with flux
  void incrementHCFluxRegisters(LevelData<FluxBox>& flux, Real fluxMult);

  /// Compute total fluxes that provide the source term for implicit enthalpy-salinity updates
  /**
   * \f$ (\mathbf{U} T + V H \mathbf{k}, \mathbf{U} S_l + V S \mathbf{k})\f$
   */
  void computeTotalAdvectiveFluxes(LevelData<FluxBox>& edgeScalTotal);

  /// Update a scalar field by advection
  /**
   * a_scalarVar is the field which is updated, whilst a_advectionVar is the field that
   * has been advected. These don't need to be the same, i.e. could do
   *    advectScalar(enthalpy, temperature, ...)
   * in order to compute an update for the equation d(enthalpy)/dt + u.grad(temperature) = 0
   */
  void advectScalar(const int a_scalarVar, const int a_advectionVar,
                    LevelData<FluxBox>& a_advVel, bool doFRupdates = true);

  /// Update a scalar field with advective terms
  void advectScalar(const int a_scalarVar, const int a_advectionVar,
                    LevelData<FluxBox>& a_advVel, bool doFRupdates,
                    LevelData<FluxBox>& flux);

  /// Update flux register for scalar field
  void updateScalarFluxRegister(int a_scalarVar, LevelData<FluxBox>& flux, Real scale);

  /// Advect lambda field (for freestream preservation)
  void advectLambda(bool doFRupdates = true);


  /// Make source term for the advection velocity solve
  void computeAdvectionVelSourceTerm(LevelData<FArrayBox>& advectionSourceTerm);

  /// Average vertical fluxes in the horizontal direction
  void horizontallyAverage(LevelData<FArrayBox>& a_averaged, LevelData<FluxBox>& a_phi);

  /// Average vertical fluxes in the horizontal direction
  void horizontallyAverage(Vector<Real>& averageVector, LevelData<FluxBox>& a_phi);

  /// Average vertical fluxes in the horizontal direction
  void horizontallyAverage(LevelData<FArrayBox>& a_averaged, LevelData<FluxBox>& a_phi, Vector<Real>& averageVector);

  /// Average cell centred fields in the horizontal direction
  void horizontallyAverage(LevelData<FArrayBox>& a_averaged, LevelData<FArrayBox>& a_phi);

  /// Average cell centred fields in the horizontal direction
  void horizontallyAverage(Vector<Real>& averageVector, LevelData<FArrayBox>& a_phi);

  /// Average cell centred values in the horizontal direction
  void horizontallyAverage(LevelData<FArrayBox>& a_averaged, LevelData<FArrayBox>& a_phi,
                           Vector<Real>& globalAveraged);

  /// Compute average value of some scalar field over the liquid region of the domain
  Real averageOverLiquidRegion(int a_var);

  /// Calculate advective + diffusive flux
  void getTotalFlux(LevelData<FluxBox>& totalFlux);

  /// Compute \f$ \nabla^2 \mathbf{U} \f$
  void computeLapVel(LevelData<FArrayBox>& a_lapVel,
                     LevelData<FArrayBox>& a_vel,
                     const LevelData<FArrayBox>* a_crseVelPtr);

  /// Compute \f$ \nabla^2 \phi \f$
  void computeScalDiffusion(const int a_var,
                            LevelData<FArrayBox>& a_lapScal,
                            Real a_time);



  /// Compute \f$ \nabla^2 (H, S) \f$
  void computeScalDiffusion(LevelData<FArrayBox>& diffusiveSrc,
                            Real a_time);

  /// Compute \f$ \nabla^2 \psi \f$ for some scalar field \f$ \psi \f$
  void computeScalDiffusion(LevelData<FArrayBox>& a_src, int a_var);

  /// Compute new cell centred velocity
  /**
   * Generally done via an implicit update due to the viscous and darcy terms.
   * Not required if just considering the Darcy equation.
   */
  void computeCCvelocity(const LevelData<FArrayBox>& advectionSourceTerm, Real a_oldTime, Real a_dt,
                         bool doFRupdates = false,
                         bool a_doProjection = true,
                         bool compute_uDelU = true,
                         bool a_MACprojection = false);

  /// Compute \f$ \mathbf{U}^* \f$, the velocity before we project it
  void computeUstar(LevelData<FArrayBox>& a_uDelu,
                    const LevelData<FArrayBox>& advectionSourceTerm,
                    Real a_oldTime, Real a_dt, bool doFRupdates, bool a_MACprojection = false, bool compute_uDelU=true);

  /// Compute the source term for the \f$\mathbf{U}^* \f$ solve
  void computeUstarSrc(LevelData<FArrayBox>& src,
                       const LevelData<FArrayBox>& advectionSourceTerm,
                       Real src_time,
                       bool a_MACprojection = false, bool compute_uDelU=true);

  /// Compute \f$ \mathbf{u} \cdot \nabla \left\( \mathbf{U}/\chi \right\)
  void computeUDelU(LevelData<FArrayBox>& U_adv_src, const LevelData<FArrayBox>& advectionSourceTerm, Real a_oldTime, Real a_dt);

  /// Check CFL based on latest velocity is safe to use for explicit updates
  bool currentCFLIsSafe(bool printWarning=false);

  /// Set velocity to zero where porosity is sufficiently small (FArrayBox version)
  void setVelZero(FArrayBox& a_vel, const FArrayBox& a_porosity, const Real a_limit, const int a_radius=0);

  /// Returns \f$ \nabla P \f$ at the specified time
  void computeGradP(LevelData<FArrayBox>& gradP, Real a_time,
                    bool a_macProjection = false);

  /// Trace the advection velocity forward in time
  void
  traceAdvectionVel(LevelData<FluxBox>& a_advVel,
                    LevelData<FArrayBox>& a_old_vel,
                    LevelData<FArrayBox>& U_chi,
                    const LevelData<FArrayBox>& a_viscousSource,
                    PatchGodunov& a_patchGodVelocity,
                    Real a_old_time, Real a_dt);

  /// Compute advection velocities at the half time step
  void computePredictedVelocities(LevelData<FluxBox>& U_chi_new,
                                  LevelData<FArrayBox>& a_traceVel,
                                  LevelData<FluxBox>& a_advVel,
                                  LevelData<FArrayBox>& U_chi,
                                  const LevelData<FArrayBox>& a_viscousSource,
                                  PatchGodunov& a_patchGodVelocity,
                                  LevelData<FluxBox>& a_grad_eLambda,
                                  LevelData<FluxBox>& a_gradPhi,
                                  LevelData<FluxBox>& porosityFace,
                                  Real a_old_time, Real a_dt);

  /// Compute the advective flux given the scalar field, advection velocity and source term specified
  /**
   * Should work for multiple components
   */
  void
  computeScalarAdvectiveFlux(LevelData<FluxBox>& a_edgeScal,
                             LevelData<FArrayBox>& a_old_scal,
                             LevelData<FluxBox>& a_adv_vel,
                             LevelData<FluxBox>& a_inflowOutflowVel,
                             LevelData<FArrayBox>& a_old_vel,
                             LevelData<FArrayBox>& a_diffusiveSrc,
                             PatchGodunov& a_patchGod,
                             Real a_old_time, Real a_dt);


  /// Refactored this so we can also use for velocities if we want
  void upwind(LevelData<FluxBox>& a_edgeScal,
              LevelData<FArrayBox>& a_old_scal,
              LevelData<FluxBox>& a_adv_vel,
              LevelData<FluxBox>& a_inflowOutflowVel,
              LevelData<FArrayBox>& a_old_vel,
              LevelData<FArrayBox>& a_diffusiveSrc,
              PatchGodunov& a_patchGod,
              Real a_old_time, Real a_dt);

  /// Compute the advective flux
  /**
   * advection var is what we advect
   * diffusion var determines the variable used to calculate the diffusive src. Set to -1 to use no source.
   */
  void computeScalarAdvectiveFlux(LevelData<FluxBox>& a_edgeScal, int a_advectionVar, int a_diffusionVar,
                                  LevelData<FluxBox>& a_advVel,
                                  Real a_old_time, Real a_dt);

  /// Compute the advective flux for the enthalpy and salinity equations,
  /**
   * edgeScal contains two components: \f$ (\mathbf{U} T, \mathbf{U} S_l) \f$
   * or \f$ (V H, V S) \f$ depending on the parameters given
   */
  void computeScalarAdvectiveFluxMultiComp(LevelData<FluxBox>& a_edgeScal,
                                           LevelData<FluxBox>& a_advVel,
                                           PatchGodunov& a_patchGod,
                                           LevelData<FArrayBox>& a_scalOld,
                                           Real a_old_time, Real a_dt);

  /// Calculate  \f$ \mathbf{u} \cdot \nabla \left( \mathbf{u}/\chi \right) \f$
  void predictVelocities(LevelData<FArrayBox>& a_uDelU,
                         LevelData<FluxBox>& a_advVel,
                         const LevelData<FArrayBox>& a_src,
                         Real old_time,
                         Real a_dt,
                         bool doFRupdates);


  /// Calculate advection velocity for momentum equations which are time independent
  void calculateTimeIndAdvectionVel(Real a_time, LevelData<FluxBox>& a_advVel);
  void fillUnprojectedDarcyVelocity(LevelData<FluxBox>& a_advVel, Real time);

  /// Get maximum level allowed in AMR hierarchy
  int getMaxLevel();

  void computeViscosity();

  /// Fill advection velocity ghost cells
  void fillAdvVel(Real time, LevelData<FluxBox>& a_advVel);

  /// Fill the velocity with some analytically determined function
  void fillAnalyticVel(LevelData<FluxBox>& a_advVel);

  /// Fill the velocity with some analytically determined function
  void fillAnalyticVel (FArrayBox& velDir, int dir, int comp, bool project);

  /// Fill the velocity with some analytically determined function
  void fillAnalyticVel(LevelData<FArrayBox>& velDir);

  /// Compute the buoyancy force
  void fillBuoyancy(LevelData<FArrayBox>& a_buoyancy, Real a_time);

  /// Compute the buoyancy force
  void fillBuoyancy(LevelData<FArrayBox>& a_buoyancy,
                    LevelData<FArrayBox>& a_temperature,
                    LevelData<FArrayBox>& a_liquidConc,
                    LevelData<FArrayBox>& a_porosity);

  /// Compute the buoyancy force
  /**
   * \f$ \mathbf{B} = \chi ( Ra_T \theta  - Ra_C \Theta_l + \mathbf{F} \cdot \mathbf{k} ) \mathbf{k} \f$
   */
  void fillBuoyancy(FArrayBox& buoyancy,FArrayBox& temperature, FArrayBox& liquidConc,
                    FArrayBox& porosity,
                    FArrayBox& bodyForce);


  /// Compute \f$ \nabla P \f$ at time \f$ t \f$
  void fillPressureSrcTerm(LevelData<FArrayBox>& gradP,
                           LevelData<FArrayBox>& pressureScale,
                           Real a_time,
                           bool a_MACprojection);

  /// Compile velocity and porosity across an AMR hierarchy
  void fillAMRVelPorosity(Vector<LevelData<FArrayBox>*> & amrVel,
                          Vector<RefCountedPtr< LevelData<FluxBox> > > &  amrPorosityFace,
                          Vector<RefCountedPtr< LevelData<FArrayBox> > >& amrPorosity);

  /// Should include correct exchange calls
  void fillHC(LevelData<FArrayBox>& a_phi, Real a_time,
              bool doInterior = true, bool quadInterp = true);

  /// Fill porosity with some fixed value if appropriate
  void fillFixedPorosity(LevelData<FArrayBox>& a_porosity);

  /// Should include correct exchange calls
  void fillTCl(LevelData<FArrayBox>& a_phi, Real a_time,
               bool doInterior = true, bool quadInterp = true);

  /// Should include correct exchange calls
  void fillMultiComp(LevelData<FArrayBox>& a_phi, Real a_time, int scal1, int scal2,
                     bool doInterior = true, bool quadInterp = true);

  /// Compile \f$ \lambda \f$ across the entire AMR hierarchy
  void fillAMRLambda(Vector<LevelData<FArrayBox>*>& amrLambda);

  /// Calculate analytic solution to this problem (if possible)
  void calculateAnalyticSolns(bool enforceSolutions = true);

  /// Set fixed porosity
  void setPorosity(FArrayBox& a_porosity);

  /// Fill frame advection velocity
  void fillFrameVelocity();

  /// Add a small perturbation to the solution
  void addPerturbation(int a_var, Real alpha, int waveNumber = -1, Real phaseShift=0);

  /// Add a melt pond if we want
  void addMeltPond();


  /// Define implicit solver for \f$ \mathbf{u}^* \f$ including timestepping
  void defineUstarSolver(Vector<RefCountedPtr<LevelBackwardEuler> >&  a_UstarBE,
                         Vector<RefCountedPtr<LevelTGA> >& a_UstarTGA);

  /// Define multigrid solver for \f$ \mathbf{u}^* \f$
  void defineUstarMultigrid();

  /// Initialize pressure on this level
  void initializeLevelPressure(Real a_currentTime, Real a_dtInit);

  /// Initialize pressure across all levels
  void initializeGlobalPressure(Real dtInit=-1, bool init=true);

  /// Compute initial advection velocities
  void computeInitAdvectionVel();

  /// Initialise phi (pressure for Darcy equation)
  void initTimeIndependentPressure(AMRLevelMushyLayer* lev);

  /// Flag if a new level has been added
//  bool newLevelAdded();

  /// Reset \f$ \lambda = 1\f$
  void resetLambda();

  /// Reflux \f$ \lambda \f$ across AMR hierarchy
  /**
   * Only call from base level
   */
  void AMRRefluxLambda();


  /// Define operator for post regrid smoothing
  void defineRegridAMROp(AMRPoissonOpFactory& a_factory,
                         const Vector<DisjointBoxLayout>& a_grids,
                         const Vector<ProblemDomain>& a_domains,
                         const Vector<Real>& a_amrDx,
                         const Vector<int>& a_refRatios,
                         const int& a_lBase);

  /// Smooth the velocity field
//  void smoothVelocityField(int a_lbase);

  /// Utility function to get pointers to coarse level objects
  void getCoarseScalarDataPointers(const int a_scalarVar,
                                   LevelData<FArrayBox>** a_coarserDataOldPtr,
                                   LevelData<FArrayBox>** a_coarserDataNewPtr,
                                   LevelFluxRegister**    a_coarserFRPtr,
                                   LevelFluxRegister**    a_finerFRPtr,
                                   Real& a_tCoarserOld,
                                   Real& a_tCoarserNew);


  /// Doesn't do any time interpolation - just fills ghost cells if necessary
  void fillVectorField(LevelData<FArrayBox>& a_vector, Real a_time, int a_var,
                       bool doInterior=false, bool quadInterp=false);


  /// Fill ghost cells of scalar field
  void fillScalars(LevelData<FArrayBox>& a_scal,
                   Real a_time, const int a_var,
                   bool doInterior = false, bool quadInterp = false);

  /// Helper function
  void fillBuoyancy(LevelData<FArrayBox>& a_buoyancy, Real a_time, bool quadInterp);


  /// Fill ghost cells of scalar field (edge-centred)
  void fillScalarFace(LevelData<FluxBox>& a_scal,
                      Real a_time, const int a_var, bool doInterior = false, bool quadInterp = false);

  /// Fill ghost cells of scalar field (edge-centred)
  void fillScalarFace(LevelData<FluxBox>& a_scal,
                      Real a_time, const int a_var, CellToEdgeAveragingMethod method,
                      bool doInterior = false,
                      bool quadInterp = false,
                      Real smoothing = 0.0);

  /// Smooth some scalar field
  /**
   * smoothing done implicitly via an elliptic solve
   */
  void smoothScalarField(LevelData<FArrayBox>& a_phi, int a_var, Real a_smoothing);

  /// Get boundary conditions for implicit problems
  void getScalarBCs(BCHolder& thisBC, int a_var, bool a_homogeneous);

  /// Refactored this to be consistent
  void computeInflowOutflowAdvVel();

  /// Get initial and boundary conditions for advection
  PhysIBC* getScalarIBCs(int a_var);

  /// Need to convert implicit dirichlet/neuman/inflow outflow BCs to advection BCs
  /*
   * For advection we're calculating div (flux) = div (u phi)
   * dirichlet -> enforce value of flux at boundary, as the flux hasn't been
   *              multiplied by u yet so is just phi.
   * neumann -> if we want normal derivative = 0, this just means flux = 0
   *            NON-ZERO DERIVATIVES NOT CURRENTLY IMPLEMENTED
   * inflow outflow -> special function to deal with this.
   *
   */
  void convertBCType(const int a_implicitBC,  const Real a_implicitVal,
                     int a_explicitBC, Real a_explicitVal);

  /// Convert implicit BCs
  int convertBCType(const int a_implicitBC);

  /// Forcing term for darcy-brinkman test
//  void stokesDarcyForcing(LevelData<FArrayBox>& T, Real time);

  /// Check the simulation hasn't crashed
  bool crashed();

  /// Operations to keep variables bounded correctly.
  /**
   * Porosity \f$ \chi > 0 \f$
   * Liquid concentration \f$ -1 <= \Theta_l <= 0 \f$
   */
  void doRegularisationOps(LevelData<FArrayBox>& a_scal,
                           int a_var);

  /// Operations to keep variables bounded correctly.
  void doRegularisationOps(LevelData<FluxBox>& a_scal, int a_var);

  /// Operations to keep variables bounded correctly.
  void doRegularisationOps(int a_var, FArrayBox& a_state);
  void doRegularisationOpsNew(int a_var, FArrayBox& a_state);



  /// Reset all flux registers
  void setFluxRegistersZero();

  /// Copy new data to old data holders
  void copyNewToOldStates();

  /// Copy old data to new data holders
//  void copyOldToNewStates();

  /// Compute advection velocity
  void computeAdvectionVelocities(LevelData<FArrayBox>& advectionSrc,  Real advVelCentering=0.5);

  /// Correct velocity (projection + volume discrepancy)
  void correctEdgeCentredVelocity(LevelData<FluxBox>& a_advVel, Real a_dt);

  /// Smoot flux in x-direction
//  void horizontallySmooth(LevelData<FluxBox>& a_flux);
  /// Compute and apply volume discrepancy correction (edge centred velocity)
//  void computeVolumeDiscrepancyCorrection(LevelData<FluxBox>& a_advVel);

  /// Compute initial \f$ \Delta t \f$
  Real computeDtInit(int finest_level);

  /// Get the next coarser level
  AMRLevelMushyLayer* getCoarserLevel() const;

  /// Get the coarsest level
  AMRLevelMushyLayer* getCoarsestLevel();

  /// Get the next finer level
  AMRLevelMushyLayer* getFinerLevel() const;

  /// Set BCs on advection physics classes
  void setAdvectionBCs();
  void
  defineIBCs ();

  /// Advection velocity
  LevelData<FluxBox>   m_advVel,

  /// Advection velocity at old time
  m_advVelOld,

  /// Advection velocity at new time
  m_advVelNew,

  /// Frame advection velocity
  m_frameAdvVel,

  /// Frame advection velocity + fluid advection velocity
  m_totalAdvVel;

  /// Total salt content in domain at new time
  Real AMRSaltSum_new,

  /// Total salt content in domain at old time
  AMRSaltSum_old,

  /// Total heat content in domain at new time
  AMREnthalpySum_new,

  /// Total heat content in domain at new time
  AMREnthalpySum_old;

  /// Salt flux at the top of the domain
  LevelData<FArrayBox> m_saltFluxTop,

  /// Salt flux at the bottom of the domain
  m_saltFluxBottom;

  /// Frame velocity (cell centred). Don't think we need this.
  LevelData<FArrayBox> m_frameVel;

  /// List of scalar vars to write out to HDF5 files
  Vector<int> m_outputScalarVars,

  /// List of vector vars to write out to HDF5 files
  m_outputVectorVars,

  /// List of vector vars to write to checkpoint files
  m_chkVectorVars,

  /// List of scalar vars to write to checkpoint files
  m_chkScalarVars;

  /// Number of component we write out to HDF5 files
  int m_numOutputComps;

  /// What to do when restarting after a failed timestep
  enum solverRestart {
    m_restartHalveDt,
    m_restartResetData
  };

  /// Whether or not to restart with a small perturbation after reaching steady state
  /*
   * Useful because the instability often isn't actually kicked off
   */
  bool m_doAutomaticRestart;

  /// Store max value of \f$ \lambda-1 \f$ from previous timestep
  /**
   * We may decide to increase \f$ \eta \f$ if we're not decreasing lambda quickly enough
   */
  Real m_maxLambda;

  /// Scalar fields at new time
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_scalarNew;

  /// Scalar fields at old time
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_scalarOld;

  /// Backed up scalar fields to be used for restarting
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_scalarRestart;

  /// Lagged \f$ \partial \chi / \partial t \f$
  LevelData<FArrayBox> m_dPorosity_dt;

  /// Vector fields at new time
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_vectorNew;

  /// Vector fields at old time
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_vectorOld;

  /// Change in vector fields
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_dVector;

  /// Backed up vector fields to be used for restarting
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_vectorRestart;

  /// Names of scalar variables
  Vector<string> m_scalarVarNames;

  /// Names of vector variables
  Vector<string> m_vectorVarNames;

  /// TGA integrators (2nd order in time)
  RefCountedPtr<LevelTGA> m_enthalpySalinityTGA;

  /// Backward Euler integrators (1st order in time)
  RefCountedPtr<LevelBackwardEuler> m_enthalpySalinityBE;

  /// AMRMultigrid operators for each scalar
  RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >      m_diffuseAMRMG[m_numScalarVars];

  /// FASAMRMultigrid operators for each scalar
  RefCountedPtr<AMRFASMultiGrid<LevelData<FArrayBox> > >      m_diffuseAMRFASMG[m_numScalarVars];

  /// Multicomponent FAS Multigrid solver (for enthalpy-bulk concentration solves)
  RefCountedPtr<AMRFASMultiGrid<LevelData<FArrayBox> > >   m_multiCompFASMG;

  /// Operator factories for scalar vars
  RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > m_HCOpFact;

  /// bottom solver for velocity solve multigrid
  static BiCGStabSolver<LevelData<FArrayBox> >  s_botSolverUStar;

  /// Bottom solve for enthalpy-bulk concentration multigrid
  static RelaxSolver<LevelData<FArrayBox> > s_botSolverHC;

  /// AMRMultigrid solver for unprojected velocity, \f$ \mathbf{u}^* \f$
  RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >      m_uStarAMRMG[SpaceDim];

  /// Operator factory for unprojected velocity, \f$ \mathbf{u}^* \f$
  RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > m_uStarOpFact[SpaceDim];

  /// Operator for computing components of \f$ \nabla^2 \mathbf{u} \f$
  RefCountedPtr<AMRPoissonOp> m_viscousOp[SpaceDim];

  /// Is object defined
  bool m_isDefined;

  /// Struct to contain all the options for this simulation
  /**
   * This is so we can easily pass all the options around,
   * and don't have to keep going into parmparse
   */
  MushyLayerOptions m_opt;

  /// The CFL number we're currently running at
  /**
   * computed from the current max velocity and the grid spacing.
   */
  Real m_computedCFL;

  /// Domain width
  Real                 m_domainWidth;

  /// Domain height
  Real                 m_domainHeight;

  /// Number of cells in each direction of domain
  //  int                  m_numCells[SpaceDim];
  IntVect               m_numCells;

  /// Maxmimum fractional growth in \f$ \Delta t \f$ between timesteps
  Real 				 m_max_dt_growth;

  /// Should we do slope limiting in advection solves?
  bool                 m_useLimiting;

  /// Cell spacing
  Real                 m_dx;

  /// Diffusion coefficients for different scalar fields
  Real m_scalarDiffusionCoeffs[m_numScalarVars];

  /// Physical parameters
  MushyLayerParams m_parameters;

  /// Should we enforce some analytic \f$ \nabla P \f$
  /**
   * Usually false, but can be useful to do this for benchmark problems
   */
  bool m_enforceGradP;

  /// Should we use previous \f$ \nabla P \f$ to calculate \f$ \mathbf{u}^* \f$
  /**
   * We then subtract off \f$ \nabla P \f$  before actually doing the projection
   */
  bool m_usePrevPressureForUStar;

  /// Which field scales the pressure
  /**
   * We either have  \f$ \chi \nabla P \f$ or  \f$ \Pi \nabla P \f$
   */
  int m_pressureScaleVar;


  /// Pointer to physics boundary condition object
  PhysBCUtil* m_physBCPtr;

  /// For interpolating scalars from coarser grids
  FineInterp m_fineInterpScalar,

  /// For interpolating vectors from coarser grids
  m_fineInterpVector;

  /// For averaging scalars to coarse grids
  CoarseAverage m_coarseAverageScalar;

  /// For 2 component enthalpy-salinity fields
  CoarseAverage m_coarseAverageHC;

  /// For averaging vectors to coarse grids
  CoarseAverage m_coarseAverageVector;

  /// For reducing the timestep after a failed timestep
  bool m_timestepReduced;

  /// Factor reduce \f$ \Delta t \f$ by after a failed timestep
  int m_dtReduction;

  /// Has this timestep failed?
  bool m_timestepFailed;

  /// If this level has only just been added due to regridding
//  bool m_newLevel;

  /// The advection velocity is calculated at old_time + m_adv_vel_centering*dt
  /**
   * Usually 0.5, but set to something smaller after regridding for stability
   */
  Real m_adv_vel_centering;

  /// Rate at which we grow m_adv_vel_centering each timestep
  Real m_adv_vel_centering_growth;

  /// For doing quadratic interpolation at coarse fine boundaries (for scalars)
  QuadCFInterp m_quadCFInterpScalar,

  /// For doing quadratic interpolation at coarse fine boundaries (for vectors)
  m_quadCFInterpVector;

  /// Scalar fill patch for one ghost cell
  PiecewiseLinearFillPatch m_piecewiseLinearFillPatchScalarOne;

  PiecewiseLinearFillPatch m_piecewiseLinearFillPatchScalarTwo;
  PiecewiseLinearFillPatch m_piecewiseLinearFillPatchScalarThree;

  /// Scalar fill patch for four ghost cells
  PiecewiseLinearFillPatch m_piecewiseLinearFillPatchScalarFour;

  /// Number of ghost cells for standard calculations
  int m_numGhost,

  /// Number of ghost cells for advection calculations
  /**
   * Generally need extra ghost cells when doing advection problems
   */
  m_numGhostAdvection;

  /// Handles advection of velocity
  PatchGodunov m_patchGodVelocity,

  /// Operator For multi component H-C advection solve
  m_patchGodHC,

  /// Operator For multi component T-Sl advection solve
  m_patchGodTSl,

  /// Operator For single component H advection solve
  m_patchGodH,

  /// Operator For single component C advection solve
  m_patchGodC,

  /// Operator For single component T advection solve
  m_patchGodT,

  /// Operator For single component Sl advection solve
  m_patchGodSl;

  /// Physics for velocity advection
  AdvectionPhysics m_advectionPhysicsVelocity;

  /// Physics for coupled H-C advection solve
  AdvectionPhysics m_advPhysHC;

  /// Physics for coupled T-Sl advection solve
  AdvectionPhysics m_advPhysTSl;

  /// Physics for H advection solve
  AdvectionPhysics m_advPhysH;

  /// Physics for C advection solve
  AdvectionPhysics m_advPhysC;

  /// Physics for T advection solve
  AdvectionPhysics m_advPhysT;

  /// Physics for Sl advection solve
  AdvectionPhysics m_advPhysSl;

  /// Handles advection of scalar vars
  RefCountedPtr<PatchGodunov> m_patchGodScalars[m_numScalarVars];

  /// Initial and Boundary Conditions for scalar vars
  PhysIBC* m_scalarIBC[m_numScalarVars];

  /// Handles pressure projection
  Projector m_projection;

  /// Contains the old pressure, in case we need to restart
  Projector m_projectionBackup;

  /// Flux registers for scalar variables
  RefCountedPtr<LevelFluxRegister> m_fluxRegisters[m_numScalarVars],

  /// Flux Register for multi component HC solves
  m_fluxRegHC,

  /// Flux registers for vector variables
  m_vectorFluxRegisters[m_numVectorVars];

  /// Flux register for salt along the domain edge
  LevelDomainFluxRegister m_saltDomainFluxRegister,

  /// Flux register for heat (enthalpy) along the domain edge
  m_heatDomainFluxRegister;

  /// Which scalar vars should we bother making flux registers for?
  bool m_makeFluxRegForScalarVar[m_numScalarVars];

  /// Which vector vars should we bother making flux registers for?
  bool m_makeFluxRegForVectorVar[m_numVectorVars];

  /// Does a coarser AMRLevelMushyLayer pointer exist?
  bool m_hasCoarser;

  /// Does a finer AMRLevelMushyLayer pointer exist?
  bool m_hasFiner;

  /// Grids on this level
  DisjointBoxLayout m_grids;

  /// Are the new grids different from current ones?
  bool m_newGrids_different;

  /// Has post regrid smoothing been done?
  bool m_regrid_smoothing_done;

  /// Coefficient for post regrid smoothing
  Real s_regrid_smoothing_coeff;

  /// Do reflux implicitly?
  bool s_implicit_reflux;

private:

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelMushyLayer&);
  AMRLevelMushyLayer(const AMRLevelMushyLayer&);
};

#include "NamespaceFooter.H"

#endif
