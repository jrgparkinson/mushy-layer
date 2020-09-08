/*
 * mushyLayerOpt.h
 *
 *  Created on: 11 Feb 2019
 *      Author: parkinsonjl
 *
 *      This structure contains all the options for the mushy layer code which shouldn't change
 *      during a simulation
 */

#ifndef SRC_MUSHYLAYEROPT_H_
#define SRC_MUSHYLAYEROPT_H_

// Comment below ensures that enums doxygen generates documentation for enums declared in this file
/*!
@file
Contains immunity flag definitions
*/

/// Different ways to add refinement
enum RefinementMethod {
  /// Tag where fluid speed exceeds threshold
  tagSpeed,

  /// Tag where channels are
  tagChannels,

  /// Tag where these is strong downflow in the mushy layer
  tagPlumeMush,

  /// Tag porosity on level 0, and strong solute rich downflow thereafter
  tagMushChannels,

  /// Tag on some scalar field given by taggingVar
  tagScalar,

  /// Tag on some vector field given by  taggingVectorVar
  tagVector,

  /// Composite criteria for tagging cells.
  /**
   * Tag where \f$ \Delta x |\nabla \chi| ( 1 - \textrm{min}(w (\mathscr{C} - \Theta), 0)   ) > r \f$
   * where \f$ r \f$ is the refinement threshold, \f$ w \f$ is the vertical velocity,
   */
  tagMushChannelsCompositeCriteria,


};

/// Different multigrid methods
enum MGmethod {
  /// Normal linear multigrid
  MGTypeStandard,

  /// Full Approximation Scheme (for nonlinear problems)
  MGTypeFAS
};

/// Strategy for tagging cells
/**
 * Set by main.refineMethod
 */
enum TaggingMethod {
  /// 0 - Tag where the undivided gradient of some field is bigger than some value
  UndividedGradient = 0,

  /// 1 - Tag where the magnitude of some field is bigger than some values
  Magnitude,

  /// 2 - Tag where the value of some field is larger than some criteria
  CompareLargerThan,

  /// 3 - Tag where the value of some field multiplied by -1 is larger than some criteria.
  /**
   * Alternatively, can think of this as when some field is smaller than some criteria.
   */
  CompareLargerThanNegative
};

/// Methods for doing implicit refluxing - i.e. what sort of equation to we solve for the correction
enum RefluxMethod {
  /// linear correction
  LinearReflux,

  /// linear variable-coefficient correction
  LinearVCReflux,

  /// nonlinear correction
  NonlinearReflux
};

/// Options for handling advection velocities, where S is the source term
enum velocityAdvectionTypes
{
  /// solve \f$ \frac{\partial (\mathbf{u}/\chi)}{\partial t} + (u/porosity) \cdot \nabla (u/porosity) = S/porosity (-(u/porosity^2)dporosity/dt) \f$
  m_porosityInAdvection,

  /// solve du/dt + (u/porosity).grad(u) = S + (u.grad(porosity)/porosity^2)u
  m_porosityOutsideAdvection,

  /// solve du/dt + u.grad(u) = S
  m_noPorosity

};


/// Different options for enforcing a fixed porosity
enum PorosityFunctions
{
  ///  Porosity = 1.0 everywhere
  constantLiquid,

  /// Porosity varies linearly with x, \f$ \chi = 0.1 + 0.9 x \f$
  linear,

  /// Porosity looks like a gaussian, \f$ \chi =  0.1 + 0.9 exp(-100 (x-0.5)^2) \f$
  gaussian,

  /// Gaussian with some sinuosoidal variation, \f$ \chi = (exp(-(x-0.5)^2))*(0.5 + 0.1 sin(2 \pi y)) \f$
  gaussianSinusoidal,

  /// Porosity = 0.1
  constantSmall,

  /// Porosity varies with the y-direction cubed, \f$ \chi = 0.01 + (y-0.01)^3 \f$
  cubic,

  /// \f$ \chi = 0.5 ( (1-tanh(20 y)) + (1-x)) \f$
  hyperbolicTan,
};

/// Cell centred vector field
  enum VectorVars {
    /// Cell centred fluid velocity, \f$ \mathbf{U} \f$
    m_fluidVel,

    /// Fluid velocity divided by porosity, \f$ \mathbf{U}/\chi \f$
    m_U_porosity,

    /// Unprojected cell centred fluid velocity, \f$ \mathbf{U}^* \f$
    m_Ustar,

    /// Unprojected face centred fluid velocity, as used for advection, which has been averaged to cell centres
    /**
     * \f$ {Av}^{F \to C} \mathbf{U}_f \f$
     */
    m_advUstar,

    /// Face centred fluid velocity, as used for advection, which has been averaged to cell centres
    /**
     * \f$ {Av}^{F \to C} \mathbf{U}_{AD} \f$
     */
    m_advectionVel,

    ///
    m_viscousSolveSrc,

    /// Inertial term \f$ \mathbf{U} \cdot \nabla \left( \mathbf{U} / \chi \right)\f$
    m_UdelU,

    /// Source term for the velocity advection
    m_advectionSrc,
    m_fluidVelAnalytic,
    m_fluidVelErr,
    m_dUdt,
    m_FsDiffusion,
    m_FsFluid,
    m_Fs,
    m_freestreamCorrection,
    m_UpreProjection,
    m_advectionImplicitSrc,
    m_MACcorrection,
    CCcorrection,
    m_advSrcLapU,
    m_advUpreProjection,
    m_bodyForce,
    m_advVelCorr,

    /// Number of vector variables
    // Make sure this comes last!
    m_numVectorVars
  };

  /// Identifiers for different scalar variables
  enum ScalarVars {
    /// Enthalpy \f$ H \f$
    m_enthalpy,

    /// Bulk concentration \f$ \Theta \f$
    m_bulkConcentration,

    /// Temperature, \f$ \theta \f$
    m_temperature,

    /// Porosity \f$ \chi \f$
    m_porosity,

    /// Liquid concentration \f$ \Theta_l \f$
    m_liquidConcentration,

    /// Solid concentration \f$ \Theta_s \f$
    m_solidConcentration,

    /// Pressure used in the face centred projection \f$ \phi \f$
    m_pressure,

    /// Permeability \f$ \Pi \f$
    m_permeability,

    /// Viscosity \f$ \nu \f$
    m_viscosity,

    /// Auxillary field lambda \f$ \lambda \f$ for computing the freestream correction
    m_lambda,

    /// \f$ \lambda / \chi \f$
    m_lambda_porosity,

    /// Solidus phase boundary, \f$ H_S \f$
    m_enthalpySolidus,

    /// Liquidus phase boundary, \f$ H_L \f$
    m_enthalpyLiquidus,

    /// Eutectic phase boundary, \f$ H_E \f$
    m_enthalpyEutectic,

    /// Some analytically calculated porosity field \f$ \chi_{analytic} \f$
    m_porosityAnalytic,

    /// Some analytically calculated temperature field \f$ \theta_{analytic} \f$
    m_temperatureAnalytic,

    /// Explicit source term for the bulk concentration update
    /**
     * \f$ - \nabla \cdot \left( \mathbf{U} \Theta_l \right) - V \frac{\partial \Theta}{\partial z} \f$
     */
    m_saltEqnSrcGodunov,

    /// Computes \f$ \theta_{analytic} - \theta \f$ if required
    m_Terr,

    /// Explicit source term for the enthalpy update
    /**
     * \f$ - \nabla \cdot \left( \mathbf{U} \theta \right) - V \frac{\partial H}{\partial z} \f$
     */
    m_enthalpySrc,

    /// Level divergence of the advection velocities \f$ \mathbf{U}_{AD} \f$
    m_divUadv,

    /// Rate of change of enthalpy with time \f$ \frac{\partial H}{\partial t} \f$
    m_dHdt,

    /// Rate of change of bulk concentration with time \f$ \frac{\partial \Theta}{\partial t} \f$
    m_dSdt,

    /// Horizontally averaged vertical solute flux
    /**
     * \f$  \int_0^L \mathbf{F}_s \cdot \mathbf{z} dx \f$
     */
    m_averageVerticalFlux,

    /// Analytic solute flux for test problems
    m_soluteFluxAnalytic,

    /// Vertical component of the solute flux \f$ \mathbf{F}_s \cdot \mathbf{z} \f$
    m_verticalFlux,

    /// Level divergence of the cell centred velocity field \f$ \nabla \cdot \mathbf{U} \f$
    m_divU,

    /// Horizontally averaged vertical heat flux
    /**
     * \f$  \int_0^L \mathbf{F}_H \cdot \mathbf{z} dx \f$
     */
    m_averageHeatFlux,

    /// Streamfunction \f$ \psi \f$ (if calculated ad hoc during post processing)
    m_streamfunction,

    /// Vorticity \f$ \omega \f$ (if calculated ad hoc during post processing)
    m_vorticity,

    ///
    m_FsVertDiffusion,

    ///
    m_FsVertFluid,

    ///
    m_FsVertFrame,

    ///
    m_divUcorr,

    /// Pressure calculated in face-centred projection \f$ \pi \f$
    m_pi,

    /// Pressure calculated in cell-centred projection \f$ \phi \f$
    m_phi,

    /// Coarse-fine boundary condition used in the face-centred projection
    m_MACBC,

    /// Right hand side for the face-centred projection solve, \f$ \nabla \cdot \mathbf{U}_f \f$
    m_MACrhs,

    /// Right hand side for the face-centred projection solve, \f$ \nabla \cdot \mathbf{U} \f$
    m_CCrhs,

    /// Passive scalar advcted by flow
    m_passiveScalar,

    /// Active scalar which is advected and also has sources/sinks
    m_activeScalar,

    /// Intensity of radiation from the surface
    m_lightIntensity,

    /// Number of scalars variables
    // Make sure this comes last!
    m_numScalarVars

  };


//TODO: these should probably all be const's
//TODO: add comments to these options so it's obvious what they do
//TODO: collate more options in here rather than having lots of parmparses in the main code, so it's obvious what all the options are
/// Contains most of the options for running the code, and how to handle the equations
struct MushyLayerOptions {

    /// Output directory for plot and checkpoint files
  string output_dir;

  /// Prefix for plot files
  string plotfile_prefix;

  /// Turn on to only produce a minimal ammount of output
  bool minimalOutput;

  /// Turn on to produce lots of output
  bool debug;

  /// Maximum value by which the solution may have changed for us to have reached steady state
  Real steadyStateCondition;

  /// Turn on to ignore changes in velocity when deciding if we've reached steady state
  bool ignoreVelocitySteadyState;

  /// Turn on to ignore changes in bulk concentration when deciding if we've reached steady state
  bool ignoreBulkConcSteadyState;

  /// Domain width
  Real domainWidth;

  /// Domain height
  Real domainHeight;

  /// Maximum allowed CFL number
  Real cfl;

  /// TODO: remove this
  bool skipUnsafeCFL;

  /// Force the code to compute a CFL condition based on the size of \f$\mathbf{U}/\chi\f$, rather than \f$ \mathbf{U} \f$
  bool forceUseUChiForCFL;

  /// Minimum time to run simulations for (even if the steady state condition has been met)
  Real min_time;

  /// Maximum allowed fractional increase in dt from one timestep to the next
  Real max_dt_growth;

  /// Multiply the dt computed during initialisation procedures by this factor before using it
  /**
   * Useful for using a smaller dt initially for stability
   * TODO: can probably remove this, it's just a duplicate of initial_dt_multiplier really
   */
  Real init_dt_scale;

  /// Enforce a fixed timestep  (if > 0)
  Real fixedDt;

  /// Scale initial dt's by this factor
  Real initial_dt_multiplier;

  /// Minimum dt allowed
  Real minDt;

  /// Maximum dt allowed
  Real max_dt;

  /// Maximum initial dt allowed
  Real max_init_dt;

  /// CFL number for acceleration. Only used if useAccelDt = true
  /**
   * When velocities are changing rapidly due to a large acceleration, the timestep
   * computed based on the current velocities may not be a appropriate. Instead,
   * we should estimate the new velocities given the acceleration:
   *
   * \f$ \Delta t = \sqrt{\sigma  \Delta x / a} \f$,
   *
   * where \f$ \sigma \f$ is the cfl number and \f$ a\f$ is the maximum acceleration.
   *
   */
  Real accelCFL;

  /// Turn on to print some diagnostics about the dt computed according to acceleration considerations
  bool printAccelDt;

  /// Use the timesteps computed according to the cfl condition applied to accelerations
  bool useAccelDt;

  /// For the first timestep, use dt computed according to the cfl condition applied to accelerations
  bool useInitAccelDt;

  /// Maximum level allowed in this AMR hierarchy
  int max_possible_level;

  /// Turn on to compute the vorticity/streamfunction as a diagnostic
  bool computeVorticityStreamFunction;

  /// Turn on to use fortran routines for regularising the solution on cell faces (i.e. ensuring porosity is not 0)
  /**
   * Turned off by default as doesn't work yet.
   */
  bool useFortranRegularisationFace;

  /// Turn on using fortran routines for regularising the solution on cell centres (i.e. ensuring porosity is not 0)
  bool useFortranRegularisation;

//  Real stokesDarcyForcingTimescale;

  /// Whether or not we want to use viscous boundary conditions
  bool viscousBCs;

  /// Number of loops of the pressure initialisation routine
  int num_init_passes;

  /// If >= 0, and we're restarting, set the simulation time to this
  Real restart_new_time;

  /// Whether to use previous estimates of the pressure to compute the initial cell centred \f$ \mathbf{U}^* \f$
  /**
   * Turned off by default
   */
  bool init_use_prev_pressure_for_Ustar;


  /// Whether to compute \f$ \mathbf{U} \cdot \nabla \mathbf{U}/\chi \f$ terms during initialisation
  /**
   * Turned on by default
   */
  bool init_compute_uDelu;


//  bool increaseDt;

  /// Distance from the bottom of the domain that a sponge region will extend
  Real spongeHeight;

//  Real postTraceSmoothing;
//  bool skipNewLevelScalars;

  /// Turn on to set the advective source term in the salt equation to 0
  bool skipSaltUpdate;

  /// Turn on to stop evolving the enthalpy and bulk concentration fields
  bool skipHCUpdate;

  /// Use a diffusive source term to compute the advective terms more accurately
  /**
   * When computing things like \f$ \mathbf{U} \cdot \nabla T \f$, we upwind temperature
   * to cell faces. As the temperature also obeys diffuses, we can add a diffusive source term
   * to estimate the upwinded temperatures more accurately. Similarly for other advected and diffused
   * quantities.
   */
  bool doDiffusionSrc;


//  bool noMultigrid;
//  int noMultigridIter;

  /// Fractional increase in the buoyancy force with each timestep
  Real rampBuoyancy;

  /// Initial RaC when we're changing the buoyancy force at each timestep
  Real initRaC;

  /// Initial RaT when we're changing the buoyancy force at each timestep
  Real initRaT;

  /// Max RaC when we're changing the buoyancy force at each timestep
  Real maxRaC;

  /// Max RaT when we're changing the buoyancy force at each timestep
  Real maxRaT;

  /// At which time to compute advection velocities during initialisation
  /**
   * Default is 0.1, i.e. at time 0 + 0.1*dt
   */
  Real initAdvVelCentering;

  /// When computing the advection velocity during multiple initialisation steps,
  /// grow MushyLayerOptions::initAdvVelCentering by this factor with each iteration
  Real adv_vel_centering_growth;

  /// Whether or not we want to ignore solvers failing
  /**
   * If false, we'll try and do some thing about the situation as dictated by
   * MushyLayerOptions::solverFailRestartMethod
   */
  bool ignoreSolveFails;

  /// What to do if a solver fails and we're not ignoring the problem
  int solverFailRestartMethod;


//  bool initiallyDoAutomaticRestart;

  /// When testing for steady state can consider different norms of the solution
  /**
   * 0 -> compute the max value
   * 1 -> compute the L1 norm of the solution (i.e. some of absolute values)
   * 2 -> compute the L2 norm of the solution (i.e. some of values squared)
   */
  int steadyStateNormType;

  /// Whether we want to try and load the advection velocity from a restart file
  /**
   * Turned off by default, as this won't always work with some of the restart files we generate
   * (e.g. after refinement/adding new cells).
   * This is only really useful for when we want to do post processing without recomputing velocity fields.
   */
  bool load_advVel;

  /// Whether to do 1st order (set to 1) or 2nd order (set to 2) interpolation at coarse-fine boundaries
  int CFinterpOrder_advection;

  /// Whether we want to try and do some dodgy stuff to make the enthalpy-bulk concentration nonlinear solver quicker
  /**
   * Turned off by default, use with extreme caution.
   */
  bool nonlinearHCOpSuperOptimised;

  /// Whether to apply BCs to temperature, porosity etc. explicitly during multigrid solves
  /**
   * If false, use BC comptued from enthalpy/bulk concentration
   */
  bool apply_diagnostic_bcs;

  /// Whether or not to do subcycling
  bool useSubcycling;

  /// How much text output to produce on a scale of 0 -> infinity
  int verbosity;

  /// Use slope limiting in advection calculations?
  bool useLimiting;

  /// Number of cells to add around tagged cells before computing new meshes
  int tagBufferSize;

  /// Refinement threshold
  Real refineThresh;

  bool doEulerPart;

  /// Whether or not to compute diagnostics ad hoc
  bool computeDiagnostics;

  // Projection stuff

  /// Whether or not to do projection.
  bool doProjection;

  /// Whether or not to use the previous pressure when solving Darcy's equation
  /**
   * If true, we then only solve for a small extra pressure correction which is usually a lot quicker
   */
  bool useIncrementalPressure;

  /// Whether or not to solve for an incremental pressure change on fine levels, rather than the full pressure
  /**
   * This may not work fully yet...
   */
  bool useIncrementalPressureRefinedLevels;

//  Real phiScale;
//  bool scaleMACBCWithChi;
//  Real MACBCscale;

  /// Whether or not to do synchronisation operations over multiple levels have reached the same time
  bool doSyncOperations;

  /// Whether or not to enforce an analytic solution at the start of the simulation
  bool enforceAnalyticSoln;

  /// Define the analytic solution to apply.
  /**
   * By default, analyticSolution = MushyLayerParams::physicalProblem
   */
  int analyticSolution;

//  bool useAnalyticSource;

  /// Try and make sure the maximum divergence of the face centred velocity is less than this
  Real maxDivUFace;

  /// Whether or not the MAC (face-centred) projection should solve for a pressure that is scaled by the porosity or permeability
  /**
   * If scaled, we solve something like \f$ \nabla \cdot \chi \nabla \phi = \nabla \cdot \mathbf{U}^* \f$,
   * else we'll solve \f$ \nabla^2 \phi = \nabla \cdot \mathbf{U}^* \f$.
   */
  bool scaleP_MAC;

  /// Whether or not the CC (cell-centred) projection should solve for a pressure that is scaled by the porosity
  /**
   * If scaled, we solve something like \f$ \nabla \cdot \chi \nabla \pi = \nabla \cdot \mathbf{U}^* \f$,
   * else we'll solve \f$ \nabla^2 \pi = \nabla \cdot \mathbf{U}^* \f$.
   */
  bool scaleP_CC;

  /// Compute viscous terms explicitly in velocity calculation
  bool explicitViscousVelSolve;


//  bool explicitDarcyTerm;

  /// Use the cell centred pressure \f$ \pi \f$ as the boundary condition when solving for the face centred pressure \f$ \phi \f$
  /**
   * Default is true.
   */
  bool usePiAdvectionBCs;

  /// How much output the projection object should generate
  int projection_verbosity;


//  bool implicitAdvectionSolve;
//  bool usePhiForImplicitAdvectionSolve;

  /// Solve for all velocity componenents at once. Not currently implemented.
  bool multiCompUStarSolve;

  /// Porosity below which we enforce zero velocity.
  /**
   * This needs to be greater than 0
   */
  Real solidPorosity;

  /// Porosity below which we ensure advection velocities are 0
  Real advPorosityLimit;

  /// If > 0, use Darcy's equation to compute velocities in cells with porosity less than this value
  Real chiLimit;

  /// Porosity below which we ensure cell-centred velocities are 0
  Real ccVelPorosityLimit;

  /// Porosity below which we ensure that the src term for computing advection velocities is 0
  Real advVelsrcChiLimit;

  /// Porosity below which we ensure that the \f$ \mathbf{U} \cdot \nabla \mathbf{U} \f$ part of the cell-centred source term is 0
  Real uDelU_porosityLimit;

  /// Porosity below which we ensure that velocities are zero when computing \f$ \mathbf{U} \cdot \nabla \mathbf{U} \f$
  Real advVelChiLimit;

  /// Use advection velocity from previous timestep for tracing advection velocities
  /**
   * Default is false, in which case we average the cell centred velocities from the previous
   * timestep to cell faces and use these (as they are calculated at a later time, they should be
   * more accurate).
   */
  bool useOldAdvVel;

  /// Whether to enforce an analytic velocity field
  bool enforceAnalyticVel;

  /// Whether to project the enforced analytic velocity field
  /**
   * Useful for testing if the projection works as expected
   */
  bool projectAnalyticVel;

  /// If enforcing an analytic velocity, decide what to use
  int analyticVelType;

  /// Whether or not to initialise the velocity to an analytically specified value
  bool initAnalyticVel;

//  int lapVelNumSmooth;
//  Real lapVelSmoothScale;

  /// Maximum number of projections to do on levle 0
  int maxProjBaseLevel;

  /// Maximum number of MAC (face cented) projections to do
  int maxNumMACProj;

  /// Accuracy for boundary conditions applied to laplacian(velocity) in the advection velocity source term
  int lapVelBCOrder;

  /// Where to compute the CC velocity source term (default = 0.5 i.e. halfway through the timestep)
  Real CCVelSrcTermCentering;

  // u.grad(u) stuff
  int uDeluMethod;
  int uDelU_grow;
  bool uDelUConservativeForm;

  // Advection source term

  /// Whether or not to include pressure in the advection velocity source term
  bool advVelPressureSrc;

  /// Whether or not to include the Darcy term in the advection velocity source
  bool advVelDarcySrc;

  /// Whether or not to include the viscous term in the advection velocity source
  bool advVelViscousSrc;

  /// Whether or not to include buoyancy in the advection velocity source term
  bool advVelBuoyancySrc;


//  bool advSrcAllowLaggedLapVel;

  /// Whether or not to include \f$ \mathbf{U} \cdot \nabla \mathbf{U}/\chi \f$ in the cell centred velocity source term
  bool CCAdvSrc;

  /// Whether or not to include buoyancy in the cell centred velocity source term
  bool CCBuoyancySrc;

  /// Whether or not to include pressure in the cell centred velocity source term
  /**
   * Only matters if MushylayerOptions::CCPressureSrcOverride is turned on
   */
  bool CCPressureSrc;

  /// Turn to make sure we use MushyLayerOptions::CCPressureSrc to decide whether or not to include pressure
  /**
   * Otherwise, we will use AMRLevelMushyLayer::m_usePrevPressureForUStar.
   * This is turned on by default if MushyLayerOptions::CCPressureSrc is defined through the inputs file
   */
  bool CCPressureSrcOverride;



  /// Whether or not to reflux momentum
  bool reflux_momentum;

  /// Whether or not to just reflux momentum in the normal direction
  bool reflux_normal_momentum;

  /// Whether or not to reflux the enthalpy field
  bool reflux_enthalpy ;

  /// Whether or not to reflux the bulk concentration field
  bool reflux_concentration ;

  /// Whether or not to reflux lambda
  /**
   * If you don't do this, you won't be able to calculate the freestream preservation correction.
   */
  bool reflux_lambda ;

  /// Whether or not to average down from fine levels to coarser levels after refluxing
  /**
   * As the solution will have changed following refluxing, doing this ensures consistency between levels
   */
  bool refluxAverageDown;

  /// Method to use for implicit refluxing
  RefluxMethod refluxMethod;

  /// Coefficient of the laplacian operator during implicit reflux solves
  /**
   * This might have to change depending on the operator being used, so left as an option.
   */
  Real refluxBetaSign;

  /// Whether to add (+1) or subtract (-1) the correction
  /**
   * This might need to be changed depending on the operator being used, so left as an option.
   */
  Real refluxCorrSign;

  /// Solver tolerance for the viscous solver, which we use for smoothing scalar fields and doing momentum reflux
  Real viscous_solver_tol;

  /// Viscous solver - number of smoothing steps during a down sweep
  int viscous_num_smooth_down;

  /// Viscous solver - number of smoothing steps during an up sweep
  int viscous_num_smooth_up;


//  int AMRMultigridRelaxMode;

  /// Verbosity for the enthalpy-bulk concentration reflux solver
  int AMRMultigridVerb;

  /// Solver tolerance
  Real AMRMultigridTolerance;

  /// Solver hang condition
  Real AMRMultigridHang;

  /// If residual norm is less than this absolute value, solver decides we've converged
  Real AMRMultigridNormThresh;

  /// Multigrid for implicit cell centred \f$ \mathbf{U}^* \f$ solve - num smoothing steps
  int velMGNumSmooth;

  /// Solver tolerance
  Real velMGTolerance;

  /// Solver tolerance
  Real velMGHang;

  /// Number of V-cycles to perform
  int velMGNumMG;

  /// Residual norm condition for convergence
  Real velMGNormThresh;

  /// Max number of multigrid iterations
  int VelMGMaxIter;

  /// Multigrid for the coupled enthalpy-bulk concentration solve: number of smooths on each up sweep
  int HCMultigridNumSmoothUp;

  /// Number of smooths on each down sweep
  int HCMultigridNumSmoothDown;

  /// Number of V-cycles to do
  int HCMultigridNumMG;

  /// Max number of multigrid iterations
  int HCMultigridMaxIter;

  /// How much output to print
  int HCMultigridVerbosity;

  /// Number of iterations of the bottom solver
  int HCMultigridBottomSolveIterations;

  /// Solver tolerance
  Real HCMultigridTolerance;

  /// Hang condition
  Real HCMultigridHang;

  /// Condition on the norm of the residual for convergence
  Real HCMultigridNormThresh;

  /// How to do relaxation
  /**
   * See AMRPoissonOp::relax in Chombo for the various options.
   * In particular, 1 = Gauss-Seidel Red-Black, 4 = Jacobi
   */
  int HCMultigridRelaxMode;

  /// Whether to use a relaxation scheme as the bottom solve
  /**
   * We use this here as it treats the nonlinearity more sensibly than linear solvers like BiCGStab
   */
  bool HCMultigridUseRelaxBottomSolverForHC;


  /// Order of the normal predictor for velocity advection solves
  int velAdvNormalPredOrder;

  /// Whether to use 4th order slopes for velocity advection solves
  bool velAdvUseFourthOrderSlopes;

  /// Whether to use a higher order flux limiter for velocity advection solves
  bool velAdvHigherOrderLimiter;

  /// Whether to use artifical viscosity for velocity advection solves
  bool velAdvUseArtVisc;

  /// Value of artifical viscosity to use for velocity advection solves
  /**
   * Only used if MushylayerOptions::velAdvUseArtVisc is turned on
   */
  Real velAdvArtVisc;

  /// Whether to use artificial viscosity for enthalpy-bulk concentration advection
  bool HCUseArtVisc;

  /// Value of artificial viscosity to use for enthalpy-bulk concentration advection
  Real HCArtVisc;

  /// Order of the normal predictor for enthalpy-bulk concentration advection
  int HCNormalPredOrder;

  /// Whether to use 4th order slopes for enthalpy-bulk concentration advection
  bool HCUseFourthOrderSlopes;

  /// Whether to use a higher order flux limiter for enthalpy-bulk concentration advection
  bool HCHigherOrderLimiter;


  /// Scalar variable to use for working out where to do refinement
  int taggingVar;

  /// Vector variable to use for working out where to do refinement
  int taggingVectorVar;

  /// Strategy for tagging cells. See ::TaggingMethod for the options.
  TaggingMethod taggingMethod;

  /// Don't start regridding until time > min_regrid_time
  Real min_regrid_time;

  /// If this is > 0, then when time > fixed_grid_time we will switch to the grid specified by 'main.regrid_gridfile'
  Real fixed_grid_time;

  /// Whether we should tag the domain boundary for refinement
  bool tagDomainBoundary;

  /// Whether we should tag the mush-liquid boundary for refinement
  bool tagMLboundary;

  /// Method for refinement
  RefinementMethod refinementMethod;

  /// Turn on to only tags cells where other refinement criteria are met and porosity \f$\chi < 1 \f$
  bool onlyTagPorousCells;

  /// Filter porous cells to remove anomalies
  /**
   * Select all cells with \f$\chi < 1\f$, then shift vertically
   * by this value and tag only the intersection of the two
   */
  int porousCellsShrink;

  /// Whether we should tag cells where the fluid velocity (magnitude) is greater than MushyLayerOptions::vel_thresh
  bool tag_velocity;
  
  /// New option for a refinement method based on a composite criteria
  bool compositeChannelTagging;

  /// If MushyLayerOptions::tag_velocity is true, specify the velocity magnitude criteria here
  Real vel_thresh;

  /// Refine mushy regions on level 1, the refine around plumes according to some empirical criteria
  /**
   * Specifically, identify plumes as regions where either
   * a) the downwards velocity is greater than MushyLayerOptions::plumeVelThreshold, and
   * b) the bulk concentration is greater than MushyLayerOptions::plumeSalinityThreshold
   */
  bool tag_plume_mush;

  /// Just refine around channels
  bool tag_channels;

  /// Velocity theshold value for identifying plumes
  Real plumeVelThreshold;

  /// Salinity threshold value for identifying plumes
  Real plumeSalinityThreshold;

  /// If less than 1, use this limit to identify cells with taggingMarginalPorosityLimit < porosity < 1.0
  /**
   * These cells are then not considered part of the mushy layer for tagging purposes, but also not liquid.
   */
  Real taggingMarginalPorosityLimit;

  /// Tag cells in the centre of the domain to make a box of this size. Leave as 0 (default) to do nothing
  int tagCenterBoxSize;

  /// If MushyLayerOptions::tagCenterBoxSize > 0, then regrid once time > tagCenterBoxRegridTime
  Real tagCenterBoxRegridTime;

  /// For testing: if turned on, refined levels will only exist for one regrid cycle before being removed
  /**
   * This let's us test levels being created then destroyed in a controlled manner
   */
  bool testRegridCoarsening;

  /// Whether to do higher order interpolationg for scalars to fill new cells following regridding
  bool scalarHOinterp;

  /// Whether to do higher order interpolation for vectors to fill new cells following regridding
  bool vectorHOinterp;

  /// Whether or not to produce lots of plot files during regridding
  /**
   * Useful for debugging
   */
  bool makeRegridPlots;

  /// Factor by which to scale dt for re-initialising pressure following regridding
  Real regrid_dt_scale;

  /// Whether or not to recalculate lambda after regridding
  bool initLambdaPostRegrid;

//  Real minEta;
//  bool computeFreestreamCorrectionSingleLevel;

  /// Whether or not to advect lambda after regridding, in order to recompute the freestream correction
  bool regrid_advect_before_freestream;

  /// Whether or not to do subcycling when advecting lambda after regridding
  /**
   * Only used if MushyLayerOptions::regrid_advect_before_freestream is true
   */
  bool regrid_freestream_subcycle;

  /// Option to scale \f$ \eta \f$ before recomputing the freestream correction after regridding
  /**
   * Default value = 1.0.
   */
  Real regrid_eta_scale;

  /// Whether to reflux lambda before recomputing freestream correction after regridding
  bool regrid_reflux_lambda;

  /// Whether or not to do some smoothing after regridding
   bool do_postRegrid_smoothing;

  /// Coefficient to determine how aggressively to do smoothing (if MushyLayerOptions::do_postRegrid_smoothing = true)
  Real regrid_smoothing_coeff;

/// Whether or not to project the initial velocity field
  /**
   * Projecting ensures it is divergence free
   */
  bool project_initial_vel;

  /// Whether or not to initialize pressures
  bool initialize_pressures;

  /// Whether or not to initialize pressure after regridding
  bool regrid_init_pressure;

  /// Whether to do bi-linear interpolation when regridding
  bool regrid_linear_interp;

  /// Set the initial value for AMRLevelMushyLayer::m_usePrevPressureForUStar
  /**
   * Determines if we should use the previous \f$ \nabla P \f$ to calculate \f$ \mathbf{u}^* \f$
   */
  bool usePrevPressureForUStar;


//  int iter_plot_interval;

  /// Time period between successive diagnostic reports
  Real diagnostics_period;

  /// Some custom options for initial data
  int customInitData;

  /// Whether to write out plotfiles during pressure initialisation
  /**
   * Useful for debugging
   */
  bool writePressureInitFields;

  /// When initialising pressure via multiple iterations, decide whether or not to reset the other fields between iterations
  bool initResetStates;


//  Real skipTrickySourceTerm;

  /// Advect both enthalpy and bulk concentration at the same time (slightly more efficient)
  bool allowMulticompAdvection;


//  Real smoothingCoeff;

  /// Whether or not to compute the volume discrepancy correction at the start of a simulation
  bool compute_initial_VD_corr;

  /// Order of accuracy for time integration (1 = backward euler, 2 = TGA)
  int timeIntegrationOrder;

  /// How much output multigrid operators should print
//  int AMRMultigrid_verbosity;


  /// We don't let the porosity go below this value
  Real lowerPorosityLimit;

  /// Size of the initial perturbation (default = 0.0)
  Real initialPerturbation;

  /// Whether or not we should do advection and diffusion of scalars (enthalpy, bulk concentration)
  bool doScalarAdvectionDiffusion;

  /// Size of initial velocity if required
  Real initVel;

  /// Phase shift for sinusoidal perturbations
  Real perturbationPhaseShift;

  /// Size of perturbation to add at time MushyLayerOptions::perturbationTime
  Real delayedPerturbation = 0.0;

  /// Specific time at which to introduce a perturbation
  Real perturbationTime = 0.0;

  /// Wavenumber of sinusoidal perturbations (default = 0.0, i.e. just a uniform value)
  Real perturbationWavenumber;

  /// If true, sinusoidal perturbations take the form of a sin wave. Otherwise a cosine wave.
  /**
   * E.g, perturbation \f$ = \alpha \sin (2 \pi N x / L) \f$
   * where \f$\alpha\f$ is MushyLayerOptions::initialPerturbation,
   * \f$ N \f$ is MushyLayerOptions::perturbationWavenumber,
   * and \f$ L \f$ is the domain width
   */
  bool perturbationSin = false;

  /// Add an initial random perturbation to enthalpy
  Real initialRandomPerturbation = 0.0;

  /// Turn on to seed random perturbations on different boxes
  /// to prevent the same 'random' number patterns in each box
  bool seedRandomPert = true;

  /// Porosity value to set if we're using a fixed porosity.
  /**
   * If < 0, use some functional form for the porosity as defined by MushyLayerOptions::porosityFunction
   */
  Real fixedPorosity = -1.0;

  /// Timescale for porosity variations if enforcing a fixed porosity
  Real porosityTimescale;

  /// Functional form to use for the porosity if using a fixed porosity (rather than a porosity computed from the solution).
  PorosityFunctions porosityFunction;

  /// When restarting, if MushyLayerOptions::restartPerturbation > 0.0 then add a perturbation which if the sum of waves with wavenumbers from 1->maxRestartWavenumber
  int maxRestartWavenumbers;

  /// Size of perturbation to add when restarting a simulation
  Real restartPerturbation;


  /// Create a porous hole of this radius during initialisation
  Real porousHoleRadius;

  /// Velocity scale for the initial conditions
  Real initVelScale;


  /// Whether to initialise with a linear enthalpy gradient in the vertical direction
  bool linearGradient;

  /// Approximate depth of mushy layer to create as part of an initial condition approximating some sea ice
  Real mushHeight;

  /// Different options for initialising data that looks a bit like year old sea ice
  /**
   * See AMRLevelMushyLayer::initialDataMushyLayer() for more details
   */
  int summerProfile;

  /// Depth of melt pond to add
  Real meltPondDepth;

  /// Whether or not to horizontally average all fields before restarting
  bool horizAverageRestart;

  /// When restarting, add a perturbation to this variable
  int restart_perturbation_var;

  // Restart

  /// Reference temperature for nondimensionalisation
  Real refTemp;

  /// Reference salinity for nondimensionalisation
  Real refSalinity;

  /// Previous reference salinity for nondimensionalisation (if restarting and different from the new value)
  Real prevRefSalinity;

  /// Previous reference temperature for nondimensionalisation (if restarting and different from the new value)
  Real prevRefTemp;

  /// Facotr by which to reduce dt after restarting: \f$ dt = dt / dtReductionFactor \f$
  Real dtReductionFactor;

  /// Fixed porosity maximum value
  Real fixedPorosityMaxChi;

  /// Standard deviation of the gaussian used to create a fixed porosity field
  Real FixedPorositySTD;

  /// For a fixed porosity which varies linearly from some constant value inside the domain to the edge, this is the radius of the inner region
  Real fixedPorosityFractionalInnerRadius;

  /// If using a time dependent enforced porosity, stop changing it after this time has passed
  Real fixedPorosityEndTime;

  /// Whether to exchange corner cells by default for scalar fields
  bool scalarExchangeCorners;

  /// Set buoyancy forces to zero after this time
  Real buoyancy_zero_time;

  /// If this is set to greater than 0,
//  Real maxEta;

  /// Whether to use linear or FAS geometric multigrid
  /**
   * Should really be FAS as we're solving nonlinear problems, but leaving the option
   * here to change back to linear multigrid.
   */
  MGmethod MGtype;

  /// Specify how to treat the \f$ \mathbf{U} \cdot \nabla \left( \mathbf{U}/\chi \right) \f$ term
  velocityAdvectionTypes advectionMethod;

  /// Whether to include tracer dynamics
  bool includeTracers;

  /// If > 0, compute penetration of irradiance down through the ice
  Real surfaceIrradiance;


};


#endif /* SRC_MUSHYLAYEROPT_H_ */
