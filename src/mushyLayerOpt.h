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

/// Different multigrid methods
enum MGmethod {
  /// Normal linear multigrid
  MGTypeStandard,

  /// Full Approximation Scheme (for nonlinear problems)
  MGTypeFAS
};

/// Strategy for tagging cells
enum TaggingMethod {
  /// Tag where the undivided gradient of some field is bigger than some value
  UndividedGradient,

  /// Tag where the magnitude of some field is bigger than some values
  Magnitude,

  /// Tag where the value of some field is larger than some criteria
  CompareLargerThan,

  /// Tag where the value of some field multiplied by -1 is larger than some criteria.
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

/// Identifiers for vector variables
  enum VectorVars {
    m_fluidVel,
    m_U_porosity,
    m_Ustar,
    m_advUstar,
    m_advectionVel,
    m_viscousSolveSrc,
    m_UdelU,
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
    m_enthalpy,
    m_bulkConcentration,
    m_temperature,
    m_porosity,
    m_liquidConcentration,
    m_solidConcentration,
    m_pressure,
    m_permeability,
    m_viscosity,
    m_lambda,
    m_lambda_porosity,
    m_enthalpySolidus,
    m_enthalpyLiquidus,
    m_enthalpyEutectic,
    m_porosityAnalytic,
    m_temperatureAnalytic,
    m_saltEqnSrcGodunov,
    m_saltEqnSrcFiniteDiff,
    m_saltEqnOp,
    m_Terr,
    m_enthalpyOp,
    m_enthalpySrc,
    m_divUadv,
    m_dHdt,
    m_dSdt,
    m_averageVerticalFlux,
    m_soluteFluxAnalytic,
    m_verticalFlux,
    m_saltResidual,
    m_divU,
    m_averageHeatFlux,
    m_streamfunction,
    m_vorticity,
    m_FsVertDiffusion,
    m_FsVertFluid,
    m_FsVertFrame,
    m_divUcorr,
    m_pi,
    m_phi,
    m_MACBC,
    m_MACrhs,
    m_CCrhs,

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

  Real initAdvVelCentering;

  Real adv_vel_centering_growth;
  int solverFailRestartMethod;
  bool ignoreSolveFails;
  bool initiallyDoAutomaticRestart;
  int steadyStateNormType;
  bool load_advVel;
  Real CFinterpOrder_advection;


  bool nonlinearHCOpSuperOptimised;

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

  /// Whether we should tag cells where the fluid velocity (magnitude) is greater than MushyLayerOptions::vel_thresh
  bool tag_velocity;

  /// If MushyLayerOptions::tag_velocity is true, specify the velocity magnitude criteria here
  Real vel_thresh;

  /// Refine mushy regions on level 1, the refine around plumes according to some empirical criteria
  /**
   * Specifically, identify plumes as regions where either
   * a) the downwards velocity is greater than MushyLayerOptions::plumeVelThreshold, and
   * b) the bulk concentration is greater than MushyLayerOptions::plumeSalinityThreshold
   */
  bool tag_plume_mush;

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

  /// Set the initial value for AMRLevelMushyLayer::m_usePrevPressureForUStar
  /**
   * Determines if we should use the previous \f$ \nabla P \f$ to calculate \f$ \mathbf{u}^* \f$
   */
  bool usePrevPressureForUStar;


//  int iter_plot_interval;

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

  bool compute_initial_VD_corr;

  int timeIntegrationOrder;
  int AMRMultigrid_verbosity;

  Real lowerPorosityLimit;
  Real initialPerturbation;

  bool doScalarAdvectionDiffusion;
  Real initVel;

  Real perturbationPhaseShift;
  Real delayedPerturbation = 0.0;
  Real perturbationTime = 0.0;
  Real perturbationWavenumber = 1.0;
  bool perturbationSin = false;
  Real fixedPorosity = -1.0;
  Real porosityTimescale;
  PorosityFunctions porosityFunction;

  int maxRestartWavenumbers;
  Real restartPerturbation;
  Real radius;
  Real initVelScale;

  bool linearGradient;
  Real mushHeight;
  int summerProfile;
  Real meltPondDepth;
  bool horizAverageRestart;
  int restart_perturbation_var;



  // Restart
  Real refTemp;
  Real refSalinity;
  Real prevRefSalinity;
  Real prevRefTemp;
  Real dtReductionFactor;

  // Porosity
  Real fixedPorosityMaxChi;
  Real FixedPorositySTD;
  Real fixedPorosityFractionalInnerRadius;
  Real fixedPorosityEndTime;

  bool scalarExchangeCorners;
  Real buoyancy_zero_time;

  /// If this is set to greater than 0,
  Real maxEta;

  /// Whether to use linear or FAS geometric multigrid
  /**
   * Should really be FAS as we're solving nonlinear problems, but leaving the option
   * here to change back to linear multigrid.
   */
  MGmethod MGtype;

  /// Specify how to treat the \f$ \mathbf{U} \cdot \nabla \left( \mathbf{U}/\chi \right) \f$ term
  velocityAdvectionTypes advectionMethod;

};


#endif /* SRC_MUSHYLAYEROPT_H_ */
