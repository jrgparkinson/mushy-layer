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

enum RefluxMethod {
  /// 0 - linear correction
  LinearReflux,

  /// 1 - linear VC correction
  LinearVCReflux,

  /// 2 - nonlinear correction
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

  int maxProjBaseLevel;
  int maxNumMACProj;

  int lapVelBCOrder;
  Real CCVelSrcTermCentering;
  bool legacyComputePredictVel;

  // u.grad(u) stuff
  int uDeluMethod;
  int uDelU_grow;
  bool uDelUConservativeForm;

  // Advection source term
  bool advVelPressureSrc;
  bool advVelDarcySrc;
  bool advVelViscousSrc;
  bool advVelBuoyancySrc;
  bool advSrcAllowLaggedLapVel;

  bool CCAdvSrc;
//  bool CCDarcySrc;
  bool CCBuoyancySrc;
  bool CCPressureSrc;
  bool CCPressureSrcOverride;

  bool do_postRegrid_smoothing;
  bool reflux_momentum;
  bool reflux_normal_momentum;

  bool reflux_enthalpy ;
  bool reflux_concentration ;
  bool reflux_lambda ;

  bool refluxAverageDown;
  RefluxMethod refluxMethod;
  Real refluxBetaSign;
  Real refluxCorrSign;

  Real viscous_solver_tol;
  int viscous_num_smooth_down;
  int viscous_num_smooth_up;

  int AMRMultigridRelaxMode;
  int AMRMultigridVerb;
  Real AMRMultigridTolerance;
  Real AMRMultigridHang;
  Real AMRMultigridNormThresh;

  int velMGNumSmooth;
  Real velMGTolerance;
  Real velMGHang;
  int velMGNumMG;
  Real velMGNormThresh;
  int VelMGMaxIter;

  int HCMultigridNumSmoothUp;
  int HCMultigridNumSmoothDown;
  int HCMultigridNumMG;
  int HCMultigridMaxIter;
  int HCMultigridVerbosity;
  int HCMultigridBottomSolveIterations;
  Real HCMultigridTolerance;
  Real HCMultigridHang;
  Real HCMultigridNormThresh;
  int HCMultigridRelaxMode; // 1=GSRB, 4=jacobi
  bool HCMultigridUseRelaxBottomSolverForHC;

  int velAdvNormalPredOrder;
  bool velAdvUseFourthOrderSlopes;
  bool velAdvHigherOrderLimiter;
  bool velAdvUseArtVisc;
  Real velAdvArtVisc;
  bool HCUseArtVisc;
  Real HCArtVisc;
  int HCNormalPredOrder;
  bool HCUseFourthOrderSlopes;
  bool HCHigherOrderLimiter;

  int taggingVar;
  int taggingVectorVar;
  TaggingMethod taggingMethod;
  Real min_regrid_time;
  Real fixed_grid_time;
  bool tagDomainBoundary;
  bool tagMLboundary;
  bool tag_velocity;
  Real vel_thresh;
  bool tag_plume_mush;
  Real plumeVelThreshold;
  Real plumeSalinityThreshold;
  Real taggingMarginalPorosityLimit;
  Real regridTime;
  // Tag cells in the centre of the domain to make a box of this size. Leave as 0 (default) to do nothing
  int tagCenterBoxSize;
  bool testRegridCoarsening;

  bool scalarHOinterp;
  bool vectorHOinterp;

  bool makeRegridPlots;
  Real regrid_dt_scale;

  bool initLambda;
//  Real variable_eta_factor;
  Real minEta;
  bool computeFreestreamCorrectionSingleLevel;
  bool regrid_advect_before_freestream;
  bool regrid_freestream_subcycle;
  Real regrid_eta_scale;

  /// default = 0: no smoothing
  Real regrid_smoothing_coeff;

  bool project_initial_vel;
  bool initialize_pressures;
  bool addSubtractGradP;

  int iter_plot_interval;

  int customInitData;
  bool writePressureInitFields;
  bool initResetStates;

  Real skipTrickySourceTerm;
  bool allowMulticompAdvection;
  Real smoothingCoeff;

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
