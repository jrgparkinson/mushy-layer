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
/// Options for handling advection velocities
enum velocityAdvectionTypes
{
  m_porosityInAdvection,
  m_porosityOutsideAdvection,
  m_noPorosity
};


enum PorosityFunctions
{
  constantLiquid,
  linear,
  gaussian,
  gaussianSinusoidal,
  constantSmall,
  cubic,
  hyperbolicTan,
};

/// Identifiers for different scalar variables
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

    // Make sure this comes last!
    m_numVectorVars
  };

  enum ScalarVars {
    m_enthalpy,
    m_bulkConcentration,
    m_temperature,
    m_porosity,
    m_liquidConcentration,
    m_solidConcentration,
    m_pressure,
    m_permeability,
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


    // Make sure this comes last!
    m_numScalarVars

  };


//TODO: these should probably all be const's
//TODO: add comments to these options so it's obvious what they do
//TODO: collate more options in here rather than having lots of parmparses in the main code, so it's obvious what all the options are
struct MushyLayerOptions {

  string output_dir;
  string plotfile_prefix;

  bool minimalOutput;
  bool debug;

  Real steadyStateCondition;
  bool ignoreVelocitySteadyState;
  bool ignoreBulkConcSteadyState;

  /// Domain width
  Real domainWidth;
  Real domainHeight;
  Real cfl;
  bool forceUseUChiForCFL;
  Real min_time;
  Real max_dt_growth;
  Real init_dt_scale;
  Real fixedDt;

  Real initial_dt_multiplier;
  Real minDt;
  Real max_dt;
  Real max_init_dt;
  Real accelCFL;
  bool printAccelDt;
  bool useAccelDt;
  bool useInitAccelDt;

  int max_possible_level;

  bool computeVorticityStreamFunction;

  bool useFortranRegularisationFace;
  bool useFortranRegularisation;

  Real stokesDarcyForcingTimescale;

  int num_init_passes;
  Real restart_new_time;
  bool init_add_subtract_grad_p;
  bool init_compute_uDelu;
  bool increaseDt;

  Real spongeHeight;
  Real postTraceSmoothing;
  bool skipNewLevelScalars;
  bool skipSaltUpdate;
  bool skipHCUpdate;
  bool doDiffusionSrc;

  bool noMultigrid;
  int noMultigridIter;

  Real rampBuoyancy;
  Real initRaC;
  Real initRaT;
  Real maxRaC;
  Real maxRaT;

  Real advVelCentering;
  Real adv_vel_centering_growth;
  int solverFailRestartMethod;
  bool ignoreSolveFails;
  int steadyStateNormType;
  Real CFinterpOrder_advection;


  bool nonlinearHCOpSuperOptimised;

  /// Whether or not to do subcycling
  bool useSubcycling;
  int verbosity;
  /// Use slope limiting in advection calculations?
  bool useLimiting;

  /// Tag buffer size
  int tagBufferSize;
  /// Refinement threshold
  Real refineThresh;

  bool doEulerPart;

  bool computeDiagnostics;

  // Projection stuff
  bool doProjection;
  bool useIncrementalPressure;
  Real phiScale;
  bool scaleMACBCWithChi;
  Real MACBCscale;

  bool doSyncOperations;
  bool enforceAnalyticSoln;
  bool useAnalyticSource;

  Real maxDivUFace;
  bool scaleP_MAC;
  bool scaleP_CC;
  bool explicitDarcyTerm;
  bool usePiAdvectionBCs;
  int projection_verbosity;

  bool implicitAdvectionSolve;
  bool usePhiForImplicitAdvectionSolve;

  bool multiCompUStarSolve;

  Real solidPorosity;
  Real advPorosityLimit;
  Real chiLimit;
  Real ccVelPorosityLimit;
  Real advVelsrcChiLimit;
  Real uDelU_porosityLimit;
  Real advVelChiLimit;

  bool useOldAdvVel;
  bool enforceAnalyticVel;
  bool projectAnalyticVel;
  int analyticVelType;

  int lapVelNumSmooth;
  Real lapVelSmoothScale;

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
  bool CCDarcySrc;
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
  Real variable_eta_factor;
  Real minEta;
  bool computeFreestreamCorrection;
  bool regrid_advect_before_freestream;
  bool regrid_freestream_subcycle;
  Real regrid_eta_scale;

  int iter_plot_interval;

  Real skipTrickySourceTerm;
  bool allowMulticompAdvection;
  Real smoothingCoeff;

  bool compute_initial_VD_corr;

  int timeIntegrationOrder;
  int verbosity_multigrid;

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

  Real maxEta;

  MGmethod MGtype;

  velocityAdvectionTypes advectionMethod;

};


#endif /* SRC_MUSHYLAYEROPT_H_ */
