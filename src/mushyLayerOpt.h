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
  /// Domain width
  Real domainWidth;
  Real domainHeight;
  Real cfl;
  Real max_dt_growth;
  Real init_dt_scale;
  Real fixedDt;

  Real initial_dt_multiplier;
  int num_init_passes;
  Real restart_new_time;
  bool init_add_subtract_grad_p;
  bool init_compute_uDelu;
  bool increaseDt;

  Real advVelCentering;
  Real adv_vel_centering_growth;
  int solverFailRestartMethod;
  bool ignoreSolveFails;
  int steadyStateNormType;
  Real CFinterpOrder_advection;

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

  bool doProjection;

  bool doSyncOperations;
  bool enforceAnalyticSoln;

  Real maxDivUFace;
  bool scaleP_MAC;
  bool scaleP_CC;
  bool explicitDarcyTerm;
  bool implicitAdvectionSolve;

  bool multiCompUStarSolve;

  Real solidPorosity;

  bool do_postRegrid_smoothing;
  bool reflux_momentum;
  bool reflux_normal_momentum;

  bool reflux_enthalpy ;
  bool reflux_concentration ;
  bool reflux_lambda ;

  Real viscous_solver_tol;
  int viscous_num_smooth_down;
  int viscous_num_smooth_up;

  Real variable_eta_factor;

  int iter_plot_interval;

  Real skipTrickySourceTerm;
  bool allowMulticompAdvection;

  bool compute_initial_VD_corr;

  int timeIntegrationOrder;
  int verbosity_multigrid;

  Real lowerPorosityLimit;
  Real initialPerturbation;

  bool doScalarAdvectionDiffusion;
  Real initVel;

  Real delayedPerturbation = 0.0;
  Real perturbationTime = 0.0;
  Real perturbationWavenumber = 1.0;
  bool perturbationSin = false;
  Real fixedPorosity = -1.0;
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

  MGmethod MGtype;

  velocityAdvectionTypes advectionMethod;

};


#endif /* SRC_MUSHYLAYEROPT_H_ */
