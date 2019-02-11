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


//TODO: these should probably all be const's
//TODO: add comments to these options so it's obvious what they do
//TODO: collate more options in here rathet than having lots of parmparses in the main code, so it's obvious what all the options are
struct MushyLayerOptions {
  /// Domain width
  Real domainWidth;
  Real cfl;
  Real max_dt_growth;
  Real fixedDt;
  Real initial_dt_multiplier;
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

  bool compute_initial_VD_corr;

  int timeIntegrationOrder;
  int verbosity_multigrid;

  Real lowerPorosityLimit;
  Real initialPerturbation;

  Real delayedPerturbation = 0.0;
  Real perturbationTime = 0.0;
  Real perturbationWavenumber = 1.0;
  bool perturbationSin = false;
  Real fixedPorosity = -1.0;
  PorosityFunctions porosityFunction;

  Real restartPerturbation;

  MGmethod MGtype;

  velocityAdvectionTypes advectionMethod;

};


#endif /* SRC_MUSHYLAYEROPT_H_ */
