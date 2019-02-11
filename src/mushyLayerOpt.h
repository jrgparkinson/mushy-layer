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


//TODO: these should probably all be const's
struct mushy_layer_options {
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
};


#endif /* SRC_MUSHYLAYEROPT_H_ */
