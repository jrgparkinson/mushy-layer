#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PHASEDIAGRAM_H_
#define _PHASEDIAGRAM_H_

#include <cmath>

#include "AMR.H"
#include "memtrack.H"
#include "parstream.H"
#include "ParmParse.H"
#include "MushyLayerParams.h"

#include "NamespaceHeader.H"

void updateEnthalpyVariables(LevelData<FArrayBox>& HC,
                             LevelData<FArrayBox>& temperature, LevelData<FArrayBox>& compositionLiquid,
                             LevelData<FArrayBox>& compositionSolid, LevelData<FArrayBox>& porosity,
                             LevelData<FArrayBox>& enthalpySolid, LevelData<FArrayBox>& enthalpyLiquid,
                             LevelData<FArrayBox>& enthalpyEutectic,
                             MushyLayerParams a_params);

/// Update thermodynamic variables \f$ \theta, \chi, \Theta_l, \Theta_s\f$ from \f$H, \Theta\f$ and the phase diagram
void updateEnthalpyVariables(LevelData<FArrayBox>& HC,
                             LevelData<FArrayBox>& theta, LevelData<FArrayBox>& compositionLiquid,
                             LevelData<FArrayBox>& compositionSolid, LevelData<FArrayBox>& porosity,
                             MushyLayerParams a_params);

/// Update thermodynamic variables \f$ \theta, \chi, \Theta_l, \Theta_s\f$ from \f$H, \Theta\f$ and the phase diagram
void updateEnthalpyVariables(LevelData<FArrayBox>& enthalpy, LevelData<FArrayBox>& composition,
                LevelData<FArrayBox>& theta, LevelData<FArrayBox>& compositionLiquid,
                LevelData<FArrayBox>& compositionSolid, LevelData<FArrayBox>& porosity,
                MushyLayerParams m_parameters);


/// Update thermodynamic variables \f$ \theta, \chi, \Theta_l, \Theta_s\f$ from \f$H, \Theta\f$ and the phase diagram
void updateEnthalpyVariables(LevelData<FArrayBox>& enthalpy, LevelData<FArrayBox>& composition,
                LevelData<FArrayBox>& theta, LevelData<FArrayBox>& compositionLiquid,
                LevelData<FArrayBox>& compositionSolid, LevelData<FArrayBox>& porosity,
                LevelData<FArrayBox>& enthalpySolid, LevelData<FArrayBox>& enthalpyLiquid,
                LevelData<FArrayBox>& enthalpyEutectic,
                MushyLayerParams m_parameters );


//Real computePorosity(Real H, Real a_C, Real H_s, Real H_l, Real H_e, MushyLayerParams a_params);
//Real computeCl(Real H, Real a_C, Real H_s, Real H_l, Real H_e, MushyLayerParams a_params);
//Real computeCs(Real H, Real a_C, Real H_s, Real H_l, Real H_e, MushyLayerParams a_params);
//Real computeTemperature(Real H, Real a_C, Real H_s, Real H_l, Real H_e, MushyLayerParams a_params);

/// Compute the porosity \f$ \chi \f$ given enthalpy, bulk concentration, and physical parameters
Real computePorosity(Real H, Real C, Real compositionRatio, Real specificHeatRatio,
                               Real stefan, Real waterDistributionCoeff, Real heatCapacityRatio,
                               Real thetaEutectic, Real ThetaEutectic);

/// Compute the porosity \f$ \chi\f$ within a mushy layer
Real computePorosityMushyLayer(Real H, Real a_C, Real compositionRatio, Real specificHeatRatio,
                               Real stefan, Real waterDistributionCoeff);

/// Compute the transitions in the phase diagram
void computeBoundingEnergy(const Real H, const Real C, Real& H_s, Real& H_l, Real& H_e,
                const Real heatCapacityRatio, const Real stefan, const Real compositionRatio,
                const Real waterDistributionCoeff,
                const Real thetaEutectic, const Real ThetaEutectic);

/// Compute the transitions in the phase diagram
void computeBoundingEnergy(const LevelData<FArrayBox>& a_H, const LevelData<FArrayBox>& a_C,
                LevelData<FArrayBox>& a_Hs, LevelData<FArrayBox>& a_He, LevelData<FArrayBox>& a_Hl,
                MushyLayerParams a_params);

/// Compute thermodynamic variables \f$ \theta, \chi, \Theta_l, \Theta_s\f$ from \f$H, \Theta\f$ and the phase diagram
void computeEnthalpyVars(const Real H, const Real Theta, Real&porosity, Real& theta, Real& ThetaL, Real& ThetaS,
                const Real H_s, const Real H_l, const Real H_e,
                const Real specificHeatRatio, const Real Stefan,
                const Real compositionRatio, const Real waterDistributionCoeff,
                const Real thetaEutectic, const Real ThetaEutectic);
#include "NamespaceFooter.H"

#endif
