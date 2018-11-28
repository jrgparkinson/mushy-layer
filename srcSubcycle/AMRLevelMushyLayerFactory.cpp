
#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRLevel.H"

#include "AMRLevelMushyLayerFactory.H"

#include "NamespaceHeader.H"

//todo - mushy layer factory should own far more properties than this
// for example can probably contain mushy layer params in this class

AMRLevelMushyLayerFactory::
AMRLevelMushyLayerFactory(
                             const Real&                 a_cfl,
                             const Real&                 a_domainLength,
                             const Real&                 a_refineThresh,
                             const int&                  a_tagBufferSize,
                             const bool&                 a_useLimiting,
                             const int a_CFInterpOrder,
                             const int a_steadyStateNormType,
                             const Real a_fixedDt,
                             const Real a_max_dt_growth,
                             const int  a_verbosity,
                             const bool a_useSubcycling,
                             const bool ignoreSolveFails,
                             const int solverFailRestartMethod,
                             const Real adv_vel_centering_growth,
                             const Real initial_dt_multiplier)
{
  m_cfl                    = a_cfl;
  m_domainWidth           = a_domainLength;
  m_refineThresh           = a_refineThresh;
  m_tagBufferSize          = a_tagBufferSize;
  m_useLimiting            = a_useLimiting;

  m_CFInterpOrder = a_CFInterpOrder;
  m_steadyStateNormType =  a_steadyStateNormType;
  m_fixedDt =  a_fixedDt;
  m_max_dt_growth=          a_max_dt_growth;
  m_verbosity=       a_verbosity;
  m_useSubcycling=        a_useSubcycling;

  m_ignoreSolverFails = ignoreSolveFails;
  m_solverFailRestartMethod = solverFailRestartMethod;
  m_adv_vel_centering_growth = adv_vel_centering_growth;
  m_initial_dt_multiplier = initial_dt_multiplier;


}

// Virtual constructor
AMRLevel* AMRLevelMushyLayerFactory::new_amrlevel() const
{
  // Create a new AMRLevelAdvectDiffuse
  AMRLevelMushyLayer* amrMLPtr = new AMRLevelMushyLayer(
      m_cfl,
      m_domainWidth,
      m_refineThresh,
      m_tagBufferSize,
      m_useLimiting,
      m_CFInterpOrder,
      m_steadyStateNormType,
      m_fixedDt,
      m_max_dt_growth,
      m_verbosity,
      m_useSubcycling,
      m_ignoreSolverFails, m_initial_dt_multiplier, m_adv_vel_centering_growth, m_solverFailRestartMethod);

  // Return it
  return (static_cast <AMRLevel*> (amrMLPtr));
}

#include "NamespaceFooter.H"
