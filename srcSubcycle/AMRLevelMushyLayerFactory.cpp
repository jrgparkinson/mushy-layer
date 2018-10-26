
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

AMRLevelMushyLayerFactory::
AMRLevelMushyLayerFactory(
                             const Real&                 a_cfl,
                             const Real&                 a_domainLength,
                             const Real&                 a_refineThresh,
                             const int&                  a_tagBufferSize,
                             const bool&                 a_useLimiting)
{
  m_cfl                    = a_cfl;
  m_domainWidth           = a_domainLength;
  m_refineThresh           = a_refineThresh;
  m_tagBufferSize          = a_tagBufferSize;
  m_useLimiting            = a_useLimiting;

}

// Virtual constructor
AMRLevel* AMRLevelMushyLayerFactory::new_amrlevel() const
{
  // Create a new AMRLevelAdvectDiffuse
  AMRLevelMushyLayer* amrGodPtr = new AMRLevelMushyLayer(
                                                               m_cfl,
                                                               m_domainWidth,
                                                               m_refineThresh,
                                                               m_tagBufferSize,
                                                               m_useLimiting);

  // Return it
  return (static_cast <AMRLevel*> (amrGodPtr));
}

#include "NamespaceFooter.H"
