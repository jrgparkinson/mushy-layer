
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
AMRLevelMushyLayerFactory( MushyLayerOptions a_opt, MushyLayerParams a_params) : m_options(a_opt), m_params(a_params) {}

// Virtual constructor
AMRLevel* AMRLevelMushyLayerFactory::new_amrlevel() const
{
  // Create a new AMRLevelAdvectDiffuse
  AMRLevelMushyLayer* amrMLPtr = new AMRLevelMushyLayer(m_options, m_params);

  // Return it
  return (static_cast <AMRLevel*> (amrMLPtr));
}

#include "NamespaceFooter.H"
