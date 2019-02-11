
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
AMRLevelMushyLayerFactory( MushyLayerOptions a_opt)
{
  m_options = a_opt;
}

// Virtual constructor
AMRLevel* AMRLevelMushyLayerFactory::new_amrlevel() const
{
  // Create a new AMRLevelAdvectDiffuse
  AMRLevelMushyLayer* amrMLPtr = new AMRLevelMushyLayer(m_options);

  // Return it
  return (static_cast <AMRLevel*> (amrMLPtr));
}

#include "NamespaceFooter.H"
