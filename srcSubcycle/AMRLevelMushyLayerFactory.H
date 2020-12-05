
#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELMUSHYLAYERFACTORY_H_
#define _AMRLEVELMUSHYLAYERFACTORY_H_

#include "AMRLevelFactory.H"

#include "AMRLevelMushyLayer.H"
#include "mushyLayerOpt.h"

#include "NamespaceHeader.H"

/// Factory to create AMRLevelMushyLayer
/**
 *  Factory to create AMRLevelMushyLayer
 */
class AMRLevelMushyLayerFactory : public AMRLevelFactory
{
public:
  /// Constructor
  /**
   */
  AMRLevelMushyLayerFactory(MushyLayerOptions opt, MushyLayerParams params);

  /// Destructor
  virtual ~AMRLevelMushyLayerFactory() {}

  /// Create new AMRLevel
  virtual AMRLevel *new_amrlevel() const;

protected:
  /// Options for how to solve the equations and run the code
  MushyLayerOptions m_options;

  /// Physical parameters for the system
  MushyLayerParams m_params;

private:
  /// Weak construction is disallowed
  AMRLevelMushyLayerFactory() { MayDay::Error("invalid operator"); }

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelMushyLayerFactory &);
  AMRLevelMushyLayerFactory(const AMRLevelMushyLayerFactory &);
};

#include "NamespaceFooter.H"

#endif
