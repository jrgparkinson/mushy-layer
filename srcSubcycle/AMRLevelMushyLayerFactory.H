
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

#include "AdvectionFunctions.H"
#include "AMRLevelMushyLayer.H"

#include "NamespaceHeader.H"

/// Factory to create AMRLevelMushyLayer
/**
 *  Factory to create AMRLevelMushyLayer
 */
class AMRLevelMushyLayerFactory: public AMRLevelFactory
{
public:

  /// Constructor
  /**
   */
  AMRLevelMushyLayerFactory(
                               /// CFL number
                               const Real&                 a_cfl,
                               /// physical width of domain
                               const Real&                 a_domainWidth,
                               /// undivided gradient size over which a cell will be tagged for refinement
                               const Real&                 a_refineThresh,
                               /// number of buffer cells around each tagged cell that will also be tagged
                               const int&                  a_tagBufferSize,
                               /// whether to use van Leer limiting
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
                                                            const Real initial_dt_multiplier);


  /// Destructor
  virtual ~AMRLevelMushyLayerFactory()
  {
  }


  /// Create new AMRLevel
  virtual AMRLevel* new_amrlevel() const;

protected:

  /// CFL number
  Real                  m_cfl;

  /// Domain width
  Real                  m_domainWidth;

  /// Refinement threshold
  Real                  m_refineThresh;

  /// Tag buffer size
  int                   m_tagBufferSize;

  /// Use slope limiting in advection calculations?
  bool                  m_useLimiting;

  int m_CFInterpOrder;
  int m_steadyStateNormType;
  Real m_fixedDt;
  Real m_max_dt_growth;
  int  m_verbosity;
  bool m_useSubcycling;

   bool m_ignoreSolverFails;
   int m_solverFailRestartMethod;
   Real m_adv_vel_centering_growth;
   Real m_initial_dt_multiplier;



private:

  /// Weak construction is disallowed
  AMRLevelMushyLayerFactory()
  {
    m_cfl = 0;
    m_domainWidth = 0;
    m_refineThresh = 0;
    m_tagBufferSize = 0;
    m_useLimiting = 0;
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelMushyLayerFactory&);
  AMRLevelMushyLayerFactory(const AMRLevelMushyLayerFactory&);
};

#include "NamespaceFooter.H"

#endif
