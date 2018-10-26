
#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTDIFFUSEUTILS_H_
#define _ADVECTDIFFUSEUTILS_H_

#include <cmath>

#include "AMR.H"
#include "AMRLevelMushyLayerFactory.H"
#include "AdvectionFunctions.H"
#include "memtrack.H"
#include "parstream.H"

#include "NamespaceHeader.H"

///
/**
 */


///
/**
 */
extern void
getAMRFactory(RefCountedPtr<AMRLevelMushyLayerFactory>&  a_fact);


///
/**
 */

///
/**
 */
extern void
defineAMR(AMR&                                          a_amr,
          RefCountedPtr<AMRLevelMushyLayerFactory>&  a_fact,
          const ProblemDomain&                          a_prob_domain,
          const Vector<int>&                            a_refRat);

///
/**
 */
extern void
setupAMRForAMRRun(AMR&                                          a_amr);


/****/

#include "NamespaceFooter.H"

#endif
