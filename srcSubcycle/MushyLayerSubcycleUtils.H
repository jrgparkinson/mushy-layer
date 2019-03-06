
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
#include "memtrack.H"
#include "parstream.H"

#include "NamespaceHeader.H"



/// Make the AMRLevelMushyLayerFactory object
/**
 * Mainly just loading all the MushyLayerOptions.
 */
extern void
getAMRFactory(RefCountedPtr<AMRLevelMushyLayerFactory>&  a_fact);


/// Create the AMR object.
/**
 */
extern void
defineAMR(AMR&                                          a_amr,
          RefCountedPtr<AMRLevelMushyLayerFactory>&  a_fact,
          const ProblemDomain&                          a_prob_domain,
          const Vector<int>&                            a_refRat);

/// Setup the AMR object ready for the run
/**
 * There are three main possibilities:
 * 1) Restart for a restart file
 * 2) Do a new run with fixed grids
 * 3) Do a new run with adaptive grids
 */
extern void
setupAMRForAMRRun(AMR&                                          a_amr);


/****/

#include "NamespaceFooter.H"

#endif
