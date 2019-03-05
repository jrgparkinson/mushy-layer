#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef AMRTOOLS_LEVELDOMAINFLUXREGISTER_H_
#define AMRTOOLS_LEVELDOMAINFLUXREGISTER_H_

#include "REAL.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "IntVectSet.H"
#include "LoHiSide.H"
#include "LevelData.H"
#include "LayoutData.H"
#include "ProblemDomain.H"

#include "NamespaceHeader.H"

/// Keep track of fluxes at domain edges on a level
/*
 * Author: Jamie Parkinson.
 */
class LevelDomainFluxRegister
{
public:
  LevelDomainFluxRegister ();
  virtual
  ~LevelDomainFluxRegister ();

  /// Full define function
  void define(const ProblemDomain a_domain,
              const DisjointBoxLayout a_grids,
              const int a_refRat,
              const Real a_dx,
              LevelDomainFluxRegister* a_fineFR,
              LevelDomainFluxRegister* a_coarseFR,
              int a_numComp = 1);

  /// Increment flux registers with flux
  void incrFlux(const LevelData<FluxBox>& a_flux,
                const Real a_scale,
                const int a_comp = 0,
                const int a_localComp = 0,
                const int a_numComps = 1);

  /// Get the flux on this level
  /**
   * Somewhat meaningless function given how all levels above the base
   * are not guaranteed to cover the entire domain edge.
   * Just included for completeness.
   */
//  Real getFluxLevel(const int a_dir,
//                 const Side::LoHiSide a_side,
//                 const Real a_scale);

  /// Set the flux registers to zero
  void setToZero();

  /// Collate fluxes over entire AMR hierarchy
  /**
   * Returns the flux through the specify domain boundary, by calculating
   * the AMR sum of the fluxes over all levels.
   * Probably more useful than getFluxLevel in most cases.
   * In particular, useful for verifying that the numerical scheme is
   * globally flux conservative.
   */
  Real getFluxHierarchy(const int a_dir,
               const Side::LoHiSide a_side,
               const Real a_scale,
               const int a_localComp = 0);
protected:

  /// Problem domain on this level
  ProblemDomain m_probDomain;

  /// Fluxes at the bottom of the domain in each direction
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_fluxHi,

  /// Fluxes at the top of the domain in each direction
  m_fluxLo;

  /// Grids on this level
  DisjointBoxLayout m_grids;

  /// Coarser  domain flux registers (if they exist)
  LevelDomainFluxRegister *m_fineFR,

  /// Finer  domain flux registers (if they exist)
  *m_coarseFR;

  /// Whether object is defined
  bool m_defined;

  /// Refinement ratio to finer level
  int m_refRat;

  /// Grid spacing on this level
  Real m_dx;

  /// Number of components to keep track of
  int m_numComp;


};


#include "NamespaceFooter.H"

#endif /* AMRTOOLS_LEVELDOMAINFLUXREGISTER_H_ */
