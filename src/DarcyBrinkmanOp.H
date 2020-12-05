#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DarcyBrinkmanOp_H_
#define _DarcyBrinkmanOp_H_

#include "AMRPoissonOp.H"
#include "CoefficientInterpolator.H"

#include "NamespaceHeader.H"

/// Operator for solving the Darcy-Brinkman equation
/**
   Operator for solving one component of
   \f[
   (alpha * aCoef(x) * I - beta * (cCoef(x) + Lap) ) phi = rho
   \f]
   over an AMR hierarchy.
 */
class DarcyBrinkmanOp : public AMRPoissonOp
{
public:
  /**
\name DarcyBrinkmanOp functions */
  /*@{*/

  ///
  /**
   */
  DarcyBrinkmanOp()
  {
    m_lambdaNeedsResetting = true;
    m_time = -1;
  }

  ///
  /**
   */
  virtual ~DarcyBrinkmanOp() {}

  ///
  virtual void residualI(LevelData<FArrayBox> &a_lhs,
                         const LevelData<FArrayBox> &a_phi,
                         const LevelData<FArrayBox> &a_rhs,
                         bool a_homogeneous = false);

  ///
  virtual void preCond(LevelData<FArrayBox> &a_correction,
                       const LevelData<FArrayBox> &a_residual);

  ///
  virtual void applyOpI(LevelData<FArrayBox> &a_lhs,
                        const LevelData<FArrayBox> &a_phi,
                        bool a_homogeneous = false);

  virtual void applyOpNoBoundary(LevelData<FArrayBox> &a_lhs,
                                 const LevelData<FArrayBox> &a_phi);

  /*@}*/

  /**
<<<<<<< DarcyBrinkmanOp.H
=======
\name MGLevelOp functions */
  /*@{*/

  /**
>>>>>>> 1.6
calculate restricted residual
a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
   */
  virtual void restrictResidual(LevelData<FArrayBox> &a_resCoarse,
                                LevelData<FArrayBox> &a_phiFine,
                                const LevelData<FArrayBox> &a_rhsFine);

  /*@}*/

  /**
\name DarcyBrinkmanOp functions */
  /*@{*/

  /// For tga stuff
  virtual void setAlphaAndBeta(const Real &a_alpha, const Real &a_beta);

  /// Also calls reset lambda
  virtual void setCoefs(const RefCountedPtr<LevelData<FArrayBox>> &a_aCoef,
                        const RefCountedPtr<LevelData<FluxBox>> &a_bCoef,
                        const Real &a_alpha, const Real &a_beta,
                        const RefCountedPtr<LevelData<FArrayBox>> &a_cCoef);

  /// Should be called before the relaxation parameter is needed.
  virtual void resetLambda();

  /// Compute lambda once alpha, aCoef, beta, bCoef are defined
  virtual void computeLambda();

  ///
  virtual void reflux(const LevelData<FArrayBox> &a_phiFine,
                      const LevelData<FArrayBox> &a_phi,
                      LevelData<FArrayBox> &a_residual,
                      AMRLevelOp<LevelData<FArrayBox>> *a_finerOp);

  /*@}*/

  //! This is called on multigrid operators when their AMR operators
  //! are altered.
  void finerOperatorChanged(const MGLevelOp<LevelData<FArrayBox>> &a_operator,
                            int a_coarseningFactor);

  //! Returns identity coefficient data.
  LevelData<FArrayBox> &identityCoef() { return *m_aCoef; }

  /// For tga stuff
  // This should be turned off for backward euler at the moment
  virtual void diagonalScale(LevelData<FArrayBox> &a_rhs,
                             bool a_kappaWeighted)
  {
    DisjointBoxLayout grids = a_rhs.disjointBoxLayout();
    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      a_rhs[dit()].mult((*m_aCoef)[dit()]);
    }
  }
  /// For tga stuff
  virtual void divideByIdentityCoef(LevelData<FArrayBox> &a_rhs)
  {
    DisjointBoxLayout grids = a_rhs.disjointBoxLayout();
    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      a_rhs[dit()].divide((*m_aCoef)[dit()]);
    }
  }

  //! Sets up a model that modifies b coefficient data when the operator's
  //! time is set.
  //! \param a_bCoefInterpolator A CoefficientInterpolator that will be used
  //!                            to compute the b coefficient at specific
  //!                            times.
  void setBCoefInterpolator(
      RefCountedPtr<
          CoefficientInterpolator<LevelData<FluxBox>, LevelData<FArrayBox>>>
          &a_bCoefInterpolator)
  {
    m_bCoefInterpolator = a_bCoefInterpolator;
  }

  //! Returns the B coefficient.
  LevelData<FluxBox> &BCoef() { return *m_bCoef; }

  /// Allows access to the B coefficient interpolator.
  RefCountedPtr<
      CoefficientInterpolator<LevelData<FluxBox>, LevelData<FArrayBox>>>
  BCoefInterpolator()
  {
    return m_bCoefInterpolator;
  }

  //! Sets the time centering of the operator. This interpolates b coefficient
  //! data at the given time if an interpolator is set.
  void setTime(Real a_time);

  /// Identity operator spatially varying coefficient storage (cell-centered)
  /// --- if you change this call resetLambda()
  RefCountedPtr<LevelData<FArrayBox>> m_aCoef;

  /// Darcy coefficient
  RefCountedPtr<LevelData<FArrayBox>> m_cCoef;

  /// Laplacian operator spatially varying coefficient storage (face-centered)
  /// --- if you change this call resetLambda()
  RefCountedPtr<LevelData<FluxBox>> m_bCoef;

  /// Reciprocal of the diagonal entry of the operator matrix
  LevelData<FArrayBox> m_lambda;

  /// getFlux function which matches interface to AMRPoissonOp
  /** assumes we want to use member-data bCoef, then calls
second getFlux function */
  virtual void getFlux(FluxBox &a_flux, const LevelData<FArrayBox> &a_data,
                       const Box &a_grid, const DataIndex &a_dit,
                       Real a_scale)
  {
    const FluxBox &bCoef = (*m_bCoef)[a_dit];
    const FArrayBox &cCoef = (*m_cCoef)[a_dit];
    getFlux(a_flux, a_data, cCoef, bCoef, a_grid, a_dit, a_scale);
  }

  /// get diffusive flux
  virtual void getFlux(FluxBox &a_flux, const LevelData<FArrayBox> &a_data,
                       const FArrayBox &a_cCoef, const FluxBox &a_bCoef,
                       const Box &a_grid, const DataIndex &a_dit,
                       Real a_scale)
  {
    const FArrayBox &data = a_data[a_dit];
    a_flux.define(a_grid, a_data.nComp());
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      getFlux(a_flux[idir], data, a_cCoef, a_bCoef, a_flux[idir].box(), idir,
              1);
      a_flux[idir] *= a_scale;
    }
  }

protected:
  /// Coarse fine intvect sets
  LayoutData<CFIVS> m_loCFIVS[SpaceDim];

  /// Coarse fine intvect sets
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];

  /// Interpolator for b coefficient data.
  RefCountedPtr<
      CoefficientInterpolator<LevelData<FluxBox>, LevelData<FArrayBox>>>
      m_bCoefInterpolator;

  /// Current time.
  Real m_time;

  /// Does the relaxation coefficient need to be reset?
  bool m_lambdaNeedsResetting;

  /// Gauss-Seidel relaxation
  virtual void levelGSRB(LevelData<FArrayBox> &a_phi,
                         const LevelData<FArrayBox> &a_rhs);

  /// Multicolor relaxation
  virtual void levelMultiColor(LevelData<FArrayBox> &a_phi,
                               const LevelData<FArrayBox> &a_rhs);

  /// Gauss-Seidel relaxation
  virtual void looseGSRB(LevelData<FArrayBox> &a_phi,
                         const LevelData<FArrayBox> &a_rhs);

  /// Gauss-Seidel relaxation
  virtual void overlapGSRB(LevelData<FArrayBox> &a_phi,
                           const LevelData<FArrayBox> &a_rhs);

  /// Gauss-Seidel relaxation
  virtual void levelGSRBLazy(LevelData<FArrayBox> &a_phi,
                             const LevelData<FArrayBox> &a_rhs);

  /// jacobi relaxation
  virtual void levelJacobi(LevelData<FArrayBox> &a_phi,
                           const LevelData<FArrayBox> &a_rhs);

  /// computes flux over face-centered a_facebox.
  virtual void getFlux(FArrayBox &a_flux, const FArrayBox &a_data,
                       const FArrayBox &a_cCoef, const FluxBox &a_bCoef,
                       const Box &a_facebox, int a_dir, int a_ref = 1) const;

#if 0
	/// utility function which computes operator after all bc's have been set
	void computeOperatorNoBCs(LevelData<FArrayBox>&       a_lhs,
			const LevelData<FArrayBox>& a_phi);
#endif
};

/// Factory to create DarcyBrinkmanOp
/**
   Factory to create DarcyBrinkmanOps
 */
class DarcyBrinkmanOpFactory : public AMRLevelOpFactory<LevelData<FArrayBox>>
{
public:
  /// Default constructor
  DarcyBrinkmanOpFactory();

  virtual ~DarcyBrinkmanOpFactory() {}

  /// define method
  /**
a_coarseDomain is the domain at the coarsest level.
a_grids is the AMR  hierarchy.
a_refRatios are the refinement ratios between levels.  The ratio lives
   with the coarser level so a_refRatios[ilev] is the ratio between
   ilev and ilev+1
a_coarseDx is the grid spacing at the coarsest level.
a_bc holds the boundary conditions.
a_alpha is the identity constant coefficient
a_beta is the laplacian constant coefficient.
a_aCoef is the identity spatially varying coefficient
a_bCoef is the laplacian spatially varying coefficient.
   */
  void define(const ProblemDomain &a_coarseDomain,
              const Vector<DisjointBoxLayout> &a_grids,
              const Vector<int> &a_refRatios, const Real &a_coarsedx,
              BCHolder a_bc, const Real &a_alpha,
              Vector<RefCountedPtr<LevelData<FArrayBox>>> &a_aCoef,
              const Real &a_beta,
              Vector<RefCountedPtr<LevelData<FluxBox>>> &a_bCoef,
              Vector<RefCountedPtr<LevelData<FArrayBox>>> &a_cCoef);

  //! Defines a factory for DarcyBrinkmanOp which allows the operators to
  //! allocate their own coefficient data. \f$\alpha\f$ and \f$\beta\f$
  //! coefficients are initialized to 1.
  //! \param a_coarseDomain The domain at the coarsest level.
  //! \param a_grids The disjoint box layouts for the various refinement levels.
  //! \param a_refRatios The refinement ratios between levels.
  //! \param a_coarsedx The grid spacing at the coarsest level.
  //! \param a_bc The boundary condition imposed on the solution.
  //! \param a_ghostVect The ghost stencil to use in the created coefficient
  //! data.
  void define(const ProblemDomain &a_coarseDomain,
              const Vector<DisjointBoxLayout> &a_grids,
              const Vector<int> &a_refRatios, const Real &a_coarsedx,
              BCHolder a_bc, const IntVect &a_ghostVect);

  /// make new MGLevelOp
  virtual MGLevelOp<LevelData<FArrayBox>> *
  MGnewOp(const ProblemDomain &a_FineindexSpace, int a_depth,
          bool a_homoOnly = true);

  /// make new AMRLevelOp
  virtual AMRLevelOp<LevelData<FArrayBox>> *
  AMRnewOp(const ProblemDomain &a_indexSpace);

  /// ref ratio to next finer level
  virtual int refToFiner(const ProblemDomain &a_domain) const;

  /// coefficient averaging method
  int m_coefficient_average_type;

private:
  void setDefaultValues();

  /// Problem domains for AMR hierarchy
  Vector<ProblemDomain> m_domains;

  /// DisjointBoxLayouts for AMR hierarchy
  Vector<DisjointBoxLayout> m_boxes;

  /// Cell spacings for AMR hierarchy
  Vector<Real> m_dx;

  /// Refinement to next coarser level
  Vector<int> m_refRatios;

  /// Physical boundary conditions
  BCHolder m_bc;

  /// Coefficient \f$ \alpha \f$
  Real m_alpha;

  /// Coefficient \f$ \beta \f$
  Real m_beta;

  /// Coefficient \f$ a(\mathbf{x}) \f$
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_aCoef;

  /// Coefficient \f$ c(\mathbf{x}) \f$
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_cCoef;

  /// Coefficient \f$ b(\mathbf{x}) \f$
  Vector<RefCountedPtr<LevelData<FluxBox>>> m_bCoef;

  /// Lambda, \f$ \lambda \f$, not currently used?
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_lambda;

  /// Copiers for exchange calls
  Vector<Copier> m_exchangeCopiers;

  /// Coarse fine regions
  Vector<CFRegion> m_cfregion;
};

#include "NamespaceFooter.H"
#endif
