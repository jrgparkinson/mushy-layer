#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRNonLinearMultiCompOp_H_
#define _AMRNonLinearMultiCompOp_H_

#include "AMRPoissonOp.H"

#include "CoefficientInterpolatorLinear.H"
#include "CoefficientInterpolatorLinearFace.H"
#include "MushyLayerParams.h"
//#include "EnthalpyVariable.H"
#include "PhysBCUtil.H"
#include "VelBCHolder.H"
//#include "ComputeEnthalpyVars.h"


#include "NamespaceHeader.H"

/// Nonlinear variable coefficient operator
/**
   Operator for solving variable-coefficient
   (alpha * aCoef(x) * I - beta * Div(bCoef(x) . Grad)) phi = rho
   over an AMR hierarchy.
 */
class AMRNonLinearMultiCompOp : public AMRPoissonOp
{
public:
  /**
     \name AMRNonLinearMultiCompOp functions */
  /*@{*/
  AMRNonLinearMultiCompOp()
{
    m_lambdaNeedsResetting = true;
    m_time = -1;
    m_FAS = false;
    m_params = nullptr;
    m_apply_bcs_to_diagnostic_var = true;
    m_superOptimised = false;
}

  ///
  /**
   */
  virtual ~AMRNonLinearMultiCompOp()
  {
  }

  /// Compute residual
  virtual void residualI(LevelData<FArrayBox>&       a_lhs,
                         const LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>& a_rhs,
                         bool                        a_homogeneous = false);

  /// Preconditioner
  virtual void preCond(LevelData<FArrayBox>&       a_correction,
                       const LevelData<FArrayBox>& a_residual);

  /// Apply operator
  virtual void applyOpI(LevelData<FArrayBox>&       a_lhs,
                        const LevelData<FArrayBox>& a_phi,
                        bool                        a_homogeneous = false);

  /// Apply operator
  virtual void applyOpNoBoundary(LevelData<FArrayBox>&       a_lhs,
                                 const LevelData<FArrayBox>& a_phi);

  /// Apply operator
  void applyOpNoBoundary2(LevelData<FArrayBox>&      a_lhs,
                          const LevelData<FArrayBox>& a_phi,
                          bool                        a_homogeneous );

  /// Scale with identity coefficient
  virtual void diagonalScale(LevelData<FArrayBox>& a_rhs, bool a_kappaWeighted)
  {
    DisjointBoxLayout grids = a_rhs.disjointBoxLayout();
    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      a_rhs[dit()].mult((*m_aCoef)[dit()]);
    }
  }

  /// Divide by identity coefficient
  virtual void divideByIdentityCoef(LevelData<FArrayBox>& a_rhs)
  {
    DisjointBoxLayout grids = a_rhs.disjointBoxLayout();
    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      a_rhs[dit()].divide((*m_aCoef)[dit()]);
    }
  }

  /*@}*/

  /**
<<<<<<< AMRNonLinearMultiCompOp.H
=======
     \name MGLevelOp functions */
  /*@{*/

  /// calculate restricted residual
  /**
     a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
   */
  virtual void restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                LevelData<FArrayBox>&       a_phiFine,
                                const LevelData<FArrayBox>& a_rhsFine);

  /// calculate restricted residual
  virtual void restrictResidual(LevelData<FArrayBox>& a_resCoarse,
                                LevelData<FArrayBox>& a_phiFine,
                                const LevelData<FArrayBox>* a_phiCoarse,
                                const LevelData<FArrayBox>& a_rhsFine,
                                bool homogeneous);
  /// Create coarser version of field
  virtual  void createCoarser(LevelData<FArrayBox>&       a_coarse,
                              const LevelData<FArrayBox>& a_fine,
                              bool                        a_ghosted);

  /// Prolong correction
  virtual void prolongIncrement(LevelData<FArrayBox>&       a_phiThisLevel,
                                const LevelData<FArrayBox>& a_correctCoarse);

  /// Restrict solution to coarser level
  virtual void restrictR(LevelData<FArrayBox>& a_phiCoarse, const LevelData<FArrayBox>& a_phiFine);

  /// Apply operator
  void applyOpMg(LevelData<FArrayBox>& a_lhs,
                 LevelData<FArrayBox>& a_phi, LevelData<FArrayBox>* a_phiCoarse, bool a_homogeneous);

  /*@}*/

  /**
     \name AMRNonLinearMultiCompOp functions */
  /*@{*/

  /// For tga stuff
  virtual void setAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta);

  /// Also calls reset lambda
  virtual void setCoefs(const RefCountedPtr<LevelData<FArrayBox> >& a_aCoef,
                        const RefCountedPtr<LevelData<FluxBox>   >& a_bCoef,
                        const Real& a_alpha,
                        const Real& a_beta);

  /// Should be called before the relaxation parameter is needed.
  virtual void resetLambda();

  /// Recalculate lambda for nonlinear calculations (not actually current used)
  void resetLambdaNonLinear(const LevelData<FArrayBox>& a_phi);

  /// Compute lambda once alpha, aCoef, beta, bCoef are defined
  virtual void computeLambda();

  /// Increment flux registers and calculate reflux
  virtual void reflux(const LevelData<FArrayBox>&        a_phiFine,
                      const LevelData<FArrayBox>&        a_phi,
                      LevelData<FArrayBox>&              a_residual,
                      AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /// getFlux function which matches interface to AMRPoissonOp
  /** assumes we want to use member-data bCoef, then calls
        second getFlux function
        Note that a_data is the thing we take the gradient of, i.e. we calculate \f$ \chi \nabla \phi \f$ where \f$ \phi \f$ is a_data.
        a_data is NOT the bulk concentration or enthalpy
   */
  virtual void getFlux(FluxBox&                    a_flux,
                       const LevelData<FArrayBox>& a_data,
                       const Box&                  a_grid,
                       const DataIndex&            a_dit,
                       Real                        a_scale)
  {
    const FluxBox& bCoef = (*m_bCoef)[a_dit];
    getFlux(a_flux,a_data,bCoef,a_grid,a_dit,a_scale);
  }

  /// get diffusive flux
  virtual void getFlux(FluxBox&                    a_flux,
                       const LevelData<FArrayBox>& a_data,
                       const FluxBox&              a_bCoef,
                       const Box&                  a_grid,
                       const DataIndex&            a_dit,
                       Real                        a_scale)
  {

    a_flux.define(a_grid, a_data.nComp());
    const FArrayBox& data = a_data[a_dit];

    // a_data is Enthalpy, Bulk Concentration etc.
    // need to convert this to temperature, liquid concentration to calculate diffusive flux
    //             LevelData<FArrayBox> derivedVar(a_data.disjointBoxLayout(), a_data.nComp(), a_data.ghostVect());
    // Now do this within getFlux itself
    //             computeDerivedVar(derivedVar, a_data);
    //       const FArrayBox& data = derivedVar[a_dit];

    for (int idir=0; idir<SpaceDim; idir++)
    {
      getFlux(a_flux[idir], data, a_bCoef, a_flux[idir].box(), idir, a_dit, 1);
      a_flux[idir] *= a_scale;
    }

  }

  /*@}*/
  /// Has finer operator changed?
  //! This is called on multigrid operators when their AMR operators
  //! are altered.
  void finerOperatorChanged(const MGLevelOp<LevelData<FArrayBox> >& a_operator,
                            int a_coarseningFactor);

  //! Returns identity coefficient data.
  LevelData<FArrayBox>& identityCoef()
      {
    return *m_aCoef;
      }

  //! Sets up a model that modifies b coefficient data when the operator's
  //! time is set.
  //! \param a_bCoefInterpolator A CoefficientInterpolator that will be used
  //!                            to compute the b coefficient at specific
  //!                            times.
  void setBCoefInterpolator(RefCountedPtr<CoefficientInterpolatorLinearFace>& a_bCoefInterpolator)
  {
    m_bCoefInterpolator = a_bCoefInterpolator;
  }


  //! Sets up a model that modifies a coefficient data when the operator's
  //! time is set.
  //! \param a_aCoefInterpolator A CoefficientInterpolator that will be used
  //!                            to compute the a coefficient at specific
  //!                            times.
  void setACoefInterpolator(RefCountedPtr<CoefficientInterpolatorLinear>& a_aCoefInterpolator)
  {
    m_aCoefInterpolator = a_aCoefInterpolator;
  }

  //! Returns the B coefficient.
  LevelData<FluxBox>& BCoef()
      {
    return *m_bCoef;
      }

  /// Allows access to the B coefficient interpolator.
  RefCountedPtr<CoefficientInterpolatorLinearFace> BCoefInterpolator()
      {
    return m_bCoefInterpolator;
      }

  /// Allows access to the A coefficient interpolator.
  RefCountedPtr<CoefficientInterpolatorLinear> ACoefInterpolator()
       {
    return m_aCoefInterpolator;
       }

  //! Sets the time centering of the operator. This interpolates b coefficient
  //! data at the given time if an interpolator is set.
  void setTime(Real a_time);

  /// Compute temperature, liquid concentration from enthalpy, bulk concentration etc.
  void computeDiffusedVar(LevelData<FArrayBox>& liquidConc,  const LevelData<FArrayBox>& bulkConc,
                         bool a_homogeneous=false);

  /// Compute temperature, liquid concentration from enthalpy, bulk concentration etc.
  virtual void computeDiffusedVar(FArrayBox& a_calculatedVar,  const FArrayBox& a_primaryVar,
                         const DataIndex dit, bool a_homogeneous=false);

  /// Identity operator spatially varying coefficient storage (cell-centered) --- if you change this call resetLambda()
  RefCountedPtr<LevelData<FArrayBox> > m_aCoef,

  /// Which of enthalpy or bulk concentration we're not solving for
    m_secondaryVar,

    /// Enthalpy solidus
    m_enthalpySolidus,

    /// Enthalpy liquidus
    m_enthalpyLiquidus,

    /// Eutectic enthalpy
    m_enthalpyEutectic;



/// Boundary conditions for the diffused field
  BCHolder m_diffusedVarBC;

  /// Should we apply BCs explicitly to the diagnostic variables?
  bool m_apply_bcs_to_diagnostic_var;



  /// Laplacian operator spatially varying coefficient storage (face-centered) --- if you change this call resetLambda()
  RefCountedPtr<LevelData<FluxBox> > m_bCoef;

  /// Reciprocal of the diagonal entry of the operator matrix
  LevelData<FArrayBox> m_lambda;

  /// Physical parameters for problems
  MushyLayerParams* m_params;

  /// true if we're running in FAS mode (much quicker for nonlinear problems)
  bool m_FAS;

  /// Do some dodgy stuff here to try and make the code quicker. Turned off by default.
  bool m_superOptimised;



protected:
  /// Intvects on coarse-fine interfaces
  LayoutData<CFIVS> m_loCFIVS[SpaceDim];

  /// Intvects on coarse-fine interfaces
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];

  /// Interpolator for b coefficient data.
  RefCountedPtr<CoefficientInterpolatorLinearFace> m_bCoefInterpolator;

  /// Interpolator for a coefficient data.
  RefCountedPtr<CoefficientInterpolatorLinear> m_aCoefInterpolator;

  /// Current time.
  Real m_time;

  /// Does the relaxation coefficient need to be reset?
  bool m_lambdaNeedsResetting;

  /// Gauss-seidel relaxation
  virtual void levelGSRB(LevelData<FArrayBox>&       a_phi,
                         const LevelData<FArrayBox>& a_rhs);

  /// Multicolor relaxation
  virtual void levelMultiColor(LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_rhs);

  /// Gauss-seidel relaxation
  virtual void looseGSRB(LevelData<FArrayBox>&       a_phi,
                         const LevelData<FArrayBox>& a_rhs);

  /// Gauss-seidel relaxation
  virtual void overlapGSRB(LevelData<FArrayBox>&       a_phi,
                           const LevelData<FArrayBox>& a_rhs);

  /// Gauss-seidel relaxation
  virtual void levelGSRBLazy(LevelData<FArrayBox>&       a_phi,
                             const LevelData<FArrayBox>& a_rhs);

  /// Jacobi relaxation
  virtual void levelJacobi(LevelData<FArrayBox>&       a_phi,
                           const LevelData<FArrayBox>& a_rhs);

  /// computes flux over face-centered a_facebox.
  virtual void getFlux(FArrayBox&       a_flux,
                       const FArrayBox& a_data,
                       const FluxBox&   a_bCoef,
                       const Box&       a_facebox,
                       int              a_dir,
                       const DataIndex& a_dit,
                       int              a_ref = 1) const ;
#if 0
  /// utility function which computes operator after all bc's have been set
  void computeOperatorNoBCs(LevelData<FArrayBox>&       a_lhs,
                            const LevelData<FArrayBox>& a_phi);
#endif
};

///Factory for nonlinear variable coefficient operator
/**
   Factory to create AMRNonLinearMultiCompOps
 */
class AMRNonLinearMultiCompOpFactory: public AMRLevelOpFactory<LevelData<FArrayBox> >
{
public:
  AMRNonLinearMultiCompOpFactory();

  virtual ~AMRNonLinearMultiCompOpFactory()
  {
  }

  /// Define object
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
  void define(const ProblemDomain&                           a_coarseDomain,
              const Vector<DisjointBoxLayout>&               a_grids,
              const Vector<int>&                             a_refRatios,
              const Real&                                    a_coarsedx,
              BCHolder                                       a_bc,
              const Real&                                    a_alpha,
              Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
              const Real&                                    a_beta,
              Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef,
              Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_enthalpySolidus,
              Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_enthalpyLiquidus,
              Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_enthalpyEutectic,
              MushyLayerParams*			   a_params,
              BCHolder  a_derivedVarBC,
              int a_relaxMode,
              EdgeVelBCHolder a_porosityEdgeBC,
              bool a_apply_bcs_to_diagnostic_var = true);


  /// New MGLevelOp
  virtual MGLevelOp<LevelData<FArrayBox> >* MGnewOp(const ProblemDomain& a_FineindexSpace,
                                                    int                  a_depth,
                                                    bool                 a_homoOnly = true);

  /// New AMRLevelOp
  virtual AMRLevelOp<LevelData<FArrayBox> >* AMRnewOp(const ProblemDomain& a_indexSpace);

  /// refinement ratio to finer level
  virtual int refToFiner(const ProblemDomain& a_domain) const;

  /// Set primary field
  void setPrimaryVarPtr(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_primaryVarPtr);

  /// Whether we should use the primary variable for solves
  void setUsePrimaryVar(bool a_usePrimaryVar);

  /// Set the boundary conditions m_bc
  void setBC(BCHolder& a_bc);

  /// Set whether or not to turn on the 'super optimised' flag which does some dodgy stuff to try and speed up calculations (and is turned off by default).
  void setSuperOptimised(bool a_val);

  /// Type of coefficient averaging
  int m_coefficient_average_type;

private:

  /// Set default values
  void setDefaultValues();

  /// Problem domains across AMR hierarchy
  Vector<ProblemDomain>     m_domains;

  /// DisjointBoxLayouts across AMR hierarchy
  Vector<DisjointBoxLayout> m_boxes;

  /// Cell spacings
  Vector<Real> m_dx;

  /// refinement to next coarser level
  Vector<int>  m_refRatios;

  /// Physical boundary condition
  BCHolder m_bc, m_diffusedVarBC;

  /// Apply BCs explicitly to the diagnostic variables?
  bool m_apply_bcs_to_diagnostic_var;

  /// Porosity boundary condition
  EdgeVelBCHolder m_porosityEdgeBC;

  /// identity coefficient
  Real m_alpha;

  /// laplacian coefficient
  Real m_beta;

  /// Spatially varying identity coefficient
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_aCoef,

  /// Which of enthalpy or bulk concentration we're not solving for
  m_secondaryVar,

  /// Enthalpy solidus
  m_enthalpySolidus,

  /// Enthalpy liquidus
  m_enthalpyLiquidus,

  /// Eutectic enthalpy
  m_enthalpyEutectic;

  /// Spatially varying diffusion coefficient
  /**
   * E.g. \f$ - \chi / Le \f$
   */
  Vector<RefCountedPtr<LevelData<FluxBox> > >   m_bCoef;

  /// Physical parameters
  MushyLayerParams* m_params;

  /// Lambda. Why is this here?
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_lambda;

  /// Copiers for exchange calls
  Vector<Copier>   m_exchangeCopiers;

  /// Coarse fine regions
  Vector<CFRegion> m_cfregion;

  /// Object to compute enthalpy variables
  //  ComputeEnthalpyVars m_computeEnthalpyVars;

  /// Type of relaxation
  /**
   * 4 is for jacobi, 1 is for GSRB
   */
  int m_relaxMode;

  /// True if we're in FAS mode
  bool m_FAS;

  int m_numComp;

  /// If turned on, we'll try and do some possibly dodgy stuff to speed up the code.
  /**
   * Turned off by default, use with care.
   */
  bool m_superOptimised;

};

#include "NamespaceFooter.H"
#endif
