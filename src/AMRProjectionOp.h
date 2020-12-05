/*
 * AMRProjectionOp.h
 *
 *  Created on: 16 Nov 2018
 *      Author: parkinsonjl
 */

#ifndef SRC_AMRPROJECTIONOP_H_
#define SRC_AMRPROJECTIONOP_H_

#include "CoarseAverage.H"
#include "VCAMRPoissonOp2.H"

#include "NamespaceHeader.H"

/// Operator for doing projector with a variable coefficient.
/**
 * This is essentially the same as VCAMRPoissonOp2, but allows us to make
 * changes if we want.
 */
class AMRProjectionOp : public VCAMRPoissonOp2
{
public:
  AMRProjectionOp();
  virtual ~AMRProjectionOp();

  /// Prolong operation for multigrid
  /**
   * Transfer data from a coarse grid to a fine grid, doing interpolation as
   * required
   */
  virtual void prolongIncrement(LevelData<FArrayBox> &a_phiThisLevel,
                                const LevelData<FArrayBox> &a_correctCoarse);

  /// Reset the relaxation coefficient
  virtual void resetLambda();
};

/// Factory for creating AMRProjectionOp's
/**
 * Would have liked to base this off VCAMRPoissonOp2Factory,
 *  but that's all private so we have to do it ourselves
 */

class AMRProjectionOpFactory : public AMRLevelOpFactory<LevelData<FArrayBox>>
{
public:
  AMRProjectionOpFactory();

  //  AMRProjectionOpFactory(int a_average_type);

  virtual ~AMRProjectionOpFactory() {}

  ///
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
     a_bCoef is the laplacian spatially varying coefficient
     a_averaging_type is the method to use for averaging cell to face centred
     variables
  */
  void define(const ProblemDomain &a_coarseDomain,
              const Vector<DisjointBoxLayout> &a_grids,
              const Vector<int> &a_refRatios, const Real &a_coarsedx,
              BCHolder a_bc, const Real &a_alpha,
              Vector<RefCountedPtr<LevelData<FArrayBox>>> &a_aCoef,
              const Real &a_beta,
              Vector<RefCountedPtr<LevelData<FluxBox>>> &a_bCoef,
              int a_averaging_type = CoarseAverage::arithmetic);

  //! Defines a factory for VCAMRPoissonOp2 which allows the operators to
  //! allocate their own coefficient data. \f$\alpha\f$ and \f$\beta\f$
  //! coefficients are initialized to 1.
  //! \param a_coarseDomain The domain at the coarsest level.
  //! \param a_grids The disjoint box layouts for the various refinement levels.
  //! \param a_refRatios The refinement ratios between levels.
  //! \param a_coarsedx The grid spacing at the coarsest level.
  //! \param a_bc The boundary condition imposed on the solution.
  //! \param a_ghostVect The ghost stencil to use in the created coefficient
  //! data. \param a_averaging_type The method to use for averaging cell to face
  //! centred variables
  void define(const ProblemDomain &a_coarseDomain,
              const Vector<DisjointBoxLayout> &a_grids,
              const Vector<int> &a_refRatios, const Real &a_coarsedx,
              BCHolder a_bc, const IntVect &a_ghostVect,
              int a_averaging_type = CoarseAverage::arithmetic);

  /// Make a new multigrid operator for a given domain
  virtual MGLevelOp<LevelData<FArrayBox>> *
  MGnewOp(const ProblemDomain &a_FineindexSpace, int a_depth,
          bool a_homoOnly = true);

  /// Make a new AMR multirid operator
  virtual AMRLevelOp<LevelData<FArrayBox>> *
  AMRnewOp(const ProblemDomain &a_indexSpace);

  /// Returns the refinement ratio to the next finer level
  virtual int refToFiner(const ProblemDomain &a_domain) const;

  /// How to do coefficient averaging (arithmetic, harmonic, etc. )
  int m_coefficient_average_type;

  /// How to do relaxation (gauss-seidel, jacobi, etc.)
  int m_relaxMode;

private:
  /// Assign default values to this objects parameters
  void setDefaultValues();

  /// Problem domain on each level this operator acts on
  Vector<ProblemDomain> m_domains;

  /// Grids on each level of refinement
  Vector<DisjointBoxLayout> m_boxes;

  /// Grid spacing on each level of refinement
  Vector<Real> m_dx;

  /// refinement between each level and the coarser level
  Vector<int> m_refRatios;

  /// Boundary conditions for this operator
  BCHolder m_bc;

  /// Coefficient of the diagonal term
  Real m_alpha;

  /// Coefficient for the grad-div term
  Real m_beta;

  /// Spatially varying coefficient of the diagonal term
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_aCoef;

  /// Spatially varying coefficient b inside the grad( b . div) phi term
  Vector<RefCountedPtr<LevelData<FluxBox>>> m_bCoef;

  /// Spatially varying relaxation parameters
  Vector<RefCountedPtr<LevelData<FArrayBox>>> m_lambda;

  /// Objects for copying data between boxes on the same level
  Vector<Copier> m_exchangeCopiers;

  /// Objects that represents the edge region around disjoint box layouts
  Vector<CFRegion> m_cfregion;
};

#include "NamespaceFooter.H"
#endif /* SRC_AMRPROJECTIONOP_H_ */
