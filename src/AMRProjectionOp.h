/*
 * AMRProjectionOp.h
 *
 *  Created on: 16 Nov 2018
 *      Author: parkinsonjl
 */

#ifndef SRC_AMRPROJECTIONOP_H_
#define SRC_AMRPROJECTIONOP_H_

#include "VCAMRPoissonOp2.H"
#include "CoarseAverage.H"

#include "NamespaceHeader.H"

class AMRProjectionOp : public VCAMRPoissonOp2
{
public:
  AMRProjectionOp ();
  virtual
  ~AMRProjectionOp ();

  virtual void prolongIncrement(LevelData<FArrayBox>&       a_phiThisLevel,
                                 const LevelData<FArrayBox>& a_correctCoarse);

  virtual void resetLambda();
};

// Would have liked to base this off VCAMRPoisosnOp2Factory,
// but that's all private so we have to do it ourselves
class AMRProjectionOpFactory: public AMRLevelOpFactory<LevelData<FArrayBox> >
{
public:
  AMRProjectionOpFactory();

//  AMRProjectionOpFactory(int a_average_type);

  virtual ~AMRProjectionOpFactory()
  {
  }

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
              int a_averaging_type = CoarseAverage::arithmetic);

  //! Defines a factory for VCAMRPoissonOp2 which allows the operators to
  //! allocate their own coefficient data. \f$\alpha\f$ and \f$\beta\f$
  //! coefficients are initialized to 1.
  //! \param a_coarseDomain The domain at the coarsest level.
  //! \param a_grids The disjoint box layouts for the various refinement levels.
  //! \param a_refRatios The refinement ratios between levels.
  //! \param a_coarsedx The grid spacing at the coarsest level.
  //! \param a_bc The boundary condition imposed on the solution.
  //! \param a_ghostVect The ghost stencil to use in the created coefficient data.
  void define(const ProblemDomain&                           a_coarseDomain,
              const Vector<DisjointBoxLayout>&               a_grids,
              const Vector<int>&                             a_refRatios,
              const Real&                                    a_coarsedx,
              BCHolder                                       a_bc,
              const IntVect&                                 a_ghostVect,
              int a_averaging_type = CoarseAverage::arithmetic);

  ///
  virtual MGLevelOp<LevelData<FArrayBox> >* MGnewOp(const ProblemDomain& a_FineindexSpace,
                                                    int                  a_depth,
                                                    bool                 a_homoOnly = true);

  ///
  virtual AMRLevelOp<LevelData<FArrayBox> >* AMRnewOp(const ProblemDomain& a_indexSpace);

  ///
  virtual int refToFiner(const ProblemDomain& a_domain) const;

  int m_coefficient_average_type;
  int m_relaxMode;

private:
  void setDefaultValues();

  Vector<ProblemDomain>     m_domains;
  Vector<DisjointBoxLayout> m_boxes;

  Vector<Real> m_dx;
  Vector<int>  m_refRatios; // refinement to next coarser level

  BCHolder m_bc;

  Real m_alpha;
  Real m_beta;

  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_aCoef;
  Vector<RefCountedPtr<LevelData<FluxBox> > >   m_bCoef;

  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_lambda;

  Vector<Copier>   m_exchangeCopiers;
  Vector<CFRegion> m_cfregion;


};

#include "NamespaceFooter.H"
#endif /* SRC_AMRPROJECTIONOP_H_ */
