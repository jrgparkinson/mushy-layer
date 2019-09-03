/*
 * AMRScalarDiffusionOp.h
 *
 *  Created on: 23 Aug 2019
 *      Author: parkinsonjl
 */

#ifndef SRC_AMRSCALARDIFFUSIONOP_H_
#define SRC_AMRSCALARDIFFUSIONOP_H_


#include "VCAMRPoissonOp2.H"
#include "CoarseAverage.H"

#include "NamespaceHeader.H"



/// Operator for doing diffusion of a scalar quantity
/**
 * Derived from VCAMRPoissonOp2. We evaluate
 *
 * \f$ D \nabla \cdot \chi \nabla (\xi/\chi) \f$
 *
 * where \f$D\f$ is the diffusion coefficient, \f$ \xi\f$ is the bulk concentration of whatever we're diffusing, and
 * \f$ \chi\f$ is the porosity.
 *
 * This class is unfinished. Just using VCAMRPOissonOp2. Writing this class turned out to be more complicated than expected.
 */
class AMRScalarDiffusionOp : public VCAMRPoissonOp2
{
public:
  AMRScalarDiffusionOp ();

  virtual
  ~AMRScalarDiffusionOp ();

};

class AMRScalarDiffusionOpFactory: public VCAMRPoissonOp2Factory
{
public:
  AMRScalarDiffusionOpFactory()
{
    // default values:
    setDefaultValues();
}


  virtual ~AMRScalarDiffusionOpFactory()
  {
  }

  void define(const ProblemDomain&                           a_coarseDomain,
                const Vector<DisjointBoxLayout>&               a_grids,
                const Vector<int>&                             a_refRatios,
                const Real&                                    a_coarsedx,
                BCHolder                                       a_bc,
                const Real&                                    a_alpha,
                Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                const Real&                                    a_beta,
                Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef,
                Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_porosity);
//  {
//
//    define(a_coarseDomain, a_grids, a_refRatios, a_coarsedx, a_bc, a_alpha, a_aCoef, a_beta, a_bCoef);
//    m_porosity = a_porosity;
//  }

protected:
  void setDefaultValues()
  {
    // Default to Laplacian operator
    m_alpha = 0.0;
    m_beta = -1.0;

    m_coefficient_average_type = CoarseAverage::arithmetic;
  }

  Vector<ProblemDomain>     m_domains;
  Vector<DisjointBoxLayout> m_boxes;

  Vector<Real> m_dx;
  Vector<int>  m_refRatios; // refinement to next coarser level

  BCHolder m_bc;

  Real m_alpha;
  Real m_beta;

  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_aCoef;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_porosity;
  Vector<RefCountedPtr<LevelData<FluxBox> > >   m_bCoef;

  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_lambda;

  Vector<Copier>   m_exchangeCopiers;
  Vector<CFRegion> m_cfregion;
};

#include "NamespaceFooter.H"
#endif /* SRC_AMRSCALARDIFFUSIONOP_H_ */
