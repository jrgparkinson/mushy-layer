/*
 * AMRScalarDiffusionOp.cpp
 *
 *  Created on: 23 Aug 2019
 *      Author: parkinsonjl
 */

#include "AMRScalarDiffusionOp.h"
#include "VCAMRPoissonOp2.H"

void AMRScalarDiffusionOpFactory::define(
    const ProblemDomain &a_coarseDomain,
    const Vector<DisjointBoxLayout> &a_grids, const Vector<int> &a_refRatios,
    const Real &a_coarsedx, BCHolder a_bc, const Real &a_alpha,
    Vector<RefCountedPtr<LevelData<FArrayBox>>> &a_aCoef, const Real &a_beta,
    Vector<RefCountedPtr<LevelData<FluxBox>>> &a_bCoef,
    Vector<RefCountedPtr<LevelData<FArrayBox>>> &a_porosity)
{
  //    VCAMRPoissonOp2Factory* vcop = <VCAMRPoissonOp2Factory*>(this);
  VCAMRPoissonOp2Factory::define(a_coarseDomain, a_grids, a_refRatios,
                                 a_coarsedx, a_bc, a_alpha, a_aCoef, a_beta,
                                 a_bCoef);
  m_porosity = a_porosity;
}
