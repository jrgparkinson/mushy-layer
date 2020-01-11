#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FORT_PROTO.H"
#include "BoxIterator.H"
#include "AverageF_F.H"
#include "InterpF_F.H"
#include "LayoutIterator.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "AMRMultiGrid.H"
#include "Misc.H"

#include "AMRPoissonOpF_F.H"

#include "AMRNonLinearMultiCompOp.H"
#include "AMRNonLinearMultiCompOpF_F.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"
#include "DebugOut.H"
#include "MushyLayerUtils.H"
#include "EnthalpyVariablesF_F.H"

#include "CoefficientInterpolatorLinear.H"

#include "NamespaceHeader.H"

void AMRNonLinearMultiCompOp::computeDiffusedVar(LevelData<FArrayBox>& a_diffusedVar, const LevelData<FArrayBox>& a_phi, bool a_homogeneous)
{
  CH_TIME("AMRNonLinearMultiCompOp::computeDerivedVar");

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
  {
    computeDiffusedVar(a_diffusedVar[dit], a_phi[dit], dit(), a_homogeneous);
  }


  {
    CH_TIME("AMRNonLinearMultiCompOp::computeDiffusedVar::exchange");
    a_diffusedVar.exchange(a_diffusedVar.interval(), m_exchangeCopier);
    // Stick a corner copier in here too?
    //	const DisjointBoxLayout dbl = a_diffusedVar.disjointBoxLayout();
    //	CornerCopier cornerCopy(dbl, dbl, m_domain, a_diffusedVar.ghostVect(), true);
    //	a_calculatedVar.exchange(a_diffusedVar.interval(), cornerCopy);
  }

}
//
void AMRNonLinearMultiCompOp::computeDiffusedVar(FArrayBox& a_diffusedVar,
                                                 const FArrayBox& a_phi, const DataIndex dit,
                                                 bool a_homogeneous)
{
  // Another attempt to stop solver fails - recalculate solidus/eutectic/liquidus each time

  FArrayBox& Hs = (*m_enthalpySolidus)[dit];
  FArrayBox& He = (*m_enthalpyEutectic)[dit];
  FArrayBox& Hl = (*m_enthalpyLiquidus)[dit];

  Box region = a_diffusedVar.box();
  region &= Hs.box();

  FORT_CALCULATE_BOUNDING_ENERGY( CHF_CONST_FRA(a_phi),
                       CHF_FRA(Hs),
                       CHF_FRA(He),
                       CHF_FRA(Hl),
                       CHF_BOX(region),
                       CHF_CONST_REAL(m_params->compositionRatio),
                       CHF_CONST_REAL(m_params->waterDistributionCoeff),
                       CHF_CONST_REAL(m_params->specificHeatRatio),
                       CHF_CONST_REAL(m_params->stefan),
                       CHF_CONST_REAL(m_params->thetaEutectic),
                       CHF_CONST_REAL(m_params->ThetaEutectic));

  // Apply BCs to H-C, then compute T, S_l, porosity BC from these values
  // a_phi is const so can't apply BCs here - assume they're already set
//  m_bc(a_phi, m_domain.domainBox(), m_domain, m_dx, a_homogeneous);

  FORT_CALCULATE_T_CL(CHF_FRA(a_diffusedVar),
                      CHF_CONST_FRA(a_phi),
                      CHF_CONST_FRA(Hs),
                      CHF_CONST_FRA(He),
                      CHF_CONST_FRA(Hl),
                      CHF_BOX(region),
                      CHF_CONST_REAL(m_params->compositionRatio),
                      CHF_CONST_REAL(m_params->waterDistributionCoeff),
                      CHF_CONST_REAL(m_params->specificHeatRatio),
                      CHF_CONST_REAL(m_params->stefan),
                      CHF_CONST_REAL(m_params->thetaEutectic),
                      CHF_CONST_REAL(m_params->ThetaEutectic));

  if (m_apply_bcs_to_diagnostic_var)
  {
    m_diffusedVarBC(a_diffusedVar, m_domain.domainBox(), m_domain, m_dx, a_homogeneous);
  }

//  this->m_bc
}

void AMRNonLinearMultiCompOp::residualI(LevelData<FArrayBox>&       a_lhs,
                                        const LevelData<FArrayBox>& a_phi,
                                        const LevelData<FArrayBox>& a_rhs,
                                        bool                        a_homogeneous)
{
  CH_TIME("AMRNonLinearMultiCompOp::residualI");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  LevelData<FArrayBox> derivedVar(phi.disjointBoxLayout(), phi.nComp(), phi.ghostVect());

  Real dx = m_dx;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  {
    CH_TIME("AMRNonLinearMultiCompOp::residualIBC");

    for (dit.begin(); dit.ok(); ++dit)
    {
      m_bc(phi[dit], dbl[dit()],m_domain, dx, a_homogeneous);
    }
  }


  phi.exchange(phi.interval(), m_exchangeCopier);

  //  computeDiffusedVar(derivedVar, phi, a_homogeneous);

  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& region = dbl[dit()];
    const FluxBox& thisBCoef = (*m_bCoef)[dit];

    //FArrayBox derivedVar(phi[dit]);
    // Is this quicker than the level data version?
    computeDiffusedVar(derivedVar[dit], phi[dit], dit(), a_homogeneous);



#if CH_SPACEDIM == 1
    FORT_NLVCCOMPUTERES1D
#elif CH_SPACEDIM == 2
    FORT_NLVCCOMPUTERES2D
    //		FORT_VCCOMPUTERES2D
#elif CH_SPACEDIM == 3
    FORT_NLVCCOMPUTERES3D
#else
    //		This_will_not_compile!
#endif
    (CHF_FRA(a_lhs[dit]),
     CHF_CONST_FRA(phi[dit]),
     CHF_CONST_FRA(derivedVar[dit]),
     CHF_CONST_FRA(a_rhs[dit]),
     CHF_CONST_REAL(m_alpha),
     CHF_CONST_FRA((*m_aCoef)[dit]),
     CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
     CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
     CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
     CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
     This_will_not_compile!
#endif
     CHF_BOX(region),
     CHF_CONST_REAL(m_dx));

  } // end loop over boxes


}

/**************************/
// this preconditioner first initializes phihat to (IA)phihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSRB
void AMRNonLinearMultiCompOp::preCond(LevelData<FArrayBox>&       a_phi,
                                      const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearMultiCompOp::preCond");

  // diagonal term of this operator in:
  //
  //       alpha * a(i)
  //     + beta  * sum_over_dir (b(i-1/2*e_dir) + b(i+1/2*e_dir)) / (dx*dx)
  //
  // The inverse of this is our initial multiplier.

  int ncomp = a_phi.nComp();

  CH_assert(m_lambda.isDefined());
  CH_assert(a_rhs.nComp()    == ncomp);
  CH_assert(m_bCoef->nComp() == ncomp);

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  // don't need to use a Copier -- plain copy will do
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    // also need to average and sum face-centered bCoefs to cell-centers
    Box gridBox = a_rhs[dit].box();

    // approximate inverse
    a_phi[dit].copy(a_rhs[dit]);
    a_phi[dit].mult(m_lambda[dit], gridBox, 0, 0, ncomp);
  }

  relax(a_phi, a_rhs, 2);
}

void AMRNonLinearMultiCompOp::applyOpMg(LevelData<FArrayBox>& a_lhs, LevelData<FArrayBox>& a_phi,
                                        LevelData<FArrayBox>* a_phiCoarse, bool a_homogeneous)
{
  // Do CF stuff if we have a coarser level that's not just a single grid cell
  if (a_phiCoarse != nullptr)
  {
    const ProblemDomain& probDomain = a_phiCoarse->disjointBoxLayout().physDomain();
    const Box& domBox = probDomain.domainBox();
    //    IntVect hi = domBox.b
    if (domBox.bigEnd() != domBox.smallEnd())
    {
      m_interpWithCoarser.coarseFineInterp(a_phi, *a_phiCoarse);
    }
  }

  applyOpI(a_lhs, a_phi, a_homogeneous);
}

void AMRNonLinearMultiCompOp::applyOpI(LevelData<FArrayBox>&      a_lhs,
                                       const LevelData<FArrayBox>& a_phi,
                                       bool                        a_homogeneous )
{
  CH_TIME("AMRNonLinearMultiCompOp::applyOpI");
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  Real dx = m_dx;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
  {
    m_bc(phi[dit], dbl[dit()],m_domain, dx, a_homogeneous);
  }

  applyOpNoBoundary2(a_lhs, a_phi,a_homogeneous );
}

void AMRNonLinearMultiCompOp::applyOpNoBoundary2(LevelData<FArrayBox>&      a_lhs,
                                                 const LevelData<FArrayBox>& a_phi,
                                                 bool                        a_homogeneous )
{
  CH_TIME("AMRNonLinearMultiCompOp::applyOpNoBoundary");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  LevelData<FArrayBox> derivedVar(phi.disjointBoxLayout(), phi.nComp(), phi.ghostVect());

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();

  phi.exchange(phi.interval(), m_exchangeCopier);

  computeDiffusedVar(derivedVar, phi, a_homogeneous);

  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& region = dbl[dit()];
    const FluxBox& thisBCoef = (*m_bCoef)[dit];

#if CH_SPACEDIM == 1
    FORT_NLVCCOMPUTEOP1D
#elif CH_SPACEDIM == 2
    FORT_NLVCCOMPUTEOP2D
    //		FORT_VCCOMPUTEOP2D
#elif CH_SPACEDIM == 3
    FORT_NLVCCOMPUTEOP3D
#else
    //		This_will_not_compile!
#endif
    (CHF_FRA(a_lhs[dit]),
     CHF_CONST_FRA(phi[dit]),
     CHF_CONST_FRA(derivedVar[dit]),
     CHF_CONST_REAL(m_alpha),
     CHF_CONST_FRA((*m_aCoef)[dit]),
     CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
     CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
     CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
     CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
     This_will_not_compile!
#endif
     CHF_BOX(region),
     CHF_CONST_REAL(m_dx));
    //		int temp=0;

  } // end loop over boxes
}

void AMRNonLinearMultiCompOp::applyOpNoBoundary(LevelData<FArrayBox>&      a_lhs,
                                                const LevelData<FArrayBox>& a_phi)
{
  CH_TIME("AMRNonLinearMultiCompOp::applyOpNoBoundary");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  LevelData<FArrayBox> derivedVar(phi.disjointBoxLayout(), phi.nComp(), phi.ghostVect());

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();

  phi.exchange(phi.interval(), m_exchangeCopier);

  computeDiffusedVar(derivedVar, phi);

  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& region = dbl[dit()];
    const FluxBox& thisBCoef = (*m_bCoef)[dit];

#if CH_SPACEDIM == 1
    FORT_NLVCCOMPUTEOP1D
#elif CH_SPACEDIM == 2
    FORT_NLVCCOMPUTEOP2D
    //		FORT_VCCOMPUTEOP2D
#elif CH_SPACEDIM == 3
    FORT_NLVCCOMPUTEOP3D
#else
    //		This_will_not_compile!
#endif
    (CHF_FRA(a_lhs[dit]),
     CHF_CONST_FRA(phi[dit]),
     CHF_CONST_FRA(derivedVar[dit]),
     CHF_CONST_REAL(m_alpha),
     CHF_CONST_FRA((*m_aCoef)[dit]),
     CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
     CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
     CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
     CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
     This_will_not_compile!
#endif
     CHF_BOX(region),
     CHF_CONST_REAL(m_dx));
  } // end loop over boxes
}

// Re implement here in a slightly different way to AMRPoissonOp
void AMRNonLinearMultiCompOp::createCoarser(LevelData<FArrayBox>&       a_coarse,
                                            const LevelData<FArrayBox>& a_fine,
                                            bool                        a_ghosted)
{
  CH_TIME("AMRPoissonOp::createCoarser");

  // CH_assert(!a_ghosted);
  IntVect ghost = a_fine.ghostVect();

  const DisjointBoxLayout& fineGrids =  a_fine.disjointBoxLayout();
  DisjointBoxLayout crseGrids;

  CH_assert(fineGrids.coarsenable(2));
  coarsen(crseGrids,fineGrids, 2); //multigrid, so coarsen by 2

  a_coarse.define(crseGrids, a_fine.nComp(), ghost);
}

void AMRNonLinearMultiCompOp::restrictResidual(LevelData<FArrayBox>& a_resCoarse,
                                               LevelData<FArrayBox>& a_phiFine,
                                               const LevelData<FArrayBox>* a_phiCoarse,
                                               const LevelData<FArrayBox>& a_rhsFine,
                                               bool homogeneous)
{
  CH_TIME("AMRNonLinearMultiCompOp::restrictResidual");

  if (m_FAS)
  {
    if (a_phiCoarse != nullptr)
    {
      // Do inhomo CF bcs
      m_interpWithCoarser.coarseFineInterp(a_phiFine, *a_phiCoarse);
    }
  }
  else
  {
    homogeneousCFInterp(a_phiFine);
  }


  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
  {
    FArrayBox& phi = a_phiFine[dit];
    m_bc(phi, dblFine[dit()], m_domain, m_dx, homogeneous);
  }

  a_phiFine.exchange(a_phiFine.interval(), m_exchangeCopier);

  LevelData<FArrayBox> derivedVar(a_phiFine.disjointBoxLayout(), a_phiFine.nComp(), a_phiFine.ghostVect());
  computeDiffusedVar(derivedVar, a_phiFine, homogeneous);

  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
  {
    FArrayBox&       phi = a_phiFine[dit];
    FArrayBox&       Cl = derivedVar[dit];
    const FArrayBox& rhs = a_rhsFine[dit];
    FArrayBox&       res = a_resCoarse[dit];

    const FArrayBox& thisACoef = (*m_aCoef)[dit];
    const FluxBox&   thisBCoef = (*m_bCoef)[dit];

    Box region = dblFine.get(dit());
    const IntVect& iv = region.smallEnd();
    IntVect civ = coarsen(iv, 2);

    res.setVal(0.0);

#if CH_SPACEDIM == 1
    FORT_RESTRICTRESNLVC1D
#elif CH_SPACEDIM == 2
    FORT_RESTRICTRESNLVC2D
    //		FORT_RESTRICTRESVC2D
#elif CH_SPACEDIM == 3
    FORT_RESTRICTRESNLVC3D
#else
    //		This_will_not_compile!
#endif
    (CHF_FRA_SHIFT(res, civ),
     CHF_CONST_FRA_SHIFT(phi, iv),
     CHF_CONST_FRA_SHIFT(Cl, iv),
     CHF_CONST_FRA_SHIFT(rhs, iv),
     CHF_CONST_REAL(m_alpha),
     CHF_CONST_FRA_SHIFT(thisACoef, iv),
     CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
     CHF_CONST_FRA_SHIFT(thisBCoef[0], iv),
#endif
#if CH_SPACEDIM >= 2
     CHF_CONST_FRA_SHIFT(thisBCoef[1], iv),
#endif
#if CH_SPACEDIM >= 3
     CHF_CONST_FRA_SHIFT(thisBCoef[2], iv),
#endif
#if CH_SPACEDIM >= 4
     This_will_not_compile!
#endif
     CHF_BOX_SHIFT(region, iv),
     CHF_CONST_REAL(m_dx));
  }

}

// ---------------------------------------------------------
void AMRNonLinearMultiCompOp::prolongIncrement(LevelData<FArrayBox>&       a_phiThisLevel,
                                               const LevelData<FArrayBox>& a_correctCoarse)
{
  CH_TIME("AMRNonLinearMultiCompOp::prolongIncrement");

  DisjointBoxLayout dbl = a_phiThisLevel.disjointBoxLayout();
  DataIterator dit = a_phiThisLevel.dataIterator();
  int nbox=dit.size();

#pragma omp parallel
  {
#pragma omp for
    for(int ibox = 0; ibox < nbox; ibox++)
    {
      FArrayBox& phi =  a_phiThisLevel[dit[ibox]];
      const FArrayBox& coarse = a_correctCoarse[dit[ibox]];
      Box region = dbl[dit[ibox]];
      const IntVect& iv = region.smallEnd();
      IntVect civ=coarsen(iv, 2);

      // refinement of two as this is multigrid
      int refinement = 2;

      FORT_PROLONG_2(CHF_FRA_SHIFT(phi, iv),
                     CHF_CONST_FRA_SHIFT(coarse, civ),
                     CHF_BOX_SHIFT(region, iv),
                     CHF_CONST_INT(refinement));

          }
  }//end pragma
}

void AMRNonLinearMultiCompOp::restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                               LevelData<FArrayBox>&       a_phiFine,
                                               const LevelData<FArrayBox>& a_rhsFine)
{
  // default implementation is homogeneous
  restrictResidual(a_resCoarse, a_phiFine, nullptr, a_rhsFine, true);
}

void AMRNonLinearMultiCompOp::restrictR(LevelData<FArrayBox>& a_phiCoarse,
                                       const LevelData<FArrayBox>& a_phiFine)
{
  //	a_phiFine.exchange(a_phiFine.interval(), m_exchangeCopier);

  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();

  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
  {
    const FArrayBox&       phiFine = a_phiFine[dit];
    FArrayBox&       phiCoarse = a_phiCoarse[dit];

    Box region = dblFine.get(dit());
    const IntVect& iv = region.smallEnd();
    IntVect civ = coarsen(iv, 2);

    phiCoarse.setVal(0.0);

    FORT_RESTRICT(CHF_FRA_SHIFT(phiCoarse, civ),
                  CHF_CONST_FRA_SHIFT(phiFine, iv),
                  CHF_BOX_SHIFT(region, iv),
                  CHF_CONST_REAL(m_dx));
  }
}
void AMRNonLinearMultiCompOp::setAlphaAndBeta(const Real& a_alpha,
                                              const Real& a_beta)
{
  m_alpha = a_alpha;
  m_beta  = a_beta;

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;
}

void AMRNonLinearMultiCompOp::setCoefs(const RefCountedPtr<LevelData<FArrayBox> >& a_aCoef,
                                       const RefCountedPtr<LevelData<FluxBox  > >& a_bCoef,
                                       const Real&                                 a_alpha,
                                       const Real&                                 a_beta)
{
  m_alpha = a_alpha;
  m_beta  = a_beta;

  m_aCoef = a_aCoef;
  m_bCoef = a_bCoef;

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;
}



void AMRNonLinearMultiCompOp::resetLambda()
{
  CH_TIME("AMRNonLinearMultiCompOp::resetLambda");
  if (m_lambdaNeedsResetting)
  {
    Real scale = 1.0 / (m_dx*m_dx);

    // Compute it box by box, point by point
    for (DataIterator dit = m_lambda.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox&       lambdaFab = m_lambda[dit];
      const FArrayBox& aCoefFab  = (*m_aCoef)[dit];
      const FluxBox&   bCoefFab  = (*m_bCoef)[dit];
      const Box& curBox = lambdaFab.box();

      // Compute the diagonal term
      lambdaFab.copy(aCoefFab);
      lambdaFab.mult(m_alpha);

      for (int dir = 0; dir < SpaceDim; dir++)
      {
        FORT_SUMFACESNLVC(CHF_FRA(lambdaFab),
                          CHF_CONST_REAL(m_beta),
                          CHF_CONST_FRA(bCoefFab[dir]),
                          CHF_BOX(curBox),
                          CHF_CONST_INT(dir),
                          CHF_CONST_REAL(scale));
      }

      // Take its reciprocal
      lambdaFab.invert(1.0);

    }

    // Lambda is reset.
    m_lambdaNeedsResetting = false;
  }
}

// Compute the reciprocal of the diagonal entry of the operator matrix
void AMRNonLinearMultiCompOp::computeLambda()
{

  CH_TIME("AMRNonLinearMultiCompOp::computeLambda");

  //	CH_assert(!m_lambda.isDefined());
  if(m_lambda.isDefined())
  {
    // Make sure we recalculate it
    m_lambdaNeedsResetting = true;
  }

  // Define lambda
  m_lambda.define(m_aCoef->disjointBoxLayout(),m_aCoef->nComp());


  resetLambda();
}

#if 1
//
// AMRNonLinearMultiCompOp::reflux()
//   There are currently the new version (first) and the old version (second)
//   in this file.  Brian asked to preserve the old version in this way for
//   now. - TJL (12/10/2007)
//
void AMRNonLinearMultiCompOp::reflux(const LevelData<FArrayBox>&        a_phiFine,
                                     const LevelData<FArrayBox>&        a_phi,
                                     LevelData<FArrayBox>&              a_residual,
                                     AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIMERS("AMRNonLinearMultiCompOp::reflux");

  m_levfluxreg.setToZero();
  Interval interv(0,a_phi.nComp()-1);

  CH_TIMER("AMRNonLinearMultiCompOp::reflux::incrementCoarse", t2);
  CH_START(t2);

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& coarfab   = a_phi[dit];
    const FluxBox&   coarBCoef = (*m_bCoef)[dit];
    const Box&       gridBox   = a_phi.getBoxes()[dit];

    if (m_levfluxreg.hasCF(dit()))
    {
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        FArrayBox coarflux;
        Box faceBox = surroundingNodes(gridBox, idir);

        getFlux(coarflux, coarfab, coarBCoef, faceBox,idir, dit());

        Real scale = 1.0;
        m_levfluxreg.incrementCoarse(coarflux, scale,dit(),
                                     interv, interv, idir);
      }
    }
  }

  CH_STOP(t2);

  // const cast:  OK because we're changing ghost cells only
  LevelData<FArrayBox>& phiFineRef = ( LevelData<FArrayBox>&)a_phiFine;

  AMRNonLinearMultiCompOp* finerAMRPOp = static_cast<AMRNonLinearMultiCompOp*>(a_finerOp);
  QuadCFInterp& quadCFI = finerAMRPOp->m_interpWithCoarser;

  quadCFI.coarseFineInterp(phiFineRef, a_phi);
  // I'm pretty sure this is not necessary. bvs -- flux calculations use
  // outer ghost cells, but not inner ones
  // phiFineRef.exchange(a_phiFine.interval());
  int ncomps = a_phiFine.nComp();

  CH_TIMER("AMRNonLinearMultiCompOp::reflux::incrementFine", t3);
  CH_START(t3);

  DataIterator ditf = a_phiFine.dataIterator();
  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf)
  {
    const FArrayBox& phifFab   = a_phiFine[ditf];
    const FluxBox&   fineBCoef = (*(finerAMRPOp->m_bCoef))[ditf];
    const Box&       gridbox   = dblFine.get(ditf());

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      //int normalGhost = phiGhost[idir];
      SideIterator sit;
      for (sit.begin(); sit.ok(); sit.next())
      {
        if (m_levfluxreg.hasCF(ditf(), sit()))
        {
          Side::LoHiSide hiorlo = sit();
          Box fluxBox = bdryBox(gridbox,idir,hiorlo,1);

          FArrayBox fineflux(fluxBox,ncomps);
          //          getFlux(fineflux, phifFab, fineBCoef, fluxBox, idir, ditf(),
          //                  m_refToFiner);
          finerAMRPOp->getFlux(fineflux, phifFab, fineBCoef, fluxBox, idir, ditf(),
                               m_refToFiner);

          Real scale = 1.0;
          m_levfluxreg.incrementFine(fineflux, scale, ditf(),
                                     interv, interv, idir, hiorlo);
        }
      }
    }
  }

  CH_STOP(t3);

  Real scale = 1.0/m_dx;
  m_levfluxreg.reflux(a_residual, scale);
}

#else

void AMRNonLinearMultiCompOp::reflux(const LevelData<FArrayBox>&        a_phiFine,
                                     const LevelData<FArrayBox>&        a_phi,
                                     LevelData<FArrayBox>&              a_residual,
                                     AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("AMRNonLinearMultiCompOp::reflux");

  int ncomp = 1;
  ProblemDomain fineDomain = refine(m_domain, m_refToFiner);
  LevelFluxRegister levfluxreg(a_phiFine.disjointBoxLayout(),
                               a_phi.disjointBoxLayout(),
                               fineDomain,
                               m_refToFiner,
                               ncomp);

  levfluxreg.setToZero();
  Interval interv(0,a_phi.nComp()-1);

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& coarfab = a_phi[dit];
    const FluxBox& coarBCoef  = (*m_bCoef)[dit];
    const Box& gridBox = a_phi.getBoxes()[dit];

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox coarflux;
      Box faceBox = surroundingNodes(gridBox, idir);
      getFlux(coarflux, coarfab, coarBCoef , faceBox,  idir, dit());

      Real scale = 1.0;
      levfluxreg.incrementCoarse(coarflux, scale,dit(),
                                 interv,interv,idir);
    }
  }
  LevelData<FArrayBox>& p = ( LevelData<FArrayBox>&)a_phiFine;

  // has to be its own object because the finer operator
  // owns an interpolator and we have no way of getting to it
  AMRNonLinearMultiCompOp* finerAMRPOp = static_cast<AMRNonLinearMultiCompOp*> (a_finerOp);
  QuadCFInterp& quadCFI = finerAMRPOp->m_interpWithCoarser;

  quadCFI.coarseFineInterp(p, a_phi);
  // p.exchange(a_phiFine.interval()); // BVS is pretty sure this is not necesary.
  IntVect phiGhost = p.ghostVect();

  DataIterator ditf = a_phiFine.dataIterator();
  const  DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf)
  {
    const FArrayBox& phifFab = a_phiFine[ditf];
    const FluxBox& fineBCoef  = (*(finerAMRPOp->m_bCoef))[ditf];
    const Box& gridbox = dblFine.get(ditf());
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      int normalGhost = phiGhost[idir];
      SideIterator sit;
      for (sit.begin(); sit.ok(); sit.next())
      {
        Side::LoHiSide hiorlo = sit();
        Box fabbox;
        Box facebox;

        // assumption here that the stencil required
        // to compute the flux in the normal direction
        // is 2* the number of ghost cells for phi
        // (which is a reasonable assumption, and probably
        // better than just assuming you need one cell on
        // either side of the interface
        // (dfm 8-4-06)
        if (sit() == Side::Lo)
        {
          fabbox = adjCellLo(gridbox,idir, 2*normalGhost);
          fabbox.shift(idir, 1);
          facebox = bdryLo(gridbox, idir,1);
        }
        else
        {
          fabbox = adjCellHi(gridbox,idir, 2*normalGhost);
          fabbox.shift(idir, -1);
          facebox = bdryHi(gridbox, idir, 1);
        }

        // just in case we need ghost cells in the transverse direction
        // (dfm 8-4-06)
        for (int otherDir=0; otherDir<SpaceDim; ++otherDir)
        {
          if (otherDir != idir)
          {
            fabbox.grow(otherDir, phiGhost[otherDir]);
          }
        }
        CH_assert(!fabbox.isEmpty());

        FArrayBox phifab(fabbox, a_phi.nComp());
        phifab.copy(phifFab);

        FArrayBox fineflux;
        getFlux(fineflux, phifab, fineBCoef, facebox,  idir, ditf(),
                m_refToFiner);

        Real scale = 1.0;
        levfluxreg.incrementFine(fineflux, scale, ditf(),
                                 interv, interv, idir, hiorlo);
      }
    }
  }

  Real scale =  1.0/m_dx;
  levfluxreg.reflux(a_residual, scale);
}

#endif


/*
 * Computes a_phi = a_phi - lambda*(op(a_phi) - rhs)
 * where op(phi) = alpha*a*phi - beta*div(b*grad(phi))
 */
void AMRNonLinearMultiCompOp::levelGSRB(LevelData<FArrayBox>&       a_phi,
                                        const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearMultiCompOp::levelGSRB");

  // Going to need these more than once
  int nComp = a_phi.nComp();
  IntVect phiGhostVect = a_phi.ghostVect();

  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(phiGhostVect >= IntVect::Unit);
  CH_assert(nComp == a_rhs.nComp());

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();

  // This define is possibly taking an awfully long time
  // However I don't think there's a quicker way.
  // It's certainly quicker than defining FABs all the time
  LevelData<FArrayBox> derivedVar;
  {
    CH_TIME("AMRNonLinearMultiCompOp::derivedVarDefine");
    derivedVar.define(dbl, nComp,phiGhostVect);
  }

  DataIterator dit = a_phi.dataIterator();

  // Should never be homogeneous as non linear!
  bool homogeneous = false;

  Vector<int> compsList;
  compsList.push_back(-1); // do both comps in one fortran routine

  // do first red, then black passes
  for (int comps_i = 0; comps_i < compsList.size(); comps_i++)
  {
    int whichComponent = compsList[comps_i];

    for (int whichPass = 0; whichPass <= 1; whichPass++)
    {

      CH_TIMERS("AMRNonLinearMultiCompOp::levelGSRB::Compute");

      // fill in intersection of ghostcells and a_phi's boxes
      // if we're not super optimised always fill ghost cells.
      // if we are super optimised, only do this on first pass (this is dodgy but saves some time)
      if (!m_superOptimised || whichPass == 0)
      {
        {
          CH_TIME("AMRNonLinearMultiCompOp::levelGSRB::homogeneousCFInterp");
          if (homogeneous)
          {
            homogeneousCFInterp(a_phi);
          }
        }

        {
          CH_TIME("AMRNonLinearMultiCompOp::levelGSRB::exchange");
          if (s_exchangeMode == 0)
            a_phi.exchange( a_phi.interval(), m_exchangeCopier );
          else if (s_exchangeMode == 1)
            a_phi.exchangeNoOverlap(m_exchangeCopier);
          else
            MayDay::Abort("exchangeMode");

        }
      }

      // Recompute the relaxation coefficient
      resetLambda();

      {
        for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& region = dbl.get(dit());
          const FluxBox& thisBCoef  = (*m_bCoef)[dit];
          FArrayBox& thisPhi = a_phi[dit];
          FArrayBox& thisDerivedVar = derivedVar[dit];


          {
            CH_TIME("AMRNonLinearMultiCompOp::levelGSRB::BCs");
            m_bc(thisPhi, region, m_domain, m_dx, homogeneous);
//            m_bc(thisPhi, region., m_domain, m_dx, homogeneous);
          }
//        }
//
//        a_phi.exchange(a_phi.interval(), m_exchangeCopier);
//
//        for (dit.begin(); dit.ok(); ++dit)
//                {
//                  const Box& region = dbl.get(dit());
//                  const FluxBox& thisBCoef  = (*m_bCoef)[dit];
//                  FArrayBox& thisPhi = a_phi[dit];
//                  FArrayBox& thisDerivedVar = derivedVar[dit];



          // Need to do this every pass or convergence is very slow
//          if (!m_superOptimised || whichPass == 0)
//          {
            computeDiffusedVar(thisDerivedVar, thisPhi, dit(), homogeneous);
//          }

#if CH_SPACEDIM == 1
          FORT_NONLINEARSMOOTHING1D
#elif CH_SPACEDIM == 2
          FORT_NONLINEARSMOOTHING2D
#elif CH_SPACEDIM == 3
          FORT_NONLINEARSMOOTHING3D
#else
          //				This_will_not_compile!
#endif
          (CHF_FRA(thisPhi),
           CHF_FRA(thisDerivedVar),
           CHF_CONST_FRA(a_rhs[dit]),
           CHF_BOX(region),
           CHF_CONST_REAL(m_dx),
           CHF_CONST_REAL(m_alpha),
           CHF_CONST_FRA((*m_aCoef)[dit]),
           CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
           CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
           CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
           CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
           This_will_not_compile!
#endif
           CHF_CONST_FRA(m_lambda[dit]),
           CHF_CONST_INT(whichPass),
           CHF_CONST_INT(whichComponent));


        } // end loop through grids

      }
    } // end loop through red-black

  } // end loop over components

}

void AMRNonLinearMultiCompOp::levelMultiColor(LevelData<FArrayBox>&       a_phi,
                                              const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearMultiCompOp::levelMultiColor");
  MayDay::Abort("AMRNonLinearMultiCompOp::levelMultiColor - Not implemented");
}

void AMRNonLinearMultiCompOp::looseGSRB(LevelData<FArrayBox>&       a_phi,
                                        const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearMultiCompOp::looseGSRB");
#if 1
  MayDay::Abort("AMRNonLinearMultiCompOp::looseGSRB - Not implemented");
#else
  // This implementation converges at half the rate of "levelGSRB" in
  // multigrid solves
  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();

  DataIterator dit = a_phi.dataIterator();

  // fill in intersection of ghostcells and a_phi's boxes
  {
    CH_TIME("AMRNonLinearMultiCompOp::looseGSRB::homogeneousCFInterp");
    homogeneousCFInterp(a_phi);
  }

  {
    CH_TIME("AMRNonLinearMultiCompOp::looseGSRB::exchange");
    a_phi.exchange(a_phi.interval(), m_exchangeCopier);
  }

  // now step through grids...
  for (dit.begin(); dit.ok(); ++dit)
  {
    // invoke physical BC's where necessary
    {
      CH_TIME("AMRNonLinearMultiCompOp::looseGSRB::BCs");
      m_bc(a_phi[dit], dbl[dit()], m_domain, m_dx, true);
    }

    const Box& region = dbl.get(dit());
    const FluxBox& thisBCoef  = (*m_bCoef)[dit];

    int whichPass = 0;

#if CH_SPACEDIM == 1
    FORT_GSRBHELMHOLTZVC1D
#elif CH_SPACEDIM == 2
    FORT_GSRBHELMHOLTZVC2D
#elif CH_SPACEDIM == 3
    FORT_GSRBHELMHOLTZVC3D
#else
    This_will_not_compile!
#endif
    (CHF_FRA(a_phi[dit]),
        CHF_CONST_FRA(a_rhs[dit]),
        CHF_BOX(region),
        CHF_CONST_REAL(m_dx),
        CHF_CONST_REAL(m_alpha),
        CHF_CONST_FRA((*m_aCoef)[dit]),
        CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
        CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
        CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
        CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
        This_will_not_compile!
#endif
        CHF_CONST_FRA(m_lambda[dit]),
        CHF_CONST_INT(whichPass));

    whichPass = 1;

#if CH_SPACEDIM == 1
    FORT_GSRBHELMHOLTZVC1D
#elif CH_SPACEDIM == 2
    FORT_GSRBHELMHOLTZVC2D
#elif CH_SPACEDIM == 3
    FORT_GSRBHELMHOLTZVC3D
#else
    This_will_not_compile!
#endif
    (CHF_FRA(a_phi[dit]),
        CHF_CONST_FRA(a_rhs[dit]),
        CHF_BOX(region),
        CHF_CONST_REAL(m_dx),
        CHF_CONST_REAL(m_alpha),
        CHF_CONST_FRA((*m_aCoef)[dit]),
        CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
        CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
        CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
        CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
        This_will_not_compile!
#endif
        CHF_CONST_FRA(m_lambda[dit]),
        CHF_CONST_INT(whichPass));
  } // end loop through grids
#endif
}

void AMRNonLinearMultiCompOp::overlapGSRB(LevelData<FArrayBox>&       a_phi,
                                          const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearMultiCompOp::overlapGSRB");
  MayDay::Abort("AMRNonLinearMultiCompOp::overlapGSRB - Not implemented");
}

void AMRNonLinearMultiCompOp::levelGSRBLazy(LevelData<FArrayBox>&       a_phi,
                                            const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearMultiCompOp::levelGSRBLazy");
  MayDay::Abort("AMRNonLinearMultiCompOp::levelGSRBLazy - Not implemented");
}

void AMRNonLinearMultiCompOp::levelJacobi(LevelData<FArrayBox>&       a_phi,
                                          const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearMultiCompOp::levelJacobi");

  resetLambda();

  LevelData<FArrayBox> resid;
  create(resid, a_rhs);

  // Get the residual
  residual(resid,a_phi,a_rhs,false);

  // Multiply by the weights
  DataIterator dit = m_lambda.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    resid[dit].mult(m_lambda[dit]);
  }

  // Do the Jacobi relaxation
  Real weight = 0.5;

  incr(a_phi, resid, weight);
}

void AMRNonLinearMultiCompOp::getFlux(FArrayBox&       a_flux,
                                      const FArrayBox& a_data,
                                      const FluxBox&   a_bCoef,
                                      const Box&       a_facebox,
                                      int              a_dir,
                                      const DataIndex&  a_dit,
                                      int              a_ref) const
{
  CH_TIME("AMRNonLinearMultiCompOp::getFlux");

  CH_assert(a_dir >= 0);
  CH_assert(a_dir <  SpaceDim);
  CH_assert(!a_data.box().isEmpty());
  CH_assert(!a_facebox.isEmpty());

  // probably the simplest way to test centering
  // a_box needs to be face-centered in the a_dir
  Box faceTestBox(IntVect::Zero, IntVect::Unit);
  faceTestBox.surroundingNodes(a_dir);
  CH_assert(a_facebox.type() == faceTestBox.type());

  const FArrayBox& bCoefDir = a_bCoef[a_dir];

  // reality check for bCoef
  CH_assert(bCoefDir.box().contains(a_facebox));

  a_flux.resize(a_facebox, a_data.nComp());
  BoxIterator bit(a_facebox);

  Real scale = m_beta * a_ref / m_dx;

  // a_data is Enthalpy, Bulk Concentration etc.
  // need to convert this to temperature, liquid concentration to calculate diffusive flux
  FArrayBox derivedVar(a_data.box(), a_data.nComp());

  AMRNonLinearMultiCompOp* op = const_cast<AMRNonLinearMultiCompOp*> (this);
  op->computeDiffusedVar(derivedVar, a_data, a_dit);

  // const FArrayBox& diffusedVar = a_data;
  FArrayBox& diffusedVar = derivedVar;// should be this, but it isn't working?

  for ( bit.begin(); bit.ok(); bit.next())
  {
    IntVect iv = bit();
    IntVect shiftiv = BASISV(a_dir);
    IntVect ivlo = iv - shiftiv;
    IntVect ivhi = iv;

    CH_assert(a_data.box().contains(ivlo));
    CH_assert(a_data.box().contains(ivhi));

    for (int ivar = 0; ivar < a_data.nComp(); ivar++)
    {
      Real phihi = diffusedVar(ivhi,ivar);
      Real philo = diffusedVar(ivlo,ivar);
      Real gradphi = (phihi - philo ) * scale;

      a_flux(iv,ivar) = -bCoefDir(iv, ivar) * gradphi;
    }
  }
}

//-----------------------------------------------------------------------
void
AMRNonLinearMultiCompOp::
setTime(Real a_time)
{
  // Jot down the time.
  m_time = a_time;

  // Interpolate the b coefficient data if necessary / possible. If
  // the B coefficient depends upon the solution, the operator is nonlinear
  // and the integrator must decide how to treat it.
  if (!m_bCoefInterpolator.isNull() &&
      !m_bCoefInterpolator->dependsUponSolution())
    m_bCoefInterpolator->interpolate(*m_bCoef, a_time);

  if (!m_aCoefInterpolator.isNull() &&
      !m_aCoefInterpolator->dependsUponSolution())
    m_aCoefInterpolator->interpolate(*m_aCoef, a_time);

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;

  // Set the time on the boundary holder.
  m_bc.setTime(a_time);

  // Notify our observers that the time has been set.
  notifyObserversOfChange();
}
//-----------------------------------------------------------------------

// Factory
AMRNonLinearMultiCompOpFactory::AMRNonLinearMultiCompOpFactory() : m_FAS(true), m_numComp(2)
{
  setDefaultValues();
}


void AMRNonLinearMultiCompOpFactory::define(const ProblemDomain& a_coarseDomain,
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
                                            bool a_apply_bcs_to_diagnostic_var)
{
  CH_TIME("AMRNonLinearMultiCompOpFactory::define");

  setDefaultValues();

  m_boxes = a_grids;

  m_refRatios = a_refRatios;

  m_bc = a_bc;
  m_diffusedVarBC = a_derivedVarBC;
  m_apply_bcs_to_diagnostic_var = a_apply_bcs_to_diagnostic_var;
  m_porosityEdgeBC = a_porosityEdgeBC;

  m_params = a_params;

  m_dx.resize(a_grids.size());
  m_dx[0] = a_coarsedx;

  m_domains.resize(a_grids.size());
  m_domains[0] = a_coarseDomain;

  m_exchangeCopiers.resize(a_grids.size());
  m_exchangeCopiers[0].exchangeDefine(a_grids[0], IntVect::Unit);
  m_exchangeCopiers[0].trimEdges(a_grids[0], IntVect::Unit);

  m_cfregion.resize(a_grids.size());
  m_cfregion[0].define(a_grids[0], m_domains[0]);

  for (int i = 1; i < a_grids.size(); i++)
  {
    m_dx[i] = m_dx[i-1] / m_refRatios[i-1];

    m_domains[i] = m_domains[i-1];
    m_domains[i].refine(m_refRatios[i-1]);

    m_exchangeCopiers[i].exchangeDefine(a_grids[i], IntVect::Unit);
    m_exchangeCopiers[i].trimEdges(a_grids[i], IntVect::Unit);

    m_cfregion[i].define(a_grids[i], m_domains[i]);
  }

  m_alpha = a_alpha;
  m_aCoef = a_aCoef;

  m_beta  = a_beta;
  m_bCoef = a_bCoef;

  m_enthalpySolidus = a_enthalpySolidus;
  m_enthalpyLiquidus = a_enthalpyLiquidus;
  m_enthalpyEutectic = a_enthalpyEutectic;

  //  m_FAS = false;
  m_FAS = true;
  ParmParse ppMG("amrmultigrid");
  ppMG.query("MGtype", m_FAS); // if MGtype = 1, do nonlinear smoothing

  m_relaxMode = a_relaxMode;

  m_numComp = m_bCoef[0]->nComp();

  m_superOptimised = false;

}

void AMRNonLinearMultiCompOpFactory::setSuperOptimised(bool a_val)
{
  m_superOptimised = a_val;
}

void AMRNonLinearMultiCompOpFactory::setBC(BCHolder& a_bc)
{
  m_bc = a_bc;
}


//-----------------------------------------------------------------------

MGLevelOp<LevelData<FArrayBox> >* AMRNonLinearMultiCompOpFactory::MGnewOp(const ProblemDomain& a_indexSpace,
                                                                          int                  a_depth,
                                                                          bool                 a_homoOnly)
{
  CH_TIME("AMRNonLinearMultiCompOpFactory::MGnewOp");

  Real dxCrse = -1.0;

  int ref;

  for (ref = 0; ref < m_domains.size(); ref++)
  {
    if (a_indexSpace.domainBox() == m_domains[ref].domainBox())
    {
      break;
    }
  }

  CH_assert(ref !=  m_domains.size()); // didn't find domain

  if (ref > 0)
  {
    dxCrse = m_dx[ref-1];
  }

  ProblemDomain domain(m_domains[ref]);
  Real dx = m_dx[ref];
  int coarsening = 1;

  for (int i = 0; i < a_depth; i++)
  {
    coarsening *= 2;
    domain.coarsen(2);
  }

  if (coarsening > 1 && !m_boxes[ref].coarsenable(coarsening*AMRNonLinearMultiCompOp::s_maxCoarse))
  {
    return nullptr;
  }

  dx *= coarsening;

  DisjointBoxLayout layout, layoutCrse;
  coarsen_dbl(layout, m_boxes[ref], coarsening);

  // If there are coarser levels
  // ref = 0 means no amr refinement (starting from amr level 0)
  if (ref > 0)
  {
    // This assumes that the ref ratio from this level n to the coarser level n-1 is the same
    // as that from level n-1 to n-2
    if (m_boxes[ref-1].coarsenable(coarsening))
    {
      coarsen_dbl(layoutCrse, m_boxes[ref-1], coarsening);
    }
    else
    {
      return nullptr;
    }
  }

  Copier ex = m_exchangeCopiers[ref];
  CFRegion cfregion = m_cfregion[ref];

  if (coarsening > 1)
  {
    ex.coarsen(coarsening);
    cfregion.coarsen(coarsening);
  }


  AMRNonLinearMultiCompOp* newOp = new AMRNonLinearMultiCompOp;

  // Can we just do AMR define with a coarse level?
  int refRatio = 2; // default MG refinement
  if (ref-1 < m_refRatios.size())
  {
    refRatio = m_refRatios[ref-1]; // AMR refinement if defined, else default back to MG
  }


  // ref > 0 means we've done AMR coarsening i.e. there are coarser amr levels
  if (ref > 0)
  {
    newOp->define(layout, layoutCrse, dx,
                  refRatio,
                  domain, m_bc,
                  ex, cfregion, m_numComp);
  }
  else
  {
    // no coarse level define
    newOp->define(layout, dx, domain, m_bc, ex, cfregion);

  }


  // Now need to generate correctly coarsened data on this level


  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;

  newOp->m_params = m_params;

  newOp->m_diffusedVarBC = m_diffusedVarBC;

  newOp->m_apply_bcs_to_diagnostic_var = m_apply_bcs_to_diagnostic_var;

  newOp->m_superOptimised = m_superOptimised;

  newOp->s_relaxMode = m_relaxMode; // 4 is for jacobi, 1 is for GSRB

  newOp->m_FAS = m_FAS;

  RefCountedPtr<CoefficientInterpolatorLinear> aCoefInterpolator;
  RefCountedPtr<CoefficientInterpolatorLinearFace> bCoefInterpolator;

  aCoefInterpolator = RefCountedPtr<CoefficientInterpolatorLinear>
  (new CoefficientInterpolatorLinear);

  bCoefInterpolator = RefCountedPtr<CoefficientInterpolatorLinearFace>
  (new CoefficientInterpolatorLinearFace);

  newOp->setACoefInterpolator(aCoefInterpolator);
  newOp->setBCoefInterpolator(bCoefInterpolator);

  if (a_depth == 0)
  {
    // don't need to coarsen anything for this
    newOp->m_aCoef = m_aCoef[ref];
    newOp->m_bCoef = m_bCoef[ref];

    newOp->m_enthalpySolidus = m_enthalpySolidus[ref];
    newOp->m_enthalpyLiquidus = m_enthalpyLiquidus[ref];
    newOp->m_enthalpyEutectic = m_enthalpyEutectic[ref];



    //    if ()

  }
  else
  {
    // need to coarsen coefficients
    RefCountedPtr<LevelData<FArrayBox> > aCoef( new LevelData<FArrayBox> );
    RefCountedPtr<LevelData<FluxBox> > bCoef( new LevelData<FluxBox> );
    RefCountedPtr<LevelData<FArrayBox> > enthalpy( new LevelData<FArrayBox> );
    RefCountedPtr<LevelData<FArrayBox> > enthalpySolidus( new LevelData<FArrayBox> );
    RefCountedPtr<LevelData<FArrayBox> > enthalpyLiquidus( new LevelData<FArrayBox> );
    RefCountedPtr<LevelData<FArrayBox> > enthalpyEutectic( new LevelData<FArrayBox> );


    aCoef->define(layout, m_aCoef[ref]->nComp(), m_aCoef[ref]->ghostVect());

    enthalpySolidus->define(layout, m_enthalpySolidus[ref]->nComp(), m_enthalpySolidus[ref]->ghostVect());
    enthalpyLiquidus->define(layout, m_enthalpyLiquidus[ref]->nComp(), m_enthalpyLiquidus[ref]->ghostVect());
    enthalpyEutectic->define(layout, m_enthalpyEutectic[ref]->nComp(), m_enthalpyEutectic[ref]->ghostVect());

    bCoef->define(layout, m_bCoef[ref]->nComp(), m_bCoef[ref]->ghostVect());


    // average coefficients to coarser level
    // for now, do this with a CoarseAverage --
    // may want to switch to harmonic averaging at some point
    CoarseAverage averager(m_aCoef[ref]->getBoxes(),
                           layout, aCoef->nComp(), coarsening);

    /// Different number of components for these
    CoarseAverage averagerBoundingEnergies(m_aCoef[ref]->getBoxes(),
                               layout, enthalpySolidus->nComp(), coarsening);

    CoarseAverageFace faceAverager(m_bCoef[ref]->getBoxes(),
                                   bCoef->nComp(), coarsening);

    if (m_coefficient_average_type == CoarseAverage::arithmetic)
    {
      averager.averageToCoarse(*aCoef, *(m_aCoef[ref]));
      faceAverager.averageToCoarse(*bCoef, *(m_bCoef[ref]));

      averagerBoundingEnergies.averageToCoarse(*enthalpySolidus, *(m_enthalpySolidus[ref]));
      averagerBoundingEnergies.averageToCoarse(*enthalpyLiquidus, *(m_enthalpyLiquidus[ref]));
      averagerBoundingEnergies.averageToCoarse(*enthalpyEutectic, *(m_enthalpyEutectic[ref]));



    }
    else if (m_coefficient_average_type == CoarseAverage::harmonic)
    {
      averager.averageToCoarseHarmonic(*aCoef, *(m_aCoef[ref]));
      faceAverager.averageToCoarseHarmonic(*bCoef, *(m_bCoef[ref]));

      averagerBoundingEnergies.averageToCoarseHarmonic(*enthalpySolidus, *(m_enthalpySolidus[ref]));
      averagerBoundingEnergies.averageToCoarseHarmonic(*enthalpyLiquidus, *(m_enthalpyLiquidus[ref]));
      averagerBoundingEnergies.averageToCoarseHarmonic(*enthalpyEutectic, *(m_enthalpyEutectic[ref]));

    }
    else
    {
      MayDay::Abort("AMRNonLinearMultiCompOpFactory::MGNewOp -- bad averagetype");
    }

    newOp->m_aCoef = aCoef;
    newOp->m_bCoef = bCoef;

    newOp->m_enthalpySolidus = enthalpySolidus;
    newOp->m_enthalpyLiquidus = enthalpyLiquidus;
    newOp->m_enthalpyEutectic = enthalpyEutectic;


  }

  newOp->computeLambda();


  newOp->m_dxCrse = dxCrse;



  return (MGLevelOp<LevelData<FArrayBox> >*)newOp;



}

AMRLevelOp<LevelData<FArrayBox> >* AMRNonLinearMultiCompOpFactory::AMRnewOp(const ProblemDomain& a_indexSpace)
{
  CH_TIME("AMRNonLinearMultiCompOpFactory::AMRnewOp");

  //	MayDay::Error("Haven't properly edited AMRNonLinearMultiCompOpFactory::AMRnewOp yet");

  AMRNonLinearMultiCompOp* newOp = new AMRNonLinearMultiCompOp;
  Real dxCrse = -1.0;

  int ref;

  for (ref = 0; ref< m_domains.size(); ref++)
  {
    if (a_indexSpace.domainBox() == m_domains[ref].domainBox())
    {
      break;
    }
  }

  if (ref == 0)
  {
    // coarsest AMR level
    if (m_domains.size() == 1)
    {
      // no finer level
      newOp->define(m_boxes[0], m_dx[0],
                    a_indexSpace, m_bc,
                    m_exchangeCopiers[0], m_cfregion[0]);
    }
    else
    {
      // finer level exists but no coarser
      int dummyRat = 1;  // argument so compiler can find right function
      int refToFiner = m_refRatios[0]; // actual refinement ratio
      newOp->define(m_boxes[0],  m_boxes[1], m_dx[0],
                    dummyRat, refToFiner,
                    a_indexSpace, m_bc,
                    m_exchangeCopiers[0], m_cfregion[0]);
    }
  }
  else if (ref ==  m_domains.size()-1)
  {
    dxCrse = m_dx[ref-1];

    // finest AMR level
    newOp->define(m_boxes[ref], m_boxes[ref-1], m_dx[ref],
                  m_refRatios[ref-1],
                  a_indexSpace, m_bc,
                  m_exchangeCopiers[ref], m_cfregion[ref], m_numComp);
  }
  else if ( ref == m_domains.size())
  {
    MayDay::Abort("Did not find a domain to match AMRnewOp(const ProblemDomain& a_indexSpace)");

  }
  else
  {
    dxCrse = m_dx[ref-1];

    // intermediate AMR level, full define
    newOp->define(m_boxes[ref], m_boxes[ref+1], m_boxes[ref-1], m_dx[ref],
                  m_refRatios[ref-1], m_refRatios[ref],
                  a_indexSpace, m_bc,
                  m_exchangeCopiers[ref], m_cfregion[ref], m_numComp);
  }

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;

  newOp->m_aCoef = m_aCoef[ref];
  newOp->m_bCoef = m_bCoef[ref];

  newOp->m_enthalpySolidus = m_enthalpySolidus[ref];
  newOp->m_enthalpyLiquidus = m_enthalpyLiquidus[ref];
  newOp->m_enthalpyEutectic = m_enthalpyEutectic[ref];

  newOp->m_params = m_params;
  //  newOp->m_calcEnthalpyVar = m_calcEnthalpyVar;
  newOp->m_diffusedVarBC = m_diffusedVarBC;
  newOp->m_apply_bcs_to_diagnostic_var = m_apply_bcs_to_diagnostic_var;
  //  newOp->m_computeEnthalpyVars =m_computeEnthalpyVars;
  //  (newOp->m_computeEnthalpyVars).define(m_computeEnthalpyVars.m_derivedVar, a_indexSpace, m_dx[ref], m_computeEnthalpyVars.m_derivedVarBC);
  //  newOp->m_porosityEdgeBC = m_porosityEdgeBC;

  newOp->m_superOptimised = m_superOptimised;

  newOp->m_FAS = m_FAS;

  newOp->computeLambda();

  newOp->m_dxCrse = dxCrse;

  newOp->s_relaxMode = m_relaxMode; // 4 is for jacobi, 1 is for GSRB

  return (AMRLevelOp<LevelData<FArrayBox> >*)newOp;
}

int AMRNonLinearMultiCompOpFactory::refToFiner(const ProblemDomain& a_domain) const
{
  int retval = -1;
  bool found = false;

  for (int ilev = 0; ilev < m_domains.size(); ilev++)
  {
    if (m_domains[ilev].domainBox() == a_domain.domainBox())
    {
      retval = m_refRatios[ilev];
      found = true;
    }
  }

  if (!found)
  {
    MayDay::Abort("Domain not found in AMR hierarchy");
  }

  return retval;
}

//-----------------------------------------------------------------------
void AMRNonLinearMultiCompOpFactory::setDefaultValues()
{
  // Default to Laplacian operator
  m_alpha = 0.0;
  m_beta = -1.0;

  m_coefficient_average_type = CoarseAverage::arithmetic;

  m_superOptimised = false;

  m_numComp = 2;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AMRNonLinearMultiCompOp::
finerOperatorChanged(const MGLevelOp<LevelData<FArrayBox> >& a_operator,
                     int a_coarseningFactor)
{
  const AMRNonLinearMultiCompOp& op =
      dynamic_cast<const AMRNonLinearMultiCompOp&>(a_operator);

  MayDay::Error("Haven't properly edited AMRNonLinearMultiCompOp::finerOperatorChanged yet");

  // Perform multigrid coarsening on the operator data.
  LevelData<FArrayBox>& acoefCoar = *m_aCoef;
  const LevelData<FArrayBox>& acoefFine = *(op.m_aCoef);
  LevelData<FluxBox>& bcoefCoar = *m_bCoef;
  const LevelData<FluxBox>& bcoefFine = *(op.m_bCoef);
  if (a_coarseningFactor != 1)
  {
    CoarseAverage cellAverage(acoefFine.disjointBoxLayout(),
                              acoefCoar.disjointBoxLayout(),
                              1, a_coarseningFactor);
    for (DataIterator dit = acoefCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
      acoefCoar[dit()].setVal(0.);
    cellAverage.averageToCoarse(acoefCoar, acoefFine);

    CoarseAverageFace faceAverage(bcoefFine.disjointBoxLayout(),
                                  1, a_coarseningFactor);
    for (DataIterator dit = bcoefCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
      bcoefCoar[dit()].setVal(0.);
    faceAverage.averageToCoarse(bcoefCoar, bcoefFine);
  }

  // Handle inter-box ghost cells.
  acoefCoar.exchange();
  bcoefCoar.exchange();

  // Mark the relaxation coefficient dirty.
  m_lambdaNeedsResetting = true;

  // Notify any observers of this change.
  notifyObserversOfChange();
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

