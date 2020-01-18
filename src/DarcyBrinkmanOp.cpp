
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

#include "DarcyBrinkmanOp.H"
#include "DarcyBrinkmanOpF_F.H"
#include "DebugOut.H"

#include "NamespaceHeader.H"

void DarcyBrinkmanOp::residualI(LevelData<FArrayBox>&       a_lhs,
                                const LevelData<FArrayBox>& a_phi,
                                const LevelData<FArrayBox>& a_rhs,
                                bool                        a_homogeneous)
{
  CH_TIME("DarcyBrinkmanOp::residualI");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  Real dx = m_dx;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  {
    CH_TIME("DarcyBrinkmanOp::residualIBC");

    for (dit.begin(); dit.ok(); ++dit)
    {
      m_bc(phi[dit], dbl[dit()],m_domain, dx, a_homogeneous);
    }
  }

  phi.exchange(phi.interval(), m_exchangeCopier);

  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& region = dbl[dit()];
    const FluxBox& thisBCoef = (*m_bCoef)[dit];

#if CH_SPACEDIM == 1
    FORT_DBCOMPUTERES1D
#elif CH_SPACEDIM == 2
    FORT_DBCOMPUTERES2D
#elif CH_SPACEDIM == 3
    FORT_DBCOMPUTERES3D
#else
    //		This_will_not_compile!
#endif
    (CHF_FRA(a_lhs[dit]),
     CHF_CONST_FRA(phi[dit]),
     CHF_CONST_FRA(a_rhs[dit]),
     CHF_CONST_REAL(m_alpha),
     CHF_CONST_FRA((*m_aCoef)[dit]),
     CHF_CONST_FRA((*m_cCoef)[dit]),
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
void DarcyBrinkmanOp::preCond(LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("DarcyBrinkmanOp::preCond");

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

void DarcyBrinkmanOp::applyOpI(LevelData<FArrayBox>&      a_lhs,
                               const LevelData<FArrayBox>& a_phi,
                               bool                        a_homogeneous )
{
  CH_TIME("DarcyBrinkmanOp::applyOpI");
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  Real dx = m_dx;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
  {
    m_bc(phi[dit], dbl[dit()],m_domain, dx, a_homogeneous);
  }

  applyOpNoBoundary(a_lhs, a_phi);
}

void DarcyBrinkmanOp::applyOpNoBoundary(LevelData<FArrayBox>&      a_lhs,
                                        const LevelData<FArrayBox>& a_phi)
{
  CH_TIME("DarcyBrinkmanOp::applyOpNoBoundary");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();

  phi.exchange(phi.interval(), m_exchangeCopier);

  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& region = dbl[dit()];
    const FluxBox& thisBCoef = (*m_bCoef)[dit];

#if CH_SPACEDIM == 1
    FORT_DBCOMPUTEOP1D
#elif CH_SPACEDIM == 2
    FORT_DBCOMPUTEOP2D
#elif CH_SPACEDIM == 3
    FORT_DBCOMPUTEOP3D
#else
    //		This_will_not_compile!
#endif
    (CHF_FRA(a_lhs[dit]),
     CHF_CONST_FRA(phi[dit]),
     CHF_CONST_REAL(m_alpha),
     CHF_CONST_FRA((*m_aCoef)[dit]),
     CHF_CONST_FRA((*m_cCoef)[dit]),
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

void DarcyBrinkmanOp::restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                       LevelData<FArrayBox>&       a_phiFine,
                                       const LevelData<FArrayBox>& a_rhsFine)
{
  CH_TIME("DarcyBrinkmanOp::restrictResidual");

  homogeneousCFInterp(a_phiFine);
  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
  {
    FArrayBox& phi = a_phiFine[dit];
    m_bc(phi, dblFine[dit()], m_domain, m_dx, true);
  }

  a_phiFine.exchange(a_phiFine.interval(), m_exchangeCopier);

  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
  {
    FArrayBox&       phi = a_phiFine[dit];
    const FArrayBox& rhs = a_rhsFine[dit];
    FArrayBox&       res = a_resCoarse[dit];

    const FArrayBox& thisACoef = (*m_aCoef)[dit];
    const FArrayBox& thisCCoef = (*m_cCoef)[dit];
    const FluxBox&   thisBCoef = (*m_bCoef)[dit];

    Box region = dblFine.get(dit());
    const IntVect& iv = region.smallEnd();
    IntVect civ = coarsen(iv, 2);

    res.setVal(0.0);

#if CH_SPACEDIM == 1
    FORT_RESTRICTRESDB1D
#elif CH_SPACEDIM == 2
    FORT_RESTRICTRESDB2D
#elif CH_SPACEDIM == 3
    FORT_RESTRICTRESDB3D
#else
    //		This_will_not_compile!
#endif
    (CHF_FRA_SHIFT(res, civ),
     CHF_CONST_FRA_SHIFT(phi, iv),
     CHF_CONST_FRA_SHIFT(rhs, iv),
     CHF_CONST_REAL(m_alpha),
     CHF_CONST_FRA_SHIFT(thisACoef, iv),
     CHF_CONST_FRA_SHIFT(thisCCoef, iv),
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

void DarcyBrinkmanOp::setAlphaAndBeta(const Real& a_alpha,
                                      const Real& a_beta)
{
  m_alpha = a_alpha;
  m_beta  = a_beta;

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;
}

void DarcyBrinkmanOp::setCoefs(const RefCountedPtr<LevelData<FArrayBox> >& a_aCoef,
                               const RefCountedPtr<LevelData<FluxBox  > >& a_bCoef,
                               const Real&                                 a_alpha,
                               const Real&                                 a_beta,
                               const RefCountedPtr<LevelData<FArrayBox> >& a_cCoef)
{
  m_alpha = a_alpha;
  m_beta  = a_beta;

  m_aCoef = a_aCoef;
  m_bCoef = a_bCoef;
  m_cCoef = a_cCoef;

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;
}

void DarcyBrinkmanOp::resetLambda()
{
  if (m_lambdaNeedsResetting)
  {
    Real scale = 1.0 / (m_dx*m_dx);

    // Compute it box by box, point by point
    for (DataIterator dit = m_lambda.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox&       lambdaFab = m_lambda[dit];
      const FArrayBox& aCoefFab  = (*m_aCoef)[dit];
      const FArrayBox& cCoefFab  = (*m_cCoef)[dit];
      const FluxBox&   bCoefFab  = (*m_bCoef)[dit];
      const Box& curBox = lambdaFab.box();

      // Compute the diagonal term
      lambdaFab.copy(aCoefFab);
      lambdaFab.mult(m_alpha);

      lambdaFab.plus(cCoefFab, -m_beta);

      for (int dir = 0; dir < SpaceDim; dir++)
      {
        FORT_SUMFACES3(CHF_FRA(lambdaFab),
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
void DarcyBrinkmanOp::computeLambda()
{
  CH_TIME("DarcyBrinkmanOp::computeLambda");

  CH_assert(!m_lambda.isDefined());

  // Define lambda
  m_lambda.define(m_aCoef->disjointBoxLayout(),m_aCoef->nComp());
  resetLambda();
}

#if 1
//
// DarcyBrinkmanOp::reflux()
//   There are currently the new version (first) and the old version (second)
//   in this file.  Brian asked to preserve the old version in this way for
//   now. - TJL (12/10/2007)
//
void DarcyBrinkmanOp::reflux(const LevelData<FArrayBox>&        a_phiFine,
                             const LevelData<FArrayBox>&        a_phi,
                             LevelData<FArrayBox>&              a_residual,
                             AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIMERS("DarcyBrinkmanOp::reflux");

  m_levfluxreg.setToZero();
  Interval interv(0,a_phi.nComp()-1);

  CH_TIMER("DarcyBrinkmanOp::reflux::incrementCoarse", t2);
  CH_START(t2);

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& coarfab   = a_phi[dit];
    const FArrayBox& coarCCoef = (*m_cCoef)[dit];
    const FluxBox&   coarBCoef = (*m_bCoef)[dit];
    const Box&       gridBox   = a_phi.getBoxes()[dit];

    if (m_levfluxreg.hasCF(dit()))
    {
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        FArrayBox coarflux;
        Box faceBox = surroundingNodes(gridBox, idir);

        getFlux(coarflux, coarfab, coarCCoef, coarBCoef, faceBox, idir);

        Real scale = 1.0;
        m_levfluxreg.incrementCoarse(coarflux, scale,dit(),
                                     interv, interv, idir);
      }
    }
  }

  CH_STOP(t2);

  // const cast:  OK because we're changing ghost cells only
  LevelData<FArrayBox>& phiFineRef = ( LevelData<FArrayBox>&)a_phiFine;

  DarcyBrinkmanOp* finerAMRPOp = static_cast<DarcyBrinkmanOp*>(a_finerOp);
  QuadCFInterp& quadCFI = finerAMRPOp->m_interpWithCoarser;

  quadCFI.coarseFineInterp(phiFineRef, a_phi);
  // I'm pretty sure this is not necessary. bvs -- flux calculations use
  // outer ghost cells, but not inner ones
  // phiFineRef.exchange(a_phiFine.interval());
  int ncomps = a_phiFine.nComp();

  CH_TIMER("DarcyBrinkmanOp::reflux::incrementFine", t3);
  CH_START(t3);

  DataIterator ditf = a_phiFine.dataIterator();
  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf)
  {
    const FArrayBox& phifFab   = a_phiFine[ditf];
    const FluxBox&   fineBCoef = (*(finerAMRPOp->m_bCoef))[ditf];
    const FArrayBox&   fineCCoef = (*(finerAMRPOp->m_cCoef))[ditf];
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
          getFlux(fineflux, phifFab, fineCCoef, fineBCoef, fluxBox, idir,
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

void DarcyBrinkmanOp::reflux(const LevelData<FArrayBox>&        a_phiFine,
                             const LevelData<FArrayBox>&        a_phi,
                             LevelData<FArrayBox>&              a_residual,
                             AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("DarcyBrinkmanOp::reflux");

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
      getFlux(coarflux, coarfab, coarBCoef , faceBox, idir);

      Real scale = 1.0;
      levfluxreg.incrementCoarse(coarflux, scale,dit(),
                                 interv,interv,idir);
    }
  }
  LevelData<FArrayBox>& p = ( LevelData<FArrayBox>&)a_phiFine;

  // has to be its own object because the finer operator
  // owns an interpolator and we have no way of getting to it
  DarcyBrinkmanOp* finerAMRPOp = (DarcyBrinkmanOp*) a_finerOp;
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
        getFlux(fineflux, phifab, fineBCoef, facebox, idir,
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

void DarcyBrinkmanOp::levelGSRB(LevelData<FArrayBox>&       a_phi,
                                const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("DarcyBrinkmanOp::levelGSRB");

  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();

  DataIterator dit = a_phi.dataIterator();

  // do first red, then black passes
  for (int whichPass = 0; whichPass <= 1; whichPass++)
  {
    CH_TIMERS("DarcyBrinkmanOp::levelGSRB::Compute");

    // fill in intersection of ghostcells and a_phi's boxes
    {
      CH_TIME("DarcyBrinkmanOp::levelGSRB::homogeneousCFInterp");
      homogeneousCFInterp(a_phi);
    }

    {
      CH_TIME("DarcyBrinkmanOp::levelGSRB::exchange");
      a_phi.exchange(a_phi.interval(), m_exchangeCopier);
    }

    {
      CH_TIME("DarcyBrinkmanOp::levelGSRB::BCs");
      // now step through grids...
      for (dit.begin(); dit.ok(); ++dit)
      {
        // invoke physical BC's where necessary
        m_bc(a_phi[dit], dbl[dit()], m_domain, m_dx, true);
      }
    }

    for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& region = dbl.get(dit());
      const FluxBox& thisBCoef  = (*m_bCoef)[dit];

#if CH_SPACEDIM == 1
      FORT_GSRBHELMHOLTZDB1D
#elif CH_SPACEDIM == 2
      FORT_GSRBHELMHOLTZDB2D
#elif CH_SPACEDIM == 3
      FORT_GSRBHELMHOLTZDB3D
#else
      //			This_will_not_compile!
#endif
      (CHF_FRA(a_phi[dit]),
       CHF_CONST_FRA(a_rhs[dit]),
       CHF_BOX(region),
       CHF_CONST_REAL(m_dx),
       CHF_CONST_REAL(m_alpha),
       CHF_CONST_FRA((*m_aCoef)[dit]),
       CHF_CONST_FRA((*m_cCoef)[dit]),
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
  } // end loop through red-black
}

void DarcyBrinkmanOp::levelMultiColor(LevelData<FArrayBox>&       a_phi,
                                      const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("DarcyBrinkmanOp::levelMultiColor");
  MayDay::Abort("DarcyBrinkmanOp::levelMultiColor - Not implemented");
}

void DarcyBrinkmanOp::looseGSRB(LevelData<FArrayBox>&       a_phi,
                                const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("DarcyBrinkmanOp::looseGSRB");
#if 1
  MayDay::Abort("DarcyBrinkmanOp::looseGSRB - Not implemented");
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
    CH_TIME("DarcyBrinkmanOp::looseGSRB::homogeneousCFInterp");
    homogeneousCFInterp(a_phi);
  }

  {
    CH_TIME("DarcyBrinkmanOp::looseGSRB::exchange");
    a_phi.exchange(a_phi.interval(), m_exchangeCopier);
  }

  // now step through grids...
  for (dit.begin(); dit.ok(); ++dit)
  {
    // invoke physical BC's where necessary
    {
      CH_TIME("DarcyBrinkmanOp::looseGSRB::BCs");
      m_bc(a_phi[dit], dbl[dit()], m_domain, m_dx, true);
    }

    const Box& region = dbl.get(dit());
    const FluxBox& thisBCoef  = (*m_bCoef)[dit];

    int whichPass = 0;

#if CH_SPACEDIM == 1
    FORT_GSRBHELMHOLTZDB1D
#elif CH_SPACEDIM == 2
    FORT_GSRBHELMHOLTZDB2D
#elif CH_SPACEDIM == 3
    FORT_GSRBHELMHOLTZDB3D
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
    FORT_GSRBHELMHOLTZDB1D
#elif CH_SPACEDIM == 2
    FORT_GSRBHELMHOLTZDB2D
#elif CH_SPACEDIM == 3
    FORT_GSRBHELMHOLTZDB3D
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

void DarcyBrinkmanOp::overlapGSRB(LevelData<FArrayBox>&       a_phi,
                                  const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("DarcyBrinkmanOp::overlapGSRB");
  MayDay::Abort("DarcyBrinkmanOp::overlapGSRB - Not implemented");
}

void DarcyBrinkmanOp::levelGSRBLazy(LevelData<FArrayBox>&       a_phi,
                                    const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("DarcyBrinkmanOp::levelGSRBLazy");
  MayDay::Abort("DarcyBrinkmanOp::levelGSRBLazy - Not implemented");
}

void DarcyBrinkmanOp::levelJacobi(LevelData<FArrayBox>&       a_phi,
                                  const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("DarcyBrinkmanOp::levelJacobi");

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  LevelData<FArrayBox> resid;
  create(resid, a_rhs);

  // Get the residual
  residual(resid,a_phi,a_rhs,true);

  // Multiply by the weights
  DataIterator dit = m_lambda.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    resid[dit].mult(m_lambda[dit]);
  }

  // Do the Jacobi relaxation
  incr(a_phi, resid, 0.5);
}

void DarcyBrinkmanOp::getFlux(FArrayBox&       a_flux,
                              const FArrayBox& a_data,
                              const FArrayBox& a_cCoef,
                              const FluxBox&   a_bCoef,
                              const Box&       a_facebox,
                              int              a_dir,
                              int              a_ref) const
{
  CH_TIME("DarcyBrinkmanOp::getFlux");

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

  Real scale = m_beta * a_ref;
  Real grad_scale = scale / m_dx;

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
      Real phihi = a_data(ivhi,ivar);
      Real philo = a_data(ivlo,ivar);
      Real gradphi = (phihi - philo ) * grad_scale;

//      Real phi = a_data(iv, ivar) * scale;

      // Note this doesn't depend on cCoef
      a_flux(iv,ivar) = -bCoefDir(iv, ivar) * gradphi;

    }
  }
}

//-----------------------------------------------------------------------
void
DarcyBrinkmanOp::
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

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;

  // Set the time on the boundary holder.
  m_bc.setTime(a_time);

  // Notify our observers that the time has been set.
  // FIXME: Must implement response of multigrid operators!
  notifyObserversOfChange();
}
//-----------------------------------------------------------------------

// Factory
DarcyBrinkmanOpFactory::DarcyBrinkmanOpFactory()
{
  setDefaultValues();
}

//-----------------------------------------------------------------------
//  AMR Factory define function
void DarcyBrinkmanOpFactory::define(const ProblemDomain&                           a_coarseDomain,
                                    const Vector<DisjointBoxLayout>&               a_grids,
                                    const Vector<int>&                             a_refRatios,
                                    const Real&                                    a_coarsedx,
                                    BCHolder                                       a_bc,
                                    const Real&                                    a_alpha,
                                    Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                                    const Real&                                    a_beta,
                                    Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef,
                                    Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_cCoef)
{
  CH_TIME("DarcyBrinkmanOpFactory::define");

  setDefaultValues();

  m_boxes = a_grids;

  m_refRatios = a_refRatios;

  m_bc = a_bc;

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

  m_cCoef = a_cCoef;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// AMR Factory define function, with coefficient data allocated automagically
// for operators.
void
DarcyBrinkmanOpFactory::
define(const ProblemDomain& a_coarseDomain,
       const Vector<DisjointBoxLayout>& a_grids,
       const Vector<int>& a_refRatios,
       const Real& a_coarsedx,
       BCHolder a_bc,
       const IntVect& a_ghostVect)
{
  // This just allocates coefficient data, sets alpha = beta = 1, and calls
  // the other define() method.
  Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(a_grids.size());
  Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(a_grids.size());
  Vector<RefCountedPtr<LevelData<FArrayBox> > > cCoef(a_grids.size());
  for (int i = 0; i < a_grids.size(); ++i)
  {
    aCoef[i] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(a_grids[i], 1, a_ghostVect));
    bCoef[i] = RefCountedPtr<LevelData<FluxBox> >(
        new LevelData<FluxBox>(a_grids[i], 1, a_ghostVect));
    cCoef[i] = RefCountedPtr<LevelData<FArrayBox> >(
        new LevelData<FArrayBox>(a_grids[i], 1, a_ghostVect));

    // Initialize the a and b coefficients to 1 for starters.
    for (DataIterator dit = aCoef[i]->dataIterator(); dit.ok(); ++dit)
    {
      (*aCoef[i])[dit()].setVal(1.0);
      (*cCoef[i])[dit()].setVal(1.0);
      for (int idir = 0; idir < SpaceDim; ++idir)
        (*bCoef[i])[dit()][idir].setVal(1.0);
    }
  }
  Real alpha = 1.0, beta = 1.0;
  define(a_coarseDomain, a_grids, a_refRatios, a_coarsedx, a_bc,
         alpha, aCoef, beta, bCoef, cCoef);
}
//-----------------------------------------------------------------------

MGLevelOp<LevelData<FArrayBox> >* DarcyBrinkmanOpFactory::MGnewOp(const ProblemDomain& a_indexSpace,
                                                                  int                  a_depth,
                                                                  bool                 a_homoOnly)
{
  CH_TIME("DarcyBrinkmanOpFactory::MGnewOp");

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

  if (coarsening > 1 && !m_boxes[ref].coarsenable(coarsening*DarcyBrinkmanOp::s_maxCoarse))
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

  DarcyBrinkmanOp* newOp = new DarcyBrinkmanOp;

  // Can we just do AMR define with a coarse level?
   int refRatio = 2; // default MG refinement
   if (ref-1 < m_refRatios.size())
   {
     refRatio = m_refRatios[ref-1]; // AMR refinement if defined, else default back to MG
   }

//  newOp->define(layout, dx, domain, m_bc, ex, cfregion);
  // ref > 0 means we've done AMR coarsening i.e. there are coarser amr levels
   if (ref > 0)
   {
     newOp->define(layout, layoutCrse, dx,
                   refRatio,
                   domain, m_bc,
                   ex, cfregion);
   }
   else
   {
     // no coarse level define
     newOp->define(layout, dx, domain, m_bc, ex, cfregion);

   }

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;

  if (a_depth == 0)
  {
    // don't need to coarsen anything for this
    newOp->m_aCoef = m_aCoef[ref];
    newOp->m_bCoef = m_bCoef[ref];
    newOp->m_cCoef = m_cCoef[ref];
  }

  else
  {
    // need to coarsen coefficients
    RefCountedPtr<LevelData<FArrayBox> > aCoef( new LevelData<FArrayBox> );
    RefCountedPtr<LevelData<FArrayBox> > cCoef( new LevelData<FArrayBox> );
    RefCountedPtr<LevelData<FluxBox> > bCoef( new LevelData<FluxBox> );
    aCoef->define(layout, m_aCoef[ref]->nComp(), m_aCoef[ref]->ghostVect());
    bCoef->define(layout, m_bCoef[ref]->nComp(), m_bCoef[ref]->ghostVect());
    cCoef->define(layout, m_cCoef[ref]->nComp(), m_cCoef[ref]->ghostVect());

    // average coefficients to coarser level
    // for now, do this with a CoarseAverage --
    // may want to switch to harmonic averaging at some point
    CoarseAverage averager(m_aCoef[ref]->getBoxes(),
                           layout, aCoef->nComp(), coarsening);

    CoarseAverageFace faceAverager(m_bCoef[ref]->getBoxes(),
                                   bCoef->nComp(), coarsening);

    if (m_coefficient_average_type == CoarseAverage::arithmetic)
    {
      averager.averageToCoarse(*aCoef, *(m_aCoef[ref]));
      averager.averageToCoarse(*cCoef, *(m_cCoef[ref]));
      faceAverager.averageToCoarse(*bCoef, *(m_bCoef[ref]));
    }
    else if (m_coefficient_average_type == CoarseAverage::harmonic)
    {
      averager.averageToCoarseHarmonic(*aCoef, *(m_aCoef[ref]));
      averager.averageToCoarseHarmonic(*cCoef, *(m_cCoef[ref]));
      faceAverager.averageToCoarseHarmonic(*bCoef, *(m_bCoef[ref]));
    }
    else
    {
      MayDay::Abort("DarcyBrinkmanOpFactory::MGNewOp -- bad averagetype");
    }

    newOp->m_aCoef = aCoef;
    newOp->m_bCoef = bCoef;
    newOp->m_cCoef = cCoef;
  }

  newOp->computeLambda();

  newOp->m_dxCrse = dxCrse;

  return (MGLevelOp<LevelData<FArrayBox> >*)newOp;
}

AMRLevelOp<LevelData<FArrayBox> >* DarcyBrinkmanOpFactory::AMRnewOp(const ProblemDomain& a_indexSpace)
{
  CH_TIME("DarcyBrinkmanOpFactory::AMRnewOp");

  DarcyBrinkmanOp* newOp = new DarcyBrinkmanOp;
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
                  m_exchangeCopiers[ref], m_cfregion[ref]);
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
                  m_exchangeCopiers[ref], m_cfregion[ref]);
  }

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;

  newOp->m_aCoef = m_aCoef[ref];
  newOp->m_bCoef = m_bCoef[ref];
  newOp->m_cCoef = m_cCoef[ref];

  newOp->computeLambda();

  newOp->m_dxCrse = dxCrse;

  return (AMRLevelOp<LevelData<FArrayBox> >*)newOp;
}

int DarcyBrinkmanOpFactory::refToFiner(const ProblemDomain& a_domain) const
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
void DarcyBrinkmanOpFactory::setDefaultValues()
{
  // Default to Laplacian operator
  m_alpha = 0.0;
  m_beta = -1.0;

  m_coefficient_average_type = CoarseAverage::arithmetic;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
DarcyBrinkmanOp::
finerOperatorChanged(const MGLevelOp<LevelData<FArrayBox> >& a_operator,
                     int a_coarseningFactor)
{
  const DarcyBrinkmanOp& op =
      dynamic_cast<const DarcyBrinkmanOp&>(a_operator);

  // Perform multigrid coarsening on the operator data.
  LevelData<FArrayBox>& acoefCoar = *m_aCoef;
  LevelData<FArrayBox>& ccoefCoar = *m_cCoef;
  const LevelData<FArrayBox>& acoefFine = *(op.m_aCoef);
  const LevelData<FArrayBox>& ccoefFine = *(op.m_cCoef);
  LevelData<FluxBox>& bcoefCoar = *m_bCoef;
  const LevelData<FluxBox>& bcoefFine = *(op.m_bCoef);
  if (a_coarseningFactor != 1)
  {
    CoarseAverage cellAverage(acoefFine.disjointBoxLayout(),
                              acoefCoar.disjointBoxLayout(),
                              1, a_coarseningFactor);
    for (DataIterator dit = acoefCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)

    {
      acoefCoar[dit()].setVal(0.);
      ccoefCoar[dit()].setVal(0.);
      cellAverage.averageToCoarse(acoefCoar, acoefFine);
      cellAverage.averageToCoarse(ccoefCoar, ccoefFine);
    }

    CoarseAverageFace faceAverage(bcoefFine.disjointBoxLayout(),
                                  1, a_coarseningFactor);
    for (DataIterator dit = bcoefCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)

    {
      bcoefCoar[dit()].setVal(0.);
      faceAverage.averageToCoarse(bcoefCoar, bcoefFine);
    }
  }

  // Handle inter-box ghost cells.
  acoefCoar.exchange();
  bcoefCoar.exchange();
  ccoefCoar.exchange();

  // Mark the relaxation coefficient dirty.
  m_lambdaNeedsResetting = true;

  // Notify any observers of this change.
  notifyObserversOfChange();
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"



