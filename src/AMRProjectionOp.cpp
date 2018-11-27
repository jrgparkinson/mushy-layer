/*
 * AMRProjectionOp.cpp
 *
 *  Created on: 16 Nov 2018
 *      Author: parkinsonjl
 */


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
#include "AMRProjectionOpF_F.H"
#include "VCAMRPoissonOpF_F.H"

#include "AMRProjectionOp.h"

#include "NamespaceHeader.H"

AMRProjectionOp::AMRProjectionOp ()
{
  // TODO Auto-generated constructor stub

}

AMRProjectionOp::~AMRProjectionOp ()
{
  // TODO Auto-generated destructor stub
}

void AMRProjectionOp::resetLambda()
{
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

//      lambdaFab.plus(-scale*m_beta*pow(2, SpaceDim));

      for (int dir = 0; dir < SpaceDim; dir++)
      {
        FORT_SUMFACES(CHF_FRA(lambdaFab),
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


void AMRProjectionOp::prolongIncrement(LevelData<FArrayBox>&       a_phiThisLevel,
                                    const LevelData<FArrayBox>& a_correctCoarse)
{
  CH_TIME("AMRProjectionOp::prolongIncrement");

  DisjointBoxLayout dbl = a_phiThisLevel.disjointBoxLayout();
  int mgref = 2; //this is a multigrid func
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

//        FORT_PROLONG_2(CHF_FRA_SHIFT(phi, iv),
//                     CHF_CONST_FRA_SHIFT(coarse, civ),
//                     CHF_BOX_SHIFT(region, iv),
//                     CHF_CONST_INT(mgref));

//        FORT_PROLONG_JP(CHF_FRA_SHIFT(phi, iv),
//                             CHF_CONST_FRA_SHIFT(coarse, civ),
//                             CHF_BOX_SHIFT(region, iv),
//                             CHF_CONST_INT(mgref));

           FORT_PROLONG(CHF_FRA_SHIFT(phi, iv),
                                     CHF_CONST_FRA_SHIFT(coarse, civ),
                                     CHF_BOX_SHIFT(region, iv),
                                     CHF_CONST_INT(mgref));
      }


  }//end pragma
}


// Factory
AMRProjectionOpFactory::AMRProjectionOpFactory()
{
  setDefaultValues();
}



//-----------------------------------------------------------------------
//  AMR Factory define function
void AMRProjectionOpFactory::define(const ProblemDomain&                           a_coarseDomain,
                                   const Vector<DisjointBoxLayout>&               a_grids,
                                   const Vector<int>&                             a_refRatios,
                                   const Real&                                    a_coarsedx,
                                   BCHolder                                       a_bc,
                                   const Real&                                    a_alpha,
                                   Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                                   const Real&                                    a_beta,
                                   Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef,
                                   int                                            a_averaging_type)
{
  CH_TIME("AMRProjectionOpFactory::define");

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

  m_coefficient_average_type = a_averaging_type;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// AMR Factory define function, with coefficient data allocated automagically
// for operators.
void
AMRProjectionOpFactory::
define(const ProblemDomain& a_coarseDomain,
       const Vector<DisjointBoxLayout>& a_grids,
       const Vector<int>& a_refRatios,
       const Real& a_coarsedx,
       BCHolder a_bc,
       const IntVect& a_ghostVect,
       int a_averaging_type)
{
  // This just allocates coefficient data, sets alpha = beta = 1, and calls
  // the other define() method.
  Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(a_grids.size());
  Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(a_grids.size());
  for (int i = 0; i < a_grids.size(); ++i)
  {
    aCoef[i] = RefCountedPtr<LevelData<FArrayBox> >(
                 new LevelData<FArrayBox>(a_grids[i], 1, a_ghostVect));
    bCoef[i] = RefCountedPtr<LevelData<FluxBox> >(
                 new LevelData<FluxBox>(a_grids[i], 1, a_ghostVect));

    // Initialize the a and b coefficients to 1 for starters.
    for (DataIterator dit = aCoef[i]->dataIterator(); dit.ok(); ++dit)
    {
      (*aCoef[i])[dit()].setVal(1.0);
      for (int idir = 0; idir < SpaceDim; ++idir)
        (*bCoef[i])[dit()][idir].setVal(1.0);
    }
  }
  Real alpha = 1.0, beta = 1.0;
  define(a_coarseDomain, a_grids, a_refRatios, a_coarsedx, a_bc,
         alpha, aCoef, beta, bCoef, a_averaging_type);
}
//-----------------------------------------------------------------------

MGLevelOp<LevelData<FArrayBox> >* AMRProjectionOpFactory::MGnewOp(const ProblemDomain& a_indexSpace,
                                                                  int                  a_depth,
                                                                  bool                 a_homoOnly)
{
  CH_TIME("AMRProjectionOpFactory::MGnewOp");

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

  if (coarsening > 1 && !m_boxes[ref].coarsenable(coarsening*AMRProjectionOp::s_maxCoarse))
  {
    return NULL;
  }

  dx *= coarsening;

  DisjointBoxLayout layout;
  coarsen_dbl(layout, m_boxes[ref], coarsening);

  Copier ex = m_exchangeCopiers[ref];
  CFRegion cfregion = m_cfregion[ref];

  if (coarsening > 1)
  {
    ex.coarsen(coarsening);
    cfregion.coarsen(coarsening);
  }

  AMRProjectionOp* newOp = new AMRProjectionOp;

  newOp->define(layout, dx, domain, m_bc, ex, cfregion);

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;

  if (a_depth == 0)
    {
      // don't need to coarsen anything for this
      newOp->m_aCoef = m_aCoef[ref];
      newOp->m_bCoef = m_bCoef[ref];
    }
  else
    {
      // need to coarsen coefficients
      RefCountedPtr<LevelData<FArrayBox> > aCoef( new LevelData<FArrayBox> );
      RefCountedPtr<LevelData<FluxBox> > bCoef( new LevelData<FluxBox> );
      aCoef->define(layout, m_aCoef[ref]->nComp(), m_aCoef[ref]->ghostVect());
      bCoef->define(layout, m_bCoef[ref]->nComp(), m_bCoef[ref]->ghostVect());

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
          faceAverager.averageToCoarse(*bCoef, *(m_bCoef[ref]));
        }
      else if (m_coefficient_average_type == CoarseAverage::harmonic)
        {
          averager.averageToCoarseHarmonic(*aCoef, *(m_aCoef[ref]));
          faceAverager.averageToCoarseHarmonic(*bCoef, *(m_bCoef[ref]));
        }
      else if (m_coefficient_average_type == CoarseAverage::geometric)
             {
               averager.averageToCoarseGeometric(*aCoef, *(m_aCoef[ref]));
               faceAverager.averageToCoarseGeometric(*bCoef, *(m_bCoef[ref]));
             }
      else
        {
          MayDay::Abort("AMRProjectionOpFactory::MGNewOp -- bad averagetype");
        }

      newOp->m_aCoef = aCoef;
      newOp->m_bCoef = bCoef;
    }

  newOp->computeLambda();

  newOp->m_dxCrse = dxCrse;

  newOp->s_relaxMode = m_relaxMode;

  return (MGLevelOp<LevelData<FArrayBox> >*)newOp;
}

AMRLevelOp<LevelData<FArrayBox> >* AMRProjectionOpFactory::AMRnewOp(const ProblemDomain& a_indexSpace)
{
  CH_TIME("AMRProjectionOpFactory::AMRnewOp");

  AMRProjectionOp* newOp = new AMRProjectionOp;
  Real dxCrse = -1.0;

  // Need to know how many components we're solving for
  int nComp = 1;

  // Find number of comps from m_bCoef, which should be defined on at least one level
  for (int lev = 0; lev <m_bCoef.size(); lev++)
  {
    if (m_bCoef[lev] != NULL)
    {
      nComp = m_bCoef[lev]->nComp();
      break;
    }

    // Should never reach this point
    CH_assert("AMRProjectionOpFactory::AMRnewOp - m_bCoef not defined");
  }

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
                        m_exchangeCopiers[0], m_cfregion[0],
                        nComp);
    }
  }
  else if (ref ==  m_domains.size()-1)
  {
    dxCrse = m_dx[ref-1];

    // finest AMR level
    newOp->define(m_boxes[ref], m_boxes[ref-1], m_dx[ref],
                  m_refRatios[ref-1],
                  a_indexSpace, m_bc,
                  m_exchangeCopiers[ref], m_cfregion[ref],
                  nComp);
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
                  m_exchangeCopiers[ref], m_cfregion[ref],
                  nComp);
  }

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;

  newOp->m_aCoef = m_aCoef[ref];
  newOp->m_bCoef = m_bCoef[ref];

  if (newOp->m_aCoef != NULL)
  {
    newOp->computeLambda();
  }

  newOp->m_dxCrse = dxCrse;

  newOp->s_relaxMode = m_relaxMode;

  return (AMRLevelOp<LevelData<FArrayBox> >*)newOp;
}

int AMRProjectionOpFactory::refToFiner(const ProblemDomain& a_domain) const
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
void AMRProjectionOpFactory::setDefaultValues()
{
  // Default to Laplacian operator
  m_alpha = 0.0;
  m_beta = -1.0;

  m_coefficient_average_type = CoarseAverage::arithmetic;
//  m_coefficient_average_type = CoarseAverage::harmonic;
  m_relaxMode = 1; // gsrb
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------


#include "NamespaceFooter.H"
