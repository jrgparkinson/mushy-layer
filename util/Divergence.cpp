#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Divergence.H"
#include "DivergenceF_F.H"
#include "LevelFluxRegister.H"
#include "CellToEdge.H"
#include "BoxIterator.H"

/// Computes divergence of fields containing multiple components
void
Divergence::levelDivergenceMACMultiComp(LevelData<FArrayBox>& a_div,
                                        const LevelData<FluxBox>& a_uEdge,
                                        const Real a_dx)
{
  int numComps = a_uEdge.nComp();

  DataIterator dit = a_div.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    a_div[dit()].setVal(0.0);

    const FluxBox& thisFluxBox = a_uEdge[dit()];
    Box cellBox(thisFluxBox.box());
    // just to be sure we don't accidentally trash memory...
    cellBox &= a_div[dit()].box();

    // now loop over coordinate directions and add to divergence
    for (int dir=0; dir<SpaceDim; dir++)
    {
      const FArrayBox& uEdgeDir = thisFluxBox[dir];

      for (int comp=0; comp<numComps; comp++)
      {

        FORT_DIVERGENCE(CHF_CONST_FRA1(uEdgeDir,comp),
                        CHF_FRA1(a_div[dit()], comp),
                        CHF_BOX(cellBox),
                        CHF_CONST_REAL(a_dx),
                        CHF_INT(dir));

//        int temp=0;

        // For debugging:
//#if DEBUG

        FArrayBox divDir(a_div[dit].box(), 1);
        divDir.setVal(0.0);
        FORT_DIVERGENCE(CHF_CONST_FRA1(uEdgeDir,comp),
                                CHF_FRA1(divDir, 0),
                                CHF_BOX(cellBox),
                                CHF_CONST_REAL(a_dx),
                                CHF_INT(dir));
//        int temp2 = 0;
//#endif


      }//  end loop over components
    } // end loop over spatial directions
  } // end loop over boxes

}
// ---------------------------------------------------------
void
Divergence::levelDivergenceMAC(LevelData<FArrayBox>& a_div,
                               const LevelData<FluxBox>& a_uEdge,
                               const Real a_dx)
{
  // for now, force this to be single-component, and
  // uEdge also single-component, with each component being the
  // normal component
  CH_assert (a_uEdge.nComp() == 1);
  CH_assert (a_div.nComp() == 1);

  int comp = 0;

  DataIterator dit = a_div.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    a_div[dit()].setVal(0.0);

    const FluxBox& thisFluxBox = a_uEdge[dit()];
    Box cellBox(thisFluxBox.box());
    // just to be sure we don't accidentally trash memory...
    cellBox &= a_div[dit()].box();

    // now loop over coordinate directions and add to divergence
    for (int dir=0; dir<SpaceDim; dir++)
    {
      const FArrayBox& uEdgeDir = thisFluxBox[dir];

      FORT_DIVERGENCE(CHF_CONST_FRA1(uEdgeDir,comp),
                      CHF_FRA1(a_div[dit()], comp),
                      CHF_BOX(cellBox),
                      CHF_CONST_REAL(a_dx),
                      CHF_INT(dir));

//      FArrayBox divDir(a_div[dit].box(), 1);
//      divDir.setVal(0.0);
//
//      FORT_DIVERGENCE(CHF_CONST_FRA1(uEdgeDir,comp),
//                      CHF_FRA1(divDir, comp),
//                      CHF_BOX(cellBox),
//                      CHF_CONST_REAL(a_dx),
//                      CHF_INT(dir));

    }
  }
}

void
Divergence::simpleDivergenceMAC( FArrayBox& a_div, const FluxBox& a_uEdge,
                                 const Real a_dx)
{

  // for now, force this to be single-component, and
  // uEdge also single-component, with each component being the
  // normal component
  CH_assert (a_uEdge.nComp() == 1);
  CH_assert (a_div.nComp() == 1);

  int comp = 0;

  a_div.setVal(0.0);
  const Box& cellBox = a_uEdge.box();

  // now loop over coordinate directions and increment divergence
  for (int dir=0; dir<SpaceDim; dir++)
  {
    const FArrayBox& uEdgeDir = a_uEdge[dir];

    FORT_DIVERGENCE(CHF_CONST_FRA1(uEdgeDir, comp),
                    CHF_FRA1(a_div, comp),
                    CHF_BOX(cellBox),
                    CHF_CONST_REAL(a_dx),
                    CHF_INT(dir));
  }
}

// ---------------------------------------------------------
void
Divergence::compDivergenceMAC(LevelData<FArrayBox>& a_div,
                              LevelData<FluxBox>& a_uEdge,
                              LevelData<FluxBox>* a_uEdgeFinePtr,
                              const Real a_dx,
                              const Real* a_dxFine,
                              const int a_nRefine,
                              const Box& a_dProblem)
{
  ProblemDomain physdomain(a_dProblem);
  compDivergenceMAC(a_div, a_uEdge, a_uEdgeFinePtr, a_dx,
                    a_dxFine, a_nRefine, physdomain);
}

// ---------------------------------------------------------
void
Divergence::compDivergenceMAC(LevelData<FArrayBox>& a_div,
                              LevelData<FluxBox>& a_uEdge,
                              LevelData<FluxBox>* a_uEdgeFinePtr,
                              const Real a_dx,
                              const Real* a_dxFine,
                              const int a_nRefine,
                              const ProblemDomain& a_dProblem)
{
  if (a_uEdgeFinePtr != nullptr)
  {
    // define a LevelFluxRegister to do coarse-fine mismatch accounting
    const DisjointBoxLayout& dblCrse = a_div.getBoxes();
    const DisjointBoxLayout& dblFine = a_uEdgeFinePtr->getBoxes();
    int ncomp = 1;

    ProblemDomain fineDomain = a_dProblem;
    fineDomain.refine(a_nRefine);
    LevelFluxRegister FR(dblFine, dblCrse, fineDomain, a_nRefine,ncomp);
    compDivergenceMAC(a_div, a_uEdge, a_uEdgeFinePtr, &FR,
                      a_dx, a_dxFine, a_nRefine, a_dProblem);
  }
  else
  {
    // do single-level version
    levelDivergenceMAC(a_div, a_uEdge, a_dx);
  }
}

// ---------------------------------------------------------
void
Divergence::compDivergenceMAC(LevelData<FArrayBox>& a_div,
                              LevelData<FluxBox>& a_uEdge,
                              LevelData<FluxBox>* a_uEdgeFinePtr,
                              LevelFluxRegister* a_fluxRegPtr,
                              const Real a_dx,
                              const Real* a_dxFine,
                              const int a_nRefine,
                              const Box& a_dProblem)
{
  ProblemDomain physdomain(a_dProblem);
  compDivergenceMAC(a_div, a_uEdge, a_uEdgeFinePtr, a_fluxRegPtr,
                    a_dx, a_dxFine, a_nRefine, physdomain);
}

// ---------------------------------------------------------
void
Divergence::compDivergenceMAC(LevelData<FArrayBox>& a_div,
                              LevelData<FluxBox>& a_uEdge,
                              LevelData<FluxBox>* a_uEdgeFinePtr,
                              LevelFluxRegister* a_fluxRegPtr,
                              const Real a_dx,
                              const Real* a_dxFine,
                              const int a_nRefine,
                              const ProblemDomain& a_dProblem)
{
  // first do simple single-level divergence, then fix up if
  // necessary
  levelDivergenceMAC(a_div, a_uEdge, a_dx);

  // for now, hardwire to simple single-component case
  CH_assert(a_div.nComp() == 1);

  // now adjust for effect of finer level (if applicable)
  if (a_uEdgeFinePtr != nullptr)
  {
    CH_assert(a_fluxRegPtr != nullptr);
    LevelFluxRegister& FR = *a_fluxRegPtr;
    FR.setToZero();

    DataIterator dit = a_div.dataIterator();

    // iterate over coarse boxes
    for (dit.reset(); dit.ok(); ++dit)
    {
      FluxBox& thisFlux = a_uEdge[dit()];
      Real scale = 1.0;

      const Interval compInterval(0,0);
      // iterate over directions
      for (int dir=0; dir<SpaceDim; dir++)
      {
        FArrayBox& thisDirFlux = thisFlux[dir];
        FR.incrementCoarse(thisDirFlux, scale, dit(),
                           compInterval,compInterval,dir);
      }
    }

    // now do fine-level data
    DataIterator ditFine = a_uEdgeFinePtr->dataIterator();
    LevelData<FluxBox>& uEdgeFine = *a_uEdgeFinePtr;

    for (ditFine.reset(); ditFine.ok(); ++ditFine)
    {
      FluxBox& thisFineFlux = uEdgeFine[ditFine()];
      Real scale = 1.0;
      Interval srcComps(0,0);
      // iterate over directions
      for (int dir=0; dir<SpaceDim; dir++)
      {
        FArrayBox& thisDirFineFlux = thisFineFlux[dir];
        // iterate over lo-hi sides (do i really need to do this?)
        SideIterator sit;
        for (sit.begin(); sit.ok(); sit.next())
        {
          Side::LoHiSide hiorlo = sit();
          FR.incrementFine(thisDirFineFlux,scale,ditFine(),srcComps,
                           srcComps,dir,hiorlo);
        }
      }
    }

    // now do reflux-divergence to put correct coarse-level divergence
    // this will only work for dx == dy.
    Real scale = 1.0/a_dx;
    FR.reflux(a_div, scale);
  } // end correction for finer level
}

// ---------------------------------------------------------
void
Divergence::levelDivergenceCC(LevelData<FArrayBox>& a_div,
                              LevelData<FArrayBox>& a_u,
                              LevelData<FArrayBox>* a_uCrsePtr,
                              const Real a_dx,
                              const int a_nRefCrse,
                              const Box& a_dProblem,
                              const bool a_quadInterp)
{
  ProblemDomain physdomain(a_dProblem);
  levelDivergenceCC(a_div, a_u, a_uCrsePtr, a_dx, a_nRefCrse,
                    physdomain, a_quadInterp);
}

// ---------------------------------------------------------
void
Divergence::levelDivergenceCC(LevelData<FArrayBox>& a_div,
                              LevelData<FArrayBox>& a_u,
                              LevelData<FArrayBox>* a_uCrsePtr,
                              const Real a_dx,
                              const int a_nRefCrse,
                              const ProblemDomain& a_dProblem,
                              const bool a_quadInterp)
{
  // compute coarse-fine boundary conditions (assume physical ones
  // already done)

  if (a_uCrsePtr != nullptr)
  {
    if (a_quadInterp)
    {
      const DisjointBoxLayout& boxes = a_div.getBoxes();
      const DisjointBoxLayout& coarseBoxes = a_uCrsePtr->getBoxes();
      int nComp = a_u.nComp();
      QuadCFInterp cfInterp(boxes, &coarseBoxes, a_dx, a_nRefCrse,
                            nComp, a_dProblem);

      levelDivergenceCC(a_div, a_u, a_uCrsePtr, a_dx,
                        a_quadInterp, cfInterp);
    }
    else
    {
      // non-quadratic CF interp bc's not implemented yet
      CH_assert(a_quadInterp);
    }
  }
  else
  {
    // for single-level, won't need C/F interpolation anyway...
    QuadCFInterp cfInterpBogus;
    levelDivergenceCC(a_div, a_u, a_uCrsePtr, a_dx, a_quadInterp,
                      cfInterpBogus);
  }
}

void
Divergence::levelDivergenceCC(LevelData<FArrayBox>& a_div,
                              LevelData<FArrayBox>& a_u,
                              const Real a_dx)
{
  QuadCFInterp cfInterp;
  levelDivergenceCC(a_div, a_u, nullptr, a_dx, false, cfInterp);
}

// ---------------------------------------------------------
void
Divergence::levelDivergenceCC(LevelData<FArrayBox>& a_div,
                              LevelData<FArrayBox>& a_u,
                              LevelData<FArrayBox>* a_uCrsePtr,
                              const Real a_dx,
                              const bool a_quadInterp,
                              QuadCFInterp& a_cfInterp)
{
  // for now, hardwire to simple single-component case
  CH_assert (a_div.nComp() == 1);
  CH_assert (a_u.nComp() == SpaceDim);


  // compute coarse-fine BC data
  if (a_uCrsePtr != nullptr)
  {
    if (a_quadInterp)
    {
      CH_assert(a_cfInterp.isDefined());

      a_cfInterp.coarseFineInterp(a_u, *a_uCrsePtr);
    }
    else
    {
      // non-quadInterp not implemented yet
      CH_assert(a_quadInterp);
    }
  } // end if coarser level exists

  // cell-centered divergence is DivMac(averageCellToEdge(uCell))
  const DisjointBoxLayout& boxes = a_u.getBoxes();

  LevelData<FluxBox> uEdge(boxes, 1);

  CellToEdge(a_u, uEdge);

  //	for (DataIterator dit = uEdge.dataIterator(); dit.ok(); ++dit)
  //	{
  //		FluxBox& edge = uEdge[dit];
  //	}

  levelDivergenceMAC(a_div, uEdge, a_dx);
}


void Divergence::levelDivergenceCCNew(LevelData<FArrayBox>& a_div,
                                      LevelData<FArrayBox>& a_u,
                                      LevelData<FArrayBox>* a_uCrsePtr,
                                      const Real a_dx,
                                      const bool a_quadInterp,
                                      QuadCFInterp& a_cfInterp)
{
  CH_TIME("Divergence::levelDivergenceCCNew");

  // for now, hardwire to simple single-component case
  CH_assert (a_div.nComp() == 1);
  CH_assert (a_u.nComp() == SpaceDim);

  // compute coarse-fine BC data
  if (a_uCrsePtr != nullptr and a_quadInterp)
  {
    CH_assert(a_cfInterp.isDefined());

    a_cfInterp.coarseFineInterp(a_u, *a_uCrsePtr);
  }


  int order = 2;

  int num_ghost = a_div.ghostVect()[0];


  for (DataIterator dit = a_div.dataIterator(); dit.ok(); ++dit)
  {

    a_div[dit].setVal(0);

    for( int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox a_uDir(a_u[dit].box(), 1);
      a_uDir.copy(a_u[dit], idir, 0, 1);

      FArrayBox a_divDir(a_div[dit].box(), 1);
      a_divDir.setVal(0);

      Box b = a_div[dit].box();


      b.grow(-num_ghost);

      b.growHi(idir, -order);

      for (SideIterator side; side.ok(); ++side)
      {
        BoxIterator bit;
        if (side() == Side::Lo)
        {
          bit = BoxIterator(b);
        }
        else if (side() == Side::Hi)
        {
          Box region = adjCellBox(b, idir, side(), order);
          bit = BoxIterator(region);
        }

        int sign = -1;
        if (side() == Side::Lo)
        {
          sign = 1;
        }

        for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          Real f_0, f_1, f_2, f_3, f_4;

          f_0 = a_u[dit](iv, idir);
          f_1 = a_u[dit](iv + sign*BASISV(idir), idir);
          f_2 = a_u[dit](iv + sign*2*BASISV(idir), idir);

          if (order > 2)
          {
            f_3 = a_u[dit](iv + sign*3*BASISV(idir), idir);
          }
          if (order > 3)
          {
            f_4 = a_u[dit](iv + sign*4*BASISV(idir), idir);
          }

          Real grad;

          if (order == 4)
          {
            grad = sign * (-25*f_0 + 48*f_1 -36*f_2 + 16*f_3 - 3*f_4)/(12*a_dx);
          }
          else if (order == 3)
          {
            grad = sign * (-11*f_0 + 18*f_1 - 9*f_2 + 2*f_3)/(6*a_dx);
          }
          else
          {
            grad = sign*(-3*f_0 + 4*f_1 - f_2)/(2*a_dx);
          }


          a_div[dit](iv, 0) += grad;
          a_divDir(iv, 0) += grad;
        } // end box iterator

      } // end side iterator

    } //direction iterator


  }

}

// ---------------------------------------------------------
void
Divergence::compDivergenceCC(LevelData<FArrayBox>& a_div,
                             LevelData<FArrayBox>& a_u,
                             LevelData<FArrayBox>* a_uCrsePtr,
                             LevelData<FArrayBox>* a_uFinePtr,
                             const Real a_dx,
                             const int a_nRefCrse,
                             const int a_nRefFine,
                             const Box& a_dProblem,
                             const bool a_quadInterp)
{
  ProblemDomain physdomain(a_dProblem);
  compDivergenceCC(a_div, a_u, a_uCrsePtr, a_uFinePtr, a_dx,
                   a_nRefCrse, a_nRefFine, physdomain, a_quadInterp);
}

// ---------------------------------------------------------
void
Divergence::compDivergenceCC(LevelData<FArrayBox>& a_div,
                             LevelData<FArrayBox>& a_u,
                             LevelData<FArrayBox>* a_uCrsePtr,
                             LevelData<FArrayBox>* a_uFinePtr,
                             const Real a_dx,
                             const int a_nRefCrse,
                             const int a_nRefFine,
                             const ProblemDomain& a_dProblem,
                             const bool a_quadInterp)
{
  QuadCFInterp cfInterpCrse;
  QuadCFInterp cfInterpFine;
  LevelFluxRegister FR;

  const DisjointBoxLayout& thisLevelBoxes = a_u.getBoxes();

  // for now, hardwire to single-component
  CH_assert(a_div.nComp() == 1);
  int nComp = 1;

  // define coarse-level C/f-BC object
  if (a_uCrsePtr != nullptr)
  {
    CH_assert(a_nRefCrse > 0);

    const DisjointBoxLayout& crseLevelBoxes = a_uCrsePtr->getBoxes();
    cfInterpCrse.define(thisLevelBoxes, &crseLevelBoxes, a_dx, a_nRefCrse,
                        SpaceDim*nComp, a_dProblem);
  }

  // define fine-level CF-BC object and flux register
  if (a_uFinePtr != nullptr)
  {
    CH_assert(a_nRefFine > 0);

    const DisjointBoxLayout& fineLevelBoxes = a_uFinePtr->getBoxes();
    ProblemDomain fineDomain(a_dProblem);
    fineDomain.refine(a_nRefFine);
    // this only works if dx == dy
    Real dxFine = a_dx/a_nRefFine;

    cfInterpFine.define(fineLevelBoxes, &thisLevelBoxes, dxFine,
                        a_nRefFine, SpaceDim*nComp, fineDomain);

    FR.define(fineLevelBoxes, thisLevelBoxes, fineDomain,
              a_nRefFine, nComp);
  }

  // now call composite divergence operator
  compDivergenceCC(a_div, a_u, a_uCrsePtr, a_uFinePtr, a_dx,
                   a_nRefFine, a_nRefCrse, a_dProblem,
                   a_quadInterp, &FR, cfInterpCrse,
                   cfInterpFine);
}

// ---------------------------------------------------------
void
Divergence::compDivergenceCC(LevelData<FArrayBox>& a_div,
                             LevelData<FArrayBox>& a_u,
                             LevelData<FArrayBox>* a_uCrsePtr,
                             LevelData<FArrayBox>* a_uFinePtr,
                             const Real a_dx,
                             const int a_nRefCrse,
                             const int a_nRefFine,
                             const Box& a_dProblem,
                             const bool a_quadInterp,
                             LevelFluxRegister* a_fluxRegFinePtr,
                             QuadCFInterp& a_cfInterpCrse,
                             QuadCFInterp& a_cfInterpFine)
{
  ProblemDomain physdomain(a_dProblem);
  compDivergenceCC(a_div, a_u, a_uCrsePtr, a_uFinePtr, a_dx,
                   a_nRefCrse, a_nRefFine, physdomain, a_quadInterp,
                   a_fluxRegFinePtr, a_cfInterpCrse, a_cfInterpFine);
}

// ---------------------------------------------------------
void
Divergence::compDivergenceCC(LevelData<FArrayBox>& a_div,
                             LevelData<FArrayBox>& a_u,
                             LevelData<FArrayBox>* a_uCrsePtr,
                             LevelData<FArrayBox>* a_uFinePtr,
                             const Real a_dx,
                             const int a_nRefCrse,
                             const int a_nRefFine,
                             const ProblemDomain& a_dProblem,
                             const bool a_quadInterp,
                             LevelFluxRegister* a_fluxRegFinePtr,
                             QuadCFInterp& a_cfInterpCrse,
                             QuadCFInterp& a_cfInterpFine)
{

  // for now, hardwire to simplest single-component case
  CH_assert(a_div.nComp() == 1);
  CH_assert(a_u.nComp() == SpaceDim);

  // first do coarse-level BC's
  if (a_uCrsePtr != nullptr)
  {
    if (a_quadInterp)
    {
      CH_assert(a_cfInterpCrse.isDefined());
      a_cfInterpCrse.coarseFineInterp(a_u, *a_uCrsePtr);
    }
    else
    {
      // non-quadInterp BC's not yet implemented
      CH_assert(a_quadInterp);
    }
  }

  // now average cells->edges
  const DisjointBoxLayout& boxes = a_u.getBoxes();
  LevelData<FluxBox> uEdge(boxes,1);
  CellToEdge(a_u, uEdge);

  // compute no-fine-level divergence
  levelDivergenceMAC(a_div, uEdge, a_dx);

  // if a fine level exists, fix up at C/F interface:
  if (a_uFinePtr != nullptr)
  {
    // dereference the pointer for convenience
    LevelData<FArrayBox>& uFine = *a_uFinePtr;
    LevelFluxRegister& FR = *a_fluxRegFinePtr;
    // initialize flux register
    FR.setToZero();

    // first do C/F BC's on finer level data
    if (a_quadInterp)
    {
      CH_assert(a_cfInterpFine.isDefined());
      a_cfInterpFine.coarseFineInterp(uFine, a_u);
    }
    else
    {
      // non-quadInterp BC's not implemented at this point
      // so fail
      CH_assert(a_quadInterp);
    }

    // subtract coarse-level data from fluxReg
    const Interval comps(0,0);

    DataIterator dit = a_div.dataIterator();
    // iterate over coarse boxes
    for (dit.reset(); dit.ok(); ++dit)
    {
      FluxBox& thisFlux = uEdge[dit()];
      Real scale = 1.0;
      // iterate over directions
      for (int dir=0; dir<SpaceDim; dir++)
      {
        FArrayBox& thisDirFlux = thisFlux[dir];
        FR.incrementCoarse(thisDirFlux, scale, dit(),
                           comps, comps, dir);
      }
    }

    // now do fine-level data...
    // loop over fine grids, only averaging cell-> edge for the
    // appropriate regions

    DataIterator ditFine = uFine.dataIterator();
    FArrayBox cellData;
    FArrayBox edgeData;
    Real scale = 1.0;

    for (ditFine.reset(); ditFine.ok(); ++ditFine)
    {
      FArrayBox& fineCCFab = uFine[ditFine()];
      for (int dir = 0; dir<SpaceDim; dir++)
      {
        SideIterator sit;
        for (sit.begin(); sit.ok(); sit.next())
        {
          Side::LoHiSide hiorlo = sit();
          Box ccEdgeBox;
          Box edgeBox;

          if (hiorlo == Side::Lo)
          {
            ccEdgeBox = adjCellLo(fineCCFab.box(),dir, 2);
            ccEdgeBox.shift(dir,2);
            edgeBox = bdryLo(fineCCFab.box(),dir,1);
            // shift to account for ghost cell in fineCCFab
            edgeBox.shift(dir,1);
          }
          else
          {
            ccEdgeBox = adjCellHi(fineCCFab.box(),dir,2);
            ccEdgeBox.shift(dir,-2);
            edgeBox = bdryHi(fineCCFab.box(),dir,1);
            // shift to account for ghost cell in fineCCFab
            edgeBox.shift(dir,-1);
          }

          CH_assert(!ccEdgeBox.isEmpty());
          CH_assert(!edgeBox.isEmpty());

          cellData.resize(ccEdgeBox,1);
          edgeData.resize(edgeBox,1);

          // now copy cell-centered fine-level data into cellData
          cellData.copy(fineCCFab,ccEdgeBox,dir,ccEdgeBox,0,1);

          // now need to average cellData->edgeData
          CellToEdge(cellData, 0, edgeData, 0, dir);

          // now increment flux register
          FR.incrementFine(edgeData,scale,ditFine(),comps,
                           comps,dir, hiorlo);

        } // end iteration over high-lo
      } // end iteration over directions
    } // end iteration over fine boxes

    // now perform reflux divergence of edge-vel mismatch
    // this only works if dx == dy
    scale = 1.0/a_dx;

    FR.reflux(a_div, scale);
  } // end correction for finer levels
}
