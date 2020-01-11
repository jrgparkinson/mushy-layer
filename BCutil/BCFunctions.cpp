#include <cassert>
#include "DataIterator.H"
#include "BCFunctions.H"
#include "BCFuncF_F.H"
#include "PhysBCUtil.H"

#include "CellToEdgeF_F.H"
#include "CoarseAverageFace.H"
#include "CoarseAverage.H"
#include "NamespaceHeader.H"


void applyCorrectScalarBC(FArrayBox&      a_state,
                          LevelData<FluxBox>* a_advVel,
                          const Box&      a_valid,
                          bool            a_homogeneous,
                          Real            a_diriValue,
                          Real            a_inflowDiriValue,
                          Real            a_neumValue,
                          int             a_idir,
                          Side::LoHiSide  a_side,
                          Real            a_dx,
                          int             a_order,
                          int             a_bcType,
                          Vector<Real>   a_plumeBounds,
                          Real a_fineDx,
                          int a_comp)
{
  CH_TIME("BCFunctions::applyCorrectScalarBC");

  if (a_bcType == PhysBCUtil::Dirichlet)
  {
    ConstantDiriBC(a_state,
                   a_valid,
                   a_homogeneous,
                   a_diriValue,
                   a_idir,
                   a_side,
                   a_order,
                   a_comp);
  }
  else if (a_bcType == PhysBCUtil::InflowOutflow)
  {
    InflowOutflowBC(a_state,
                    a_advVel,
                    a_valid,
                    a_homogeneous,
                    a_diriValue,
                    a_idir,
                    a_side,
                    a_order,
                    a_dx,
                    a_fineDx,
                    a_comp);
  }
  else if (a_bcType == PhysBCUtil::Neumann)
  {
    ConstantNeumBC(a_state,
                   a_valid,
                   a_homogeneous,
                   a_neumValue,
                   a_idir,
                   a_side,
                   a_dx,
                   a_comp);
  }
  else if (a_bcType == PhysBCUtil::OnlyInflow)
  {
    OnlyInflowBC(a_state,
                 a_valid,
                 a_homogeneous,
                 a_inflowDiriValue,
                 a_diriValue,
                 a_plumeBounds,
                 a_dx,
                 a_idir,
                 a_side,
                 a_order,
                 a_comp);

  }
  else if (a_bcType == PhysBCUtil::VariableFlux)
    {
    VariableFluxBC(a_state,
                       a_valid,
                       a_homogeneous,
                       a_diriValue,
                       a_idir,
                       a_side,
                       a_dx,
                       a_comp);
  }
  else
  {
    MayDay::Error("Unknown BCs");
  }
}

// For when each side has a fixed value along the entire side
// Currently just for the first component of a_state
void ConstantDiriBC(FArrayBox&      a_state,
                    const Box&      a_valid,
                    bool            a_homogeneous,
                    Real   			a_value,
                    int             a_dir,
                    Side::LoHiSide  a_side,
                    int             a_order,
                    int a_comp)
{
  CH_TIME("BCFunctions::ConstantDiriBC");
  int isign = sign(a_side);

  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
  toRegion &= a_state.box();

  if (a_homogeneous)
  {
    a_value = 0;
  }


  FORT_DIRIBC(CHF_FRA(a_state),
              CHF_REAL(a_value),
              CHF_BOX(toRegion),
              CHF_INT(a_order),
              CHF_INT(isign),
              CHF_INT(a_dir),
              CHF_INT(a_comp));

}

// For when each side has a fixed gradient along the entire side
// Currently just for the first component of a_state
void ConstantNeumBC(FArrayBox&      a_state,
                    const Box&      a_valid,
                    bool            a_homogeneous,
                    Real   			a_value,
                    int             a_dir,
                    Side::LoHiSide  a_side,
                    Real			a_dx,
                    int a_comp)
{
  CH_TIME("BCFunctions::ConstantNeumBC");
  int isign = sign(a_side);

  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
  toRegion &= a_state.box();

  if (a_homogeneous)
  {
    a_value = 0;
  }


  FORT_NEUMBC(CHF_FRA(a_state),
              CHF_REAL(a_value),
              CHF_BOX(toRegion),
              CHF_INT(isign),
              CHF_INT(a_dir),
              CHF_REAL(a_dx),
              CHF_INT(a_comp));

}

// For when each side has a fixed gradient along the entire side
// Currently just for the first component of a_state
void VariableFluxBC(FArrayBox&      a_state,
                    const Box&      a_valid,
                    bool            a_homogeneous,
                    Real                        a_value,
                    int             a_dir,
                    Side::LoHiSide  a_side,
                    Real                        a_dx,
                    int a_comp)
{
  CH_TIME("BCFunctions::VariableFluxBC");
  int isign = sign(a_side);

  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
  toRegion &= a_state.box();

  if (a_homogeneous)
  {
    a_value = 0;
  }

//  FORT_NEUMBC(CHF_FRA(a_state),
//              CHF_REAL(a_value),
//              CHF_BOX(toRegion),
//              CHF_INT(isign),
//              CHF_INT(a_dir),
//              CHF_REAL(a_dx),
//              CHF_INT(a_comp));

  RealVect loc;
//  pout() << endl;

  for (BoxIterator bit = BoxIterator(toRegion); bit.ok(); ++bit)
  {
    IntVect iv = bit();
//    ::getLocation(iv, loc, a_dx);

    loc = iv;
    loc *= a_dx;
//    RealVect ccOffset = (0.5, 0.5);
    RealVect ccOffset = 0.5*RealVect::Unit;
    ccOffset *= a_dx;
    loc += ccOffset; // cell centred offset

    Real flux = a_value*0.5*(1+tanh(50*(loc[SpaceDim-1]-0.75)));

//    pout() << iv[1] << ",";

//    if (loc[1] > 0.5)
//    {
//      flux = a_value;
//    }
//    else
//    {
//      flux = 0.0;
//    }

//    flux = 200;

    IntVect ivFrom = iv - isign*BASISV(a_dir);

    a_state(iv, a_comp) = a_state(ivFrom, a_comp) + a_dx*flux;


  }

//  pout() << endl;

}


void ExtrapBC(FArrayBox&      a_state,
              const Box&      a_valid,
              int             a_dir,
              Side::LoHiSide  a_side,
              int             a_order,
              int a_comp)
{
  CH_TIME("BCFunctions::ExtrapBC");
  int isign = sign(a_side);

  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
  Box thisBox = a_state.box();
  IndexType::CellIndex ix = thisBox.type(a_dir);

//  int perpDir = (a_dir == 0) ? 1 : 0;

  if (ix == IndexType::NODE)
  {
    Box enclosedCells = thisBox.enclosedCells(a_dir);
    toRegion &= enclosedCells;
  }
  // Should check both perpendicular directions
  else
  {
    for (int dir=0; dir<SpaceDim; dir++)
    {
      // Skip this direction
      if (dir == a_dir)
      {
        continue;
      }

      // This must be a perpendicular dir to get this far
      if (thisBox.type(dir) == IndexType::NODE)
      {
        Box enclosedCells = thisBox.enclosedCells(dir);
        toRegion &= enclosedCells;
      }
    }
  }

  // If user asks for second order accuracy but we only have one interior point then don't do what they ask for
  int boxSize = a_state.box().size(a_dir);
  a_order = min(a_order, boxSize);

  CH_assert(a_comp < a_state.nComp());

  FORT_EXTRAPBC(CHF_FRA(a_state),
                CHF_BOX(toRegion),
                CHF_INT(isign),
                CHF_INT(a_dir),
                CHF_INT(a_order),
                CHF_INT(a_comp));

}

void InflowOutflowBC(FArrayBox&      a_state,
                     const LevelData<FluxBox>*      a_advVel,
                     const Box&      a_valid,
                     bool            a_homogeneous,
                     Real                        a_DiriValue,
                     int             a_dir,
                     Side::LoHiSide  a_side,
                     int             a_order,
                     Real a_dx,
                     Real fine_dx,
                     int a_comp)
{
  CH_TIME("BCFunctions::InflowOutflowBC");
  int isign = sign(a_side);
  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
  toRegion &= a_state.box();

  // If no velocity, then we should throw an error
  // however just apply diri BC
  if (a_advVel == nullptr)
  {
    pout() << "InflowOutflowBC error - haven't specified the velocity field";
    //    ConstantDiriBC(a_state, a_valid, a_homogeneous, a_DiriValue, a_dir, a_side, a_order);
    //    return;
  }


  if (a_homogeneous)
  {
    a_DiriValue = 0;
  }

  // adv vel is never going to match in this case, so don't bother searching
  if (a_dx > fine_dx)
  {
    ConstantDiriBC(a_state, a_valid, a_homogeneous, a_DiriValue, a_dir, a_side, a_order, a_comp);
    return;
  }

  FArrayBox advVelDir;

  // Don't currently do this
//  bool tryCoarsening = false;
//  if (tryCoarsening)
//  {
//    int refRat = a_dx/fine_dx;
//
//    DataIterator dit = a_advVel->dataIterator();
//    //    pout() << "valid box: " << a_valid << endl;
//
//    for (dit.reset(); dit.ok(); ++dit)
//    {
//      const Box& advVelBox = (*a_advVel)[dit].box();
//      Box coarsenableBox(advVelBox);
//      coarsenableBox.coarsen(refRat);
//      //      pout() << "advVeBox: " << advVelBox << endl;
//
//      if (coarsenableBox.contains(a_valid))
//      {
//        //        pout() << "advVel contains valid" << endl;
//        break;
//      }
//    }

    advVelDir.define(a_valid, 1);
    advVelDir.setVal(0.0);

    CoarseAverage average;
    for (BoxIterator bit((*a_advVel)[dit].box()); bit.ok(); ++bit)
    {
      IntVect ivFine = bit();
      IntVect ivCoarse = ivFine;
      ivCoarse.coarsen(refRat);
      if (a_valid.contains(ivCoarse))
      {
        advVelDir(ivCoarse) += 0.5*(*a_advVel)[dit][a_dir](ivFine);
      }
    }

  }

  else
  {
    CH_TIME("InflowOutflowBC::findAdvVel");

    DataIterator dit = a_advVel->dataIterator();

    const Box& stateBox = a_state.box();

    for (dit.reset(); dit.ok(); ++dit)
    {
      const FluxBox& flux = (*a_advVel)[dit];
      Box advVelBox = flux.box();

      // Adv vel box might have more ghost cells than a_state.box - so remove them for comparison
      advVelBox &= a_state.box();

      // If this advection velocity box contains information relating to this state, apply BCs
      if (advVelBox.contains(stateBox))
      {
        //        matched = true;

        // Extrapolation for outflow, dirichlet for inflow

        FORT_INFLOWOUTFLOWBC(CHF_FRA(a_state),
                             CHF_CONST_FRA(flux[a_dir]),
                             CHF_REAL(a_DiriValue),
                             CHF_BOX(toRegion),
                             CHF_INT(a_order),
                             CHF_INT(isign),
                             CHF_INT(a_dir),
                             CHF_INT(a_comp));
        return;
      }
    }
  } // end if not trying coarsening

  // If we made it this far, we didn't find a data index
  // Can't enforce inflow/outflow, just do dirichlet BCs
  /* There is an issue with the current inflow/outflow BCs - they don't work with multigrid as
   * advVel is defined on the finest level. We should coarsen advVel to match grid,
   * or could store a vector identifying which regions are inflow and which are
   * outflow by their position along the face.
   *
   * However, we seem to survive just doing dirichlet BCs on the coarse grid then enforcing inflow/outflow
   * when we get back to the finest grid
   *
   */

  ConstantDiriBC(a_state, a_valid, a_homogeneous, a_DiriValue, a_dir, a_side, a_order, a_comp);
  //    ConstantNeumBC(a_state, a_valid, a_homogeneous, 0.0, a_dir, a_side, a_dx);
  return;

}

void OnlyInflowBC(FArrayBox&      a_state,
                  const Box&      a_valid,
                  bool            a_homogeneous,
                  Real                        a_inflowValue,
                  Real                        a_otherValue,
                  Vector<Real> a_plumeBounds,
                  Real a_dx,
                  int             a_dir,
                  Side::LoHiSide  a_side,
                  int             a_order,
                  int a_comp)
{

  CH_TIME("BCFunctions::OnlyInflowBC");

  // if homogeneous, apply homogeneous dirichlet BCs
  if (a_homogeneous)
  {
    ConstantDiriBC(a_state, a_valid, a_homogeneous, a_otherValue, a_dir, a_side, a_order, a_comp);
    return;
  }

  int isign = sign(a_side);
  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
  toRegion &= a_state.box();

  FORT_INFLOWBC(CHF_FRA(a_state),
                CHF_REAL(a_inflowValue),
                CHF_REAL(a_otherValue),
                CHF_REAL(a_plumeBounds[0]),
                CHF_REAL(a_plumeBounds[1]),
                CHF_REAL(a_dx),
                CHF_BOX(toRegion),
                CHF_INT(a_order),
                CHF_INT(isign),
                CHF_INT(a_dir),
                CHF_INT(a_comp));


}

void PressureInflowOutflow(FArrayBox&      a_state,
                           const Box&      a_valid,
                           Real            a_dx,
                           bool            a_homogeneous,
                           LevelData<FluxBox>* a_advVel,
                           Real   a_diriValue,
                           Real  a_neumValue,
                           int             a_dir,
                           Side::LoHiSide  a_side,
                           int             a_order,
                           Real a_fineDx,
                           int a_comp)
{
  //  Interval stateInterval = a_state.interval();
  //    DiriBC(a_state,a_valid, a_dx, a_homogeneous,
  //           a_diriValue,a_dir, a_side, stateInterval, a_order);
  if (a_neumValue == 0)
  {
    InflowOutflowBC(a_state, a_advVel, a_valid, a_homogeneous, a_diriValue,
                    a_dir, a_side, a_order, a_dx, a_fineDx, a_comp);
  }
  else
  {
    MayDay::Error("PressureInflowOutflow - can't enforce bcs for non-zero neumann value");
  }
}


void PressurePlumeInflow(FArrayBox&      a_state,
                         const Box&      a_valid,
                         Real            a_dx,
                         bool            a_homogeneous,
                         Vector<Real> a_plumeBounds,
                         Real   a_diriValue,
                         Real  a_neumValue,
                         int             a_dir,
                         Side::LoHiSide  a_side,
                         int             a_order,
                         int a_comp)
{
  int isign = sign(a_side);
  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
  toRegion &= a_state.box();

  FORT_PRESSUREPLUMEBC(CHF_FRA(a_state),
                       CHF_REAL(a_neumValue),
                       CHF_REAL(a_diriValue),
                       CHF_REAL(a_plumeBounds[0]),
                       CHF_REAL(a_plumeBounds[1]),
                       CHF_REAL(a_dx),
                       CHF_BOX(toRegion),
                       CHF_INT(a_order),
                       CHF_INT(isign),
                       CHF_INT(a_dir),
                       CHF_INT(a_comp));
}
