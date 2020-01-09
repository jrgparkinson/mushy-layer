
#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MushyLayerUtils.H"
#include "ParmParse.H"
#include "CoarseAverage.H"
#include "computeNorm.H"
#include "BCFunc.H"
#include "MushyLayerParams.h"
#include "Gradient.H"
#include "SetValLevel.H"
#include "CellToEdge.H"
#include "CellToEdge2.H"

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>


#include "NamespaceHeader.H"



Real computeSoluteFlux(const LevelData<FArrayBox>& S, const LevelData<FArrayBox>& Sl, const LevelData<FArrayBox>& porosity, const LevelData<FluxBox>& vel, Real a_dx, MushyLayerParams m_parameters, int a_dir, Side::LoHiSide a_side)
{
  //  Real Fs = 0;

  DisjointBoxLayout grids = S.disjointBoxLayout();
  LevelData<FluxBox> flux(grids, 1);
  LevelData<FluxBox> fluxFrameAdv(grids, 1);
  LevelData<FluxBox> fluxFluidAdv(grids, 1);
  LevelData<FluxBox> fluxDiffusion(grids, 1);
  LevelData<FluxBox> porosityEdge(grids, 1);
  LevelData<FArrayBox> shiftedS(grids, 1, IntVect::Unit);
  LevelData<FArrayBox> shiftedSl(grids, 1, IntVect::Unit);
  DataIterator dit = grids.dataIterator();

  setValLevel(fluxDiffusion, 0.0);
  setValLevel(fluxFrameAdv, 0.0);
  setValLevel(flux, 0.0);
  setValLevel(fluxFluidAdv, 0.0);

  Gradient::levelGradientMAC(fluxDiffusion, Sl, a_dx);

  for (dit.reset(); dit.ok(); ++dit)
  {
    // Make sure S is positive, so flux has correct sign (-ve flux = flux out of domain)
    shiftedS[dit].copy(S[dit]);
    //      shiftedS[dit].plus(m_parameters.compositionRatio);
    shiftedSl[dit].copy(Sl[dit]);
    //      shiftedSl[dit].plus(m_parameters.compositionRatio);
  }

  // Do cell to face averaging
  CellToEdge(shiftedS, fluxFrameAdv);
  CellToEdge(S, fluxFluidAdv);
  CellToEdge2(porosity, porosityEdge, geometricAveraging);

  // Compute fluxes
  for (dit.reset(); dit.ok(); ++dit)
  {
    flux[dit].setVal(0.0);

    fluxFrameAdv[dit][0].mult(0.0);
    fluxFrameAdv[dit][1].mult(m_parameters.nonDimVel);


    for (int dir=0; dir<SpaceDim; dir++)
    {
      fluxFluidAdv[dit][dir].mult(vel[dit][dir]);
      fluxDiffusion[dit][dir].mult(porosityEdge[dit][dir], 0, 0, 1);

      fluxDiffusion[dit][dir].mult(m_parameters.m_saltDiffusionCoeff);

      flux[dit][dir].plus(fluxFluidAdv[dit][dir]);
      flux[dit][dir].plus(fluxFrameAdv[dit][dir]);
      flux[dit][dir].minus(fluxDiffusion[dit][dir]);

    }
  }

  Box domBox = grids.physDomain().domainBox();

  // This is the length perpendicular to the direction we want the flux in
  int perpDir = 1;
  if (a_dir == 1)
  {
    perpDir = 0;
  }
  Real length = domBox.hiVect()[a_dir] + 1 - domBox.loVect()[a_dir];
  Real width = (domBox.hiVect()[perpDir] + 1 - domBox.loVect()[perpDir])*a_dx;



  // Get the box containing cells at the required edge of the domain
  /*

  Box ghost;
  if (a_side == Side::Lo)
  {
    ghost = ::adjCellLo(domBox, a_dir, 1);
  }
  else
  {
    ghost = ::adjCellHi(domBox, a_dir, 1);
  }

  Box bottomBox = ghost;

  // shift towards domain
  if (a_side == Side::Lo)
  {
    // Shift up
    bottomBox.shift(a_dir, 1);
  }
  else
  {
    //Shift down
    bottomBox.shift(a_dir, -1);
  }

  // Convert Side:Lo to -1, Side::Hi to 1
  //      int shiftDir = int(a_side);
  //      if (shiftDir == 0)
  //      {
  //        shiftDir = -1;
  //      }
  //      bottomBox.shift(a_dir, -shiftDir); // shift towards domain





  for (dit.reset(); dit.ok(); ++dit)
  {

    // Find the portion of the edgeBox that is covered by this dataIterator
    Box thisBottomBox = bottomBox;
    Box thisBox = flux[dit].box();
    thisBottomBox &= thisBox;

    for (BoxIterator bit = BoxIterator(thisBottomBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      Fs += a_dx*(flux[dit][a_dir](iv));
    }
  }

  Fs = Fs/width;

   */

  // New method - calculate flux at each point in the domain
  Box strip = ::adjCellLo(domBox, a_dir, 1);

  domBox.grow(a_dir, 1);
  //  length = length + 1;

  // This is now the strip of cells at the bottom of the domain
  //  strip.shift(a_dir, 1);

  Vector<Real> fluxes(length+1, 0.0);
  int i=0;

  while(i <= length)
  {
    Real Fs = 0;

    for (dit.reset(); dit.ok(); ++dit)
    {

      // Find the portion of the edgeBox that is covered by this dataIterator
      Box thisStripBox = strip;
      Box thisBox = flux[dit].box();
      thisStripBox &= thisBox;


      for (BoxIterator bit = BoxIterator(thisStripBox); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        //          Fs += a_dx*(flux[dit][a_dir](iv));
      }

      thisStripBox = strip;
      thisStripBox &= thisBox;

      for (BoxIterator bit = BoxIterator(thisStripBox); bit.ok(); ++bit)
      {
        IntVect iv = bit();

        // Can also try calculating by finite diffs
        IntVect ivUp = iv + BASISV(a_dir);

        Real averageSl = (Sl[dit](ivUp) + Sl[dit](iv)) / 2;
        Real averageS = (S[dit](ivUp) + S[dit](iv)) / 2;
        Real w =  vel[dit][a_dir](ivUp);
        //          Real chi = (porosity[dit](iv) + porosity[dit](ivUp))/2;
        Real chi = sqrt(porosity[dit](iv) * porosity[dit](ivUp));

        Real Fframe = m_parameters.nonDimVel * averageS;
        Real Ffluid = w * averageSl;
        Real Fdiffusive = (m_parameters.m_saltDiffusionCoeff) * (chi * (Sl[dit](ivUp) - Sl[dit](iv)) / a_dx);


        Fs += a_dx * (Fframe + Ffluid - Fdiffusive);
      }
    }

    Fs = Fs/width;
    fluxes[i] = Fs;
    i++;
    strip.shift(a_dir, 1);
  }

  Real thisFlux = 0;



  // Another attempt to get the flux - from the operator


  // Choose which flux to return
  if (a_side == Side::Lo)
  {
    thisFlux =  fluxes[1];
  }
  else
  {
    thisFlux = fluxes[length-1];
  }
  int srcProc = 0;
  broadcast(thisFlux, srcProc);

  return thisFlux;

  //  return Fs;
}

Real computeNusselt(const LevelData<FArrayBox>& T, const LevelData<FArrayBox>& vel,
                    Real a_dx, MushyLayerParams a_params, Real a_domWidth, Real a_domHeight)
{
  Real Nu = 0;
  DisjointBoxLayout grids = T.disjointBoxLayout();
  LevelData<FArrayBox> gradT(grids, SpaceDim);
  LevelData<FluxBox> gradTedge(grids, 1);

  // 0 = global average, 1 = sides
  int NusseltType = 1;

  //which direction to compute nu in
  // i.e. x direction means compute dT/dx
  int dir = 0;

  Real sideLength;
  if (dir == 0)
  {
    sideLength = a_domHeight;
  }
  else
  {
    sideLength = a_domWidth;
  }

  Box domBox = grids.physDomain().domainBox();

  // Calculate grad(T)
  if (NusseltType == 1)
  {
    Gradient::levelGradientMAC(gradTedge, T, a_dx);

    Real scale = 0.5*a_dx/sideLength;

    for (SideIterator sit = SideIterator(); sit.ok(); ++sit)
    {

      Side::LoHiSide side = sit();

      Box sideBox;

      if (side == Side::Lo)
      {
        sideBox = ::adjCellLo(domBox, dir, 1);
        sideBox.shift(dir, 1);
      }
      else
      {
        sideBox = ::adjCellHi(domBox, dir, 1);
        sideBox.shift(dir, -1);
      }


      for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
      {

        Box thisBox = grids[dit];

        thisBox &= sideBox;

        for (BoxIterator bit(thisBox); bit.ok(); ++bit)
        {
          IntVect iv = bit();

          Nu += gradTedge[dit][dir](iv)*scale;
        }
      }



    }

//    int temp=0;


  }
  else
  {
    Gradient::levelGradientCC(gradT, T, a_dx);
    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      Box thisBox = grids[dit];
      for (BoxIterator bit(thisBox); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        // Add (V*T + dTdz) * dx to Nusselt
        Nu += (-vel[dit](iv, 1)*T[dit](iv) + gradT[dit](iv, 1))*a_dx*a_dx/(a_domWidth*a_domHeight);
      }
    }
  }


  // Scale with conductive temperature profile
  //  Real height = (domBox.hiVect()[1] + 1 - domBox.loVect()[1])*a_dx;
  Real deltaT = a_params.bcValTemperatureHi[0] - a_params.bcValTemperatureLo[0]; // Just fix this for now
  Nu /= deltaT/a_domHeight;

  // Nu is only calculated by proc 0, apparently
  // need to broadcast to all procs
  int srcProc = 0;
  broadcast(Nu, srcProc);


  return Nu;
}

void
getProblemDomain(ProblemDomain& a_domain)
{
  CH_TIME("getProblemDomain");
  ParmParse pp("main");
  Vector<int> ncell(SpaceDim);
  bool symmetric = false;
  pp.query("symmetric_domain", symmetric);
  pp.getarr("num_cells", ncell, 0, SpaceDim);

  IntVect loEnd, hiEnd;

  if (symmetric)
  {
    hiEnd = IntVect(D_DECL((ncell[0]/2)-1,(ncell[1]/2)-1,(ncell[2]/2)-1));
    loEnd = IntVect(D_DECL(-(ncell[0]/2),-(ncell[1]/2),-(ncell[2]/2)));
  }
  else
  {
    hiEnd = IntVect(D_DECL(ncell[0]-1,ncell[1]-1,ncell[2]-1));
    loEnd = IntVect::Zero;
  }
  //	IntVect hiEnd
  //	Box level0Domain(IntVect::Zero, hiEnd);
  //	IntVect
  //	IntVect loEnd(
  Box level0Domain(loEnd, hiEnd);

  Vector<int> v_is_periodic(SpaceDim);
  pp.getarr("periodic_bc", v_is_periodic, 0, SpaceDim);

  bool is_periodic[SpaceDim];
  for (int idir = 0; idir < SpaceDim; idir++) is_periodic[idir] = (v_is_periodic[idir]==1);

  ProblemDomain prob_domain(level0Domain.smallEnd(),
                            level0Domain.bigEnd(),
                            &(is_periodic[0]));

  a_domain = prob_domain;
}

//(non-dimensional) Calculate the liquidus temperature for a given concentration
//Real
//liquidusTemp(Real C)
//{
//  return C;
//}

void getScalarVarNames(Vector<string>& varNames,
                       const int numScalarVars,
                       const Vector<string> scalarVarNames)
{
  for (int a_var = 0; a_var<numScalarVars; a_var++)
  {
    varNames[a_var] = scalarVarNames[a_var];
  }
}

void getVectorVarNames(Vector<string>& varNames,
                       const int numVectorVars,
                       const Vector<string> vectorVarNames,
                       const int startPos)
{
  for (int a_var = 0; a_var<numVectorVars; a_var++)
  {
    int startComp = startPos + (a_var*2);
    Vector<string> dirString(SpaceDim);
    dirString[0] = string("x");
    if (SpaceDim == 2)
    {
      dirString[1] = string("y");
    }
    else if(SpaceDim == 3)
    {
      dirString[1] = string("y");
      dirString[2] = string("z");
    }

    for (int dir = 0; dir<SpaceDim; dir++)
    {
      //varNames[startComp+dir] = m_vectorVarNames[a_var] + string(" ") + dirString[dir];
      varNames[startComp+dir] = dirString[dir] + vectorVarNames[a_var];
    }
  }
}

void getChkVectorVarNames(Vector<string>& chkVectorVarNames,
                          const Vector<int> chkVectorVars,
                          const Vector<string> vectorVarNames
)
{
  for (int i = 0; i<chkVectorVars.size(); i++)
  {
    int a_var = chkVectorVars[i];

    int startComp = (i*SpaceDim);
    Vector<string> dirString(SpaceDim);
    dirString[0] = string("x");
    if (SpaceDim == 2)
    {
      dirString[1] = string("y");
    }
    else if(SpaceDim == 3)
    {
      dirString[1] = string("y");
      dirString[2] = string("z");
    }

    for (int dir = 0; dir<SpaceDim; dir++)
    {
      //varNames[startComp+dir] = m_vectorVarNames[a_var] + string(" ") + dirString[dir];
      chkVectorVarNames[startComp+dir] = dirString[dir] + vectorVarNames[a_var];
    }
  }
}

void getVarNames(Vector<string>& varNames,
                 const int numScalarVars, const int numVectorVars,
                 const Vector<string> scalarVarNames, const Vector<string> vectorVarNames)
{
  //First do scalar vars
  //  for (int a_var = 0; a_var<numScalarVars; a_var++)
  //  {
  //    varNames[a_var] = scalarVarNames[a_var];
  //  }
  getScalarVarNames(varNames, numScalarVars, scalarVarNames);

  //Then vector vars
  getVectorVarNames(varNames, numVectorVars, vectorVarNames,
                    numScalarVars); // initial offset
  //  Interval vecSrcComps(0, SpaceDim-1);
  //  for (int a_var = 0; a_var<numVectorVars; a_var++)
  //  {
  //    int startComp = numScalarVars + (a_var*2);
  //    Vector<string> dirString(SpaceDim);
  //    dirString[0] = string("x");
  //    if (SpaceDim == 2)
  //    {
  //      dirString[1] = string("y");
  //    }
  //    else if(SpaceDim == 3)
  //    {
  //      dirString[1] = string("y");
  //      dirString[2] = string("z");
  //    }
  //
  //    for (int dir = 0; dir<SpaceDim; dir++)
  //    {
  //      //varNames[startComp+dir] = m_vectorVarNames[a_var] + string(" ") + dirString[dir];
  //      varNames[startComp+dir] = dirString[dir] + vectorVarNames[a_var];
  //    }
  //  }
}

void getVarNames(Vector<string>& varNames,
                 const Vector<int> scalarVars, const Vector<int> vectorVars,
                 const Vector<string> scalarVarNames, const Vector<string> vectorVarNames)
{
  //First do scalar vars
  int numScalarVars = scalarVars.size();
  int numVectorVars = vectorVars.size();

  varNames.resize(numScalarVars + SpaceDim*numVectorVars);

  for (int a_var = 0; a_var<numScalarVars; a_var++)
  {
    varNames[a_var] = scalarVarNames[scalarVars[a_var]];
  }

  //Then vector vars
  Interval vecSrcComps(0, SpaceDim-1);
  for (int a_var = 0; a_var<numVectorVars; a_var++)
  {
    int startComp = numScalarVars + (a_var*SpaceDim);
    Vector<string> dirString(SpaceDim);
    dirString[0] = string("x");
    if (SpaceDim == 2)
    {
      dirString[1] = string("y");
    }
    else if(SpaceDim == 3)
    {
      dirString[1] = string("y");
      dirString[2] = string("z");
    }

    for (int dir = 0; dir<SpaceDim; dir++)
    {
      //varNames[startComp+dir] = m_vectorVarNames[a_var] + string(" ") + dirString[dir];
      varNames[startComp+dir] = dirString[dir] + vectorVarNames[vectorVars[a_var]];
    }
  }
}

void printRepoVersion()
{
  // Get current code revision and write it out so we
  // know exactly what code created the output
  //      system(' hg identify --num');

  const char* cmd = " hg identify --num";
  std::array<char, 128> buffer;
  std::string result;
  std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
  if (!pipe) throw std::runtime_error("popen() failed!");
  while (!feof(pipe.get())) {
    if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
      result += buffer.data();
  }

  pout() << "Mercurial repository revision: " << result << endl;
}


void calculatePermeability(FArrayBox& permeabilityFAB, FArrayBox& solidFractionFAB,
                           MushyLayerParams& params, Real a_dx)
{

  BoxIterator bit(permeabilityFAB.box());
  for (bit.begin(); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    RealVect loc;
    getLocation(iv, loc, a_dx);
    Real x = loc[0];
    Real z = loc[1];

    Real solidFraction = solidFractionFAB(iv,0);

    // Make sure this isn't 0
    //    Real minSolidFractionAllowed = 1e-8;
    //    solidFraction = max(solidFraction, minSolidFractionAllowed);

    Real liquidFraction = (1-solidFraction);
    Real permeability;
//    Real referencePerm = params.referencePermeability;

    if(params.permeabilityFunction == PermeabilityFunctions::m_permeabilityXSquared)
    {
//      permeability = x*x;

      Real scale = 0.1;

      // Creates a channel of high permeability in the middle of the domain
//      permeability = exp(-pow(x-0.5,2)/scale);

      // Create a hole of low permeabilty in the middle of the domain
      //					permeability = 1-exp(-pow(x-0.5,2)/scale)*exp(-pow(z-0.5,2)/scale);

      // Two permeability holes in the left and right of the domain
      scale = 0.03;
      Real zScale = 0.12;
      permeability = exp(-pow(x-0.2,2)/scale)*exp(-pow(z-0.5,2)/zScale);
      permeability = permeability + exp(-pow(x-0.8,2)/scale)*exp(-pow(z-0.5,2)/zScale);
      permeability = 1-permeability;

      //Block flow at the boundaries
      scale = 0.02;
      permeability = exp(-pow(x-0.5,2)/scale);

      if (params.heleShaw)
      {
        permeability = 1 / (1/params.heleShamPermeability + 1.0/permeability);
      }
    }

    else
    {
      //			MayDay::Error("amrMushyLayer::calculatePermeability() - Unknown permeability function");
      permeability = params.calculatePermeability(liquidFraction);
    }


    permeabilityFAB(iv, 0) = permeability;


  } //box iterator
}

void getLocation(const IntVect iv, RealVect& loc, const Real a_dx, const RealVect ccOffset)
{
  //Default value for ccOffsetScale is 0.5 i.e. cell centred.

  // I think this was wrong before? As ccoffset was already scaled by m_dx
  //  RealVect scale = RealVect(ccOffsetX, ccOffsetY);
  //  RealVect ccOffset = scale*a_dx*RealVect::Unit;

  RealVect scaledOffset = ccOffset*a_dx;

  loc = iv;
  loc *= a_dx;
  loc += scaledOffset;
}


Real burgersPeriodicInit(Real x, Real y, int dir, MushyLayerParams params)
{


  if (dir == 1)
  {
    return 0.0;
  } else {



    Real phi = exp(-x*x/(4*params.darcy)) + exp(-pow(x-2*M_PI, 2)/(4*params.darcy));
    Real dPhidx = (-1/(2*params.darcy))* (x*exp(-x*x/(4*params.darcy)) + (x-2*M_PI)*exp(-pow(x-2*M_PI, 2)/(4*params.darcy)));
    Real val = 4 - 2*params.darcy * dPhidx/phi;

    return val;
  }
}

Real burgersSinInit(Real x, Real y, int dir, MushyLayerParams params)
{
  if (dir == 1)
  {
    return 0.0;
  } else {

    Real val = - sin(M_PI*x);

    return val;
  }
}

/*
 *
 * Benchmark solution for Stokes flow with darcy term:
 * grad(P) = 1
 * porosity = 0.5, karmen-kozeny permeability
 * 0<x<1
 * periodic in y
 * See latex file for more details
 */
Real stokesDarcyInit(Real x, Real y, int dir, MushyLayerParams params)
{

  if (dir == 1)
  {

    Real val = 1 - cosh(sqrt(0.5)*(0.5-x))/cosh(0.5*sqrt(0.5));
    return val;


  } else {
    return 0.0;
  }

}


#include "NamespaceFooter.H"

