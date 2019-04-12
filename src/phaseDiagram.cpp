#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "phaseDiagram.H"
#include "ParmParse.H"
#include "CoarseAverage.H"
#include "computeNorm.H"
#include "BCFunc.H"
#include "MushyLayerParams.h"
#include "Gradient.H"
#include "SetValLevel.H"
#include "CellToEdge.H"
#include "EnthalpyVariablesF_F.H"

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>


#include "NamespaceHeader.H"

void updateEnthalpyVariables(LevelData<FArrayBox>& HC,
                             LevelData<FArrayBox>& temperature, LevelData<FArrayBox>& compositionLiquid,
                             LevelData<FArrayBox>& compositionSolid, LevelData<FArrayBox>& porosity,
                             LevelData<FArrayBox>& enthalpySolid, LevelData<FArrayBox>& enthalpyLiquid,
                             LevelData<FArrayBox>& enthalpyEutectic,
                             MushyLayerParams a_params)
{
    DisjointBoxLayout grids = HC.disjointBoxLayout();
    IntVect ghostVect = HC.ghostVect();

    LevelData<FArrayBox> enthalpy(grids, 1, ghostVect);
    LevelData<FArrayBox> composition(grids, 1, ghostVect);

    HC.copyTo(Interval(0,0), enthalpy, Interval(0,0));
    HC.copyTo(Interval(1,1), composition, Interval(0,0));

    updateEnthalpyVariables(enthalpy, composition, temperature, compositionLiquid, compositionSolid,
                              porosity, enthalpySolid,enthalpyLiquid, enthalpyEutectic, a_params);
}

void updateEnthalpyVariables(LevelData<FArrayBox>& HC,
                             LevelData<FArrayBox>& theta, LevelData<FArrayBox>& compositionLiquid,
                             LevelData<FArrayBox>& compositionSolid, LevelData<FArrayBox>& porosity,
                             MushyLayerParams a_params)
{
  DisjointBoxLayout grids = theta.disjointBoxLayout();
  IntVect ghostVect = theta.ghostVect();

  LevelData<FArrayBox> enthalpy(grids, 1, ghostVect);
  LevelData<FArrayBox> composition(grids, 1, ghostVect);

  HC.copyTo(Interval(0,0), enthalpy, Interval(0,0));
  HC.copyTo(Interval(1,1), composition, Interval(0,0));

  updateEnthalpyVariables(enthalpy, composition, theta, compositionLiquid, compositionSolid,
                            porosity, a_params);
}


void updateEnthalpyVariables(LevelData<FArrayBox>& enthalpy, LevelData<FArrayBox>& composition,
                             LevelData<FArrayBox>& theta, LevelData<FArrayBox>& compositionLiquid,
                             LevelData<FArrayBox>& compositionSolid, LevelData<FArrayBox>& porosity,
                             MushyLayerParams a_params)
{
  DisjointBoxLayout grids = theta.disjointBoxLayout();
  IntVect ghostVect = theta.ghostVect();

  LevelData<FArrayBox> enthalpySolid(grids, 1, ghostVect);
  LevelData<FArrayBox> enthalpyLiquid(grids, 1, ghostVect);
  LevelData<FArrayBox> enthalpyEutectic(grids, 1, ghostVect);

  updateEnthalpyVariables(enthalpy, composition, theta, compositionLiquid, compositionSolid,
                          porosity, enthalpySolid, enthalpyLiquid, enthalpyEutectic, a_params);
}

void updateEnthalpyVariables(LevelData<FArrayBox>& enthalpy, LevelData<FArrayBox>& composition,
                             LevelData<FArrayBox>& temperature, LevelData<FArrayBox>& compositionLiquid,
                             LevelData<FArrayBox>& compositionSolid, LevelData<FArrayBox>& porosity,
                             LevelData<FArrayBox>& enthalpySolid, LevelData<FArrayBox>& enthalpyLiquid,
                             LevelData<FArrayBox>& enthalpyEutectic,
                             MushyLayerParams a_params)
{
  DisjointBoxLayout grids = temperature.disjointBoxLayout();
  IntVect ghostVect = temperature.ghostVect();

//  Real NaN = 0; // Sometimes we need to set a variable to 'nan'

  LevelData<FArrayBox> tempSolidus(grids, 1, ghostVect);
  LevelData<FArrayBox> porosityEutectic(grids, 1, ghostVect);
  LevelData<FArrayBox> enthalpyComposition(grids, 2, ghostVect);

  enthalpy.copyTo(Interval(0,0), enthalpyComposition, Interval(0,0));
  composition.copyTo(Interval(0,0), enthalpyComposition, Interval(1,1));

  //Iterate over the domain on this level

  for (DataIterator dit = temperature.dataIterator(); dit.ok(); ++dit)
  {

    FArrayBox& Hs = (enthalpySolid)[dit];
    FArrayBox& He = (enthalpyEutectic)[dit];
    FArrayBox& Hl = (enthalpyLiquid)[dit];

    FArrayBox& H = (enthalpy)[dit];
    FArrayBox& C = (composition)[dit];
    FArrayBox& T = (temperature)[dit];
    FArrayBox& Cl = (compositionLiquid)[dit];
    FArrayBox& Cs = (compositionSolid)[dit];
    FArrayBox& chi = (porosity)[dit];

    FArrayBox& HC = enthalpyComposition[dit];

    Box region = HC.box();
    region &= Hs.box();

    FORT_CALCULATE_BOUNDING_ENERGY( CHF_CONST_FRA(HC),
                                    CHF_FRA(Hs),
                                    CHF_FRA(He),
                                    CHF_FRA(Hl),
                                    CHF_BOX(region),
                                    CHF_CONST_REAL(a_params.compositionRatio),
                                    CHF_CONST_REAL(a_params.waterDistributionCoeff),
                                    CHF_CONST_REAL(a_params.specificHeatRatio),
                                    CHF_CONST_REAL(a_params.stefan),
                                    CHF_CONST_REAL(a_params.thetaEutectic),
                                    CHF_CONST_REAL(a_params.ThetaEutectic));

    //Now we calculate porosity, temperature and liquid/solid composition at each grid point based on what enthalpy regime we are in
    FORT_CALCULATEPOROSITY( CHF_FRA(chi),
                             CHF_CONST_FRA(H),
                             CHF_CONST_FRA(C),
                             CHF_CONST_FRA(Hs),
                             CHF_CONST_FRA(He),
                             CHF_CONST_FRA(Hl),
                             CHF_BOX(region),
                             CHF_CONST_REAL(a_params.compositionRatio),
                             CHF_CONST_REAL(a_params.waterDistributionCoeff),
                             CHF_CONST_REAL(a_params.specificHeatRatio),
                             CHF_CONST_REAL(a_params.stefan),
                             CHF_CONST_REAL(a_params.thetaEutectic));

    FORT_CALCULATECL( CHF_FRA(Cl),
                                 CHF_CONST_FRA(H),
                                 CHF_CONST_FRA(C),
                                 CHF_CONST_FRA(Hs),
                                 CHF_CONST_FRA(He),
                                 CHF_CONST_FRA(Hl),
                                 CHF_BOX(region),
                                 CHF_CONST_REAL(a_params.compositionRatio),
                                 CHF_CONST_REAL(a_params.waterDistributionCoeff),
                                 CHF_CONST_REAL(a_params.specificHeatRatio),
                                 CHF_CONST_REAL(a_params.stefan),
                                 CHF_CONST_REAL(a_params.ThetaEutectic));

    FORT_CALCULATECS( CHF_FRA(Cs),
                                 CHF_CONST_FRA(H),
                                 CHF_CONST_FRA(C),
                                 CHF_CONST_FRA(Hs),
                                 CHF_CONST_FRA(He),
                                 CHF_CONST_FRA(Hl),
                                 CHF_BOX(region),
                                 CHF_CONST_REAL(a_params.compositionRatio),
                                 CHF_CONST_REAL(a_params.waterDistributionCoeff),
                                 CHF_CONST_REAL(a_params.specificHeatRatio),
                                 CHF_CONST_REAL(a_params.stefan));

    FORT_CALCULATET( CHF_FRA(T),
                                 CHF_CONST_FRA(H),
                                 CHF_CONST_FRA(C),
                                 CHF_CONST_FRA(Hs),
                                 CHF_CONST_FRA(He),
                                 CHF_CONST_FRA(Hl),
                                 CHF_BOX(region),
                                 CHF_CONST_REAL(a_params.compositionRatio),
                                 CHF_CONST_REAL(a_params.waterDistributionCoeff),
                                 CHF_CONST_REAL(a_params.specificHeatRatio),
                                 CHF_CONST_REAL(a_params.stefan),
                                 CHF_CONST_REAL(a_params.thetaEutectic));



//    int tempVar = 0;

  } // end loop over boxes

}

// Need this method so we can calculate boundary conditions
void computeEnthalpyVars(const Real H, const Real C, Real& porosity, Real& theta, Real& C_l, Real& C_s,
                         const Real H_s, const Real H_l, const Real H_e,
                         const Real specificHeatRatio, const Real stefan, Real const compositionRatio,
                         const Real waterDistributionCoeff, const Real thetaEutectic, const Real ThetaEutectic)
{
  if (H <= H_s)
  {
    C_l = 0;
    C_s = C;
    theta = H/specificHeatRatio;
    porosity = 0.0;
  }
  else if (H > H_s && H <= H_e)
  {
    porosity = (H-thetaEutectic*specificHeatRatio)/(stefan + thetaEutectic*(1-specificHeatRatio));
    C_l = ThetaEutectic;
    C_s = C/(1-porosity);
    theta = thetaEutectic;
  }
  else if (H > H_e && H < H_l)
  {


    porosity = computePorosityMushyLayer( H,  C,  compositionRatio,  specificHeatRatio,
                                          stefan,  waterDistributionCoeff);


    C_l = (C + compositionRatio*(1-porosity)) / (porosity + waterDistributionCoeff*(1-porosity));
    C_s = (waterDistributionCoeff*C - compositionRatio*porosity) / (porosity + waterDistributionCoeff*(1-porosity));
    theta = -C_l;
  }
  else
  {
    C_l = C;
    porosity = 1.0;
    theta = H - stefan;
    C_s = 0;
  }
}


Real computePorosity(Real H, Real C, Real compositionRatio, Real specificHeatRatio,
                               Real stefan, Real waterDistributionCoeff, Real heatCapacityRatio,
                               Real thetaEutectic, Real ThetaEutectic)
{
  Real H_e, H_s, H_l, porosity;

  ::computeBoundingEnergy(H_e, C, H_s, H_l, H_e, heatCapacityRatio, stefan, compositionRatio, waterDistributionCoeff, thetaEutectic, ThetaEutectic);

  if (H <= H_s)
  {

    porosity = 0.0;
  }
  else if (H > H_s && H <= H_e)
  {
    porosity = (H-thetaEutectic*specificHeatRatio)/(stefan + thetaEutectic*(1-specificHeatRatio));

  }
  else if (H > H_e && H < H_l)
  {


    porosity = computePorosityMushyLayer( H,  C,  compositionRatio,  specificHeatRatio,
                                          stefan,  waterDistributionCoeff);



  }
  else
  {

    porosity = 1.0;
  }

  return porosity;
}


// Refactored this out so I can reuse it elsewhere
Real computePorosityMushyLayer(Real H, Real a_C, Real compositionRatio, Real specificHeatRatio,
                               Real stefan, Real waterDistributionCoeff)
{
  CH_TIME("MushyLayerUtils::computePorosityMushyLayer");

  Real porosity;

  if (stefan == 0)
  {
    porosity = 1; //Code can't handle this
  }
  else
  {
    Real a = compositionRatio*(specificHeatRatio -1) + stefan * (waterDistributionCoeff-1);
    Real b = compositionRatio * (1-2*specificHeatRatio) + H*(1-waterDistributionCoeff)
                                        - a_C * (specificHeatRatio-1) - waterDistributionCoeff * stefan;
    Real c = (a_C + compositionRatio)*specificHeatRatio +
        waterDistributionCoeff * H;

    porosity = (-b -sqrt(b*b - 4*a*c))/(2*a);
  }

  return porosity;
}



void computeBoundingEnergy(const Real H, const Real C, Real& H_s, Real& H_l, Real& H_e,
                           const Real heatCapacityRatio, const Real stefan,
                           const Real compositionRatio, const Real waterDistributionCoeff,
                           const Real thetaEutectic, const Real ThetaEutectic)
{
  Real zero = 0.0;
  H_s = heatCapacityRatio*(thetaEutectic + max(zero, (-C-compositionRatio)/waterDistributionCoeff));

  Real chi_e = (compositionRatio+C)/(ThetaEutectic + compositionRatio);
  H_e = chi_e*(stefan + thetaEutectic*(1- heatCapacityRatio)) + heatCapacityRatio*thetaEutectic;
  H_l = stefan - C + thetaEutectic + ThetaEutectic;
}


#include "NamespaceFooter.H"
