/*
 * MushyLayerParams.cpp
 *
 *  Created on: 23 May 2016
 *      Author: parkinsonjl
 */

#include "MushyLayerParams.h"
#include "MushyLayerUtils.H"
#include "PhysBCUtil.H"
#include "phaseDiagram.H"

//This is for printing out variable names, because c++ doesn't do reflection
#define SHOW(a) pout() << #a << ": " << (a) << std::endl

MushyLayerParams::MushyLayerParams() {
  //Initialise everything

  // Use some ridiculous defaults to highlight any parameters we forget to set

  physicalProblem = -999;
  viscosity = -999 ;
  heatConductivityLiquid = -999;
  heatConductivitySolid = -999 ;
  specificHeatLiquid = -999 ;
  specificHeatSolid = -999 ;
  liquidDensity = -999 ;
  latentHeatDissolution = -999 ;
  thermalExpansivity = -999 ;
  solutalExpansivity = -999 ;
  eutecticTemp = -999 ;
  //  bottomTemp = -999 ;
  //  topTemp = -999 ;
  eutecticComposition = -999 ;
  initialComposition = -999 ;
  liquidusSlope = -999 ;
  waterDistributionCoeff = -999 ;
  heleShawCoolingCoeff = -999 ;
  liquidSoluteDiffusivity = -999 ;
  d = -999 ;
  height = -999 ;
  referencePermeability = -999 ;
  gravitationalAcceleration = -999 ;
  V = -999 ;

  deltaTemp = -999 ;
  deltaSalt = -999;
  stefan = -999 ;
  compositionRatio = -999 ;
  liquidHeatDiffusivity = -999 ;
  heatConductivityRatio = -999 ;
  specificHeatRatio = -999 ;
  lewis = -999 ;
  darcy = -999;
  reynolds = -999;
  rayleighTemp = -999 ;
  rayleighComposition = -999 ;
  nonDimVel = -999 ;
  nonDimHeleShawCooling = -999 ;
  timescale = -999 ;

  thetaEutectic = -999 ;
  thetaInf = -999 ;
  thetaInitialLiquidus = -999 ;
  thetaInitial = -999 ;
  thetaInterface = -999 ;
  ThetaEutectic = -999 ;
  ThetaInitial = -999 ;
  ThetaInf = -999 ;
  ThetaLInitial = -999 ;
  ThetaSInitial = -999;
  Hinitial = -999 ;


  fixedTempDirection = -999;
  heleShaw = false;
  permeabilityFunction = m_kozenyCarman;
  m_porosityFunction = -1;
  inflowVelocity = -999;


  thetaPlumeInflow = -999;
  HPlumeInflow = -999;
  ThetaPlumeInflow = -999;
  ThetaLPlumeInflow = -999;
  porosityPlume = -999;
  ThetaSPlumeInflow = -999;
  permeabilityPlume = -999;
  HSolidusPlume = -999;
  HLiquidusPlume = -999;
  HEutecticPlume = -999;

  nonDimReluctance = -999.0;
  width = -999.0;


  referenceTemperature = -999;
  referenceSalinity = -999;

  m_BCamplitude = 0;
  m_BCtimescale = 1;
  m_time = -999;
  m_timeDependentBC = m_constant;

  prandtl = 0.0;

  m_advectionCoeff = 1.0;
  m_buoyancyTCoeff = 1.0;
  m_nondimensionalisation = 0;
 m_heatDiffusionCoeff= 1;
    m_saltDiffusionCoeff = 1;
           m_viscosityCoeff = 1;
           m_buoyancyTCoeff = 1;
           m_buoyancySCoeff = 1;
           m_darcyCoeff = 1;
           m_advectionCoeff = 1;

}

MushyLayerParams::~MushyLayerParams() {
}


void MushyLayerParams::getParameters()
{
  CH_TIME("MushyLayerParams::getParameters()");

  ParmParse pp("parameters");

  pp.get("problem_type", physicalProblem);

  pp.query("permeabilityFunction", permeabilityFunction);
  pp.query("porosityFunction", m_porosityFunction);
  pp.query("heleShaw", heleShaw);

  // I don't think I necessarily need these, if I'm specifying non dimensional parameters
  pp.query("viscosity", viscosity);
  pp.query("heatConductivityLiquid", heatConductivityLiquid);
  pp.query("heatConductivitySolid", heatConductivitySolid);
  pp.query("specificHeatLiquid", specificHeatLiquid);
  pp.query("specificHeatSolid", specificHeatSolid);
  pp.query("liquidDensity", liquidDensity);
  pp.query("latentHeatDissolution", latentHeatDissolution);
  pp.query("thermalExpansivity", thermalExpansivity);
  pp.query("solutalExpansivity", solutalExpansivity);

  pp.query("heleShawCoolingCoeff", heleShawCoolingCoeff);
  pp.query("liquidSoluteDiffusivity", liquidSoluteDiffusivity);
  pp.query("d", d);
  pp.query("height", height);
  pp.query("referencePermeability", referencePermeability);
  pp.query("gravitationalAcceleration", gravitationalAcceleration);
  pp.query("V", V);

  // I do need these, to define the phase diagram
  pp.get("eutecticTemp", eutecticTemp);
  pp.get("eutecticComposition", eutecticComposition);
  pp.get("initialComposition", initialComposition);
  pp.get("liquidusSlope", liquidusSlope);
  pp.get("waterDistributionCoeff", waterDistributionCoeff);

  pp.query("fixedTempDirection", fixedTempDirection);
  pp.query("inflowVelocity", inflowVelocity);

  pp.query("timeDependentBC", m_timeDependentBC);

  //Derived parameters. Can enforce these if needed (e.g. for benchmarking with a reduced model)
  if (pp.contains("deltaSalt"))
  {
    pp.get("deltaSalt", deltaSalt);
  }
  else
  {
    deltaSalt = (eutecticComposition - initialComposition);
  }


  if (pp.contains("deltaTemp"))
  {
    pp.get("deltaTemp", deltaTemp);
  }
  else
  {
    deltaTemp = - liquidusSlope * deltaSalt;
  }

  if (pp.contains("stefan"))
  {
    pp.get("stefan", stefan);
  }
  else
  {
    stefan = latentHeatDissolution / (specificHeatLiquid * deltaTemp);
  }

  if (pp.contains("compositionRatio"))
  {
    pp.get("compositionRatio", compositionRatio);
  }
  else
  {
    compositionRatio = (1-waterDistributionCoeff) * eutecticComposition / deltaSalt;
  }

  if (pp.contains("liquidHeatDiffusivity"))
  {
    pp.get("liquidHeatDiffusivity", liquidHeatDiffusivity);
  }
  else
  {
    liquidHeatDiffusivity = heatConductivityLiquid/(liquidDensity * specificHeatLiquid);
  }

  if (pp.contains("heatConductivityRatio"))
  {
    pp.get("heatConductivityRatio", heatConductivityRatio);
  }
  else
  {
    heatConductivityRatio = heatConductivitySolid/heatConductivityLiquid;
  }

  if (pp.contains("specificHeatRatio"))
  {
    pp.get("specificHeatRatio", specificHeatRatio);
  }
  else
  {
    specificHeatRatio = specificHeatSolid/specificHeatLiquid;
  }

  if (pp.contains("lewis"))
  {
    pp.get("lewis", lewis);
  }
  else
  {
    lewis = liquidHeatDiffusivity/liquidSoluteDiffusivity;
  }

  if (pp.contains("reynolds"))
  {
    pp.get("reynolds", reynolds);
  }
  else
  {
    reynolds = liquidDensity*liquidHeatDiffusivity/viscosity;
  }

  if (pp.contains("prandtl"))
  {
    pp.get("prandtl", prandtl);
  }
  else
  {
    prandtl = viscosity/liquidDensity*liquidHeatDiffusivity;
  }


  if (pp.contains("darcy"))
  {
    pp.get("darcy", darcy);
  }
  else
  {
    darcy = referencePermeability / (height*height);
  }

  if (pp.contains("nonDimReluctance"))
  {
    pp.get("nonDimReluctance", nonDimReluctance);
  }
  else
  {
    nonDimReluctance = referencePermeability*12/(d*d);
  }


  if (pp.contains("nonDimHeleShawCooling"))
  {
    pp.get("nonDimHeleShawCooling", nonDimHeleShawCooling);
  }
  else
  {
    nonDimHeleShawCooling = heleShawCoolingCoeff * (height*height)/heatConductivityLiquid;
  }

  if (pp.contains("rayleighTemp"))
  {
    pp.get("rayleighTemp", rayleighTemp);
  }
  else
  {
    rayleighTemp = thermalExpansivity*liquidDensity*gravitationalAcceleration*height*
        deltaTemp*referencePermeability / (liquidHeatDiffusivity*viscosity);
  }

  if (pp.contains("rayleighComp"))
  {
    pp.get("rayleighComp", rayleighComposition);
  }
  else
  {
    rayleighComposition = solutalExpansivity*liquidDensity*gravitationalAcceleration*height*
        deltaSalt*referencePermeability /
        (liquidHeatDiffusivity*viscosity);
  }

  timescale = height*height/liquidHeatDiffusivity;

  // Optional parameters
  //  bottomTemp = -999;
  //  topTemp = -999;
  //  pp.query("bottomTemp", bottomTemp);
  //  pp.query("topTemp", topTemp);


  //Some nondimensional parameters and boundary conditions
  if (pp.contains("nonDimVel"))
  {
    pp.get("nonDimVel", nonDimVel);
  }
  else
  {
    nonDimVel = (height/liquidHeatDiffusivity) * V;
  }


  // Default nondimensionalisation is about eutectic
  referenceTemperature = eutecticTemp;
  referenceSalinity = eutecticComposition;

  pp.query("referenceTemperature", referenceTemperature);
  pp.query("referenceSalinity", referenceSalinity);

  thetaEutectic = tempTotheta(eutecticTemp);
  ThetaEutectic = concToTheta(eutecticComposition);

  ThetaInitial = concToTheta(initialComposition);
  ThetaInf = ThetaInitial;


  // Given how we do non dimensionalisation, this is always 1. concToTheta(initialComposition);
  // however we want to be able to set the salinity to be 0 in some case, so allow it like this
  //  if (initialComposition == 0)
  //  {
  //    ThetaInitial = 0.0;
  //  }
  //  else
  //  {
  //    ThetaInitial = -1.0;
  //  }

  //  HTop = -999; HBottom = -99; ThetaTop = -999; ThetaBottom = -999;

  //See if H and C have been specified
  //  pp.query("topEnthalpy", HTop);
  //  pp.query("bottomEnthalpy", HBottom);
  //  pp.query("topBulkConc", ThetaTop);
  //  pp.query("bottomBulkConc", ThetaBottom);

  //  if (HBottom <= -99)
  //  {
  //    Real TBottom = -99;
  //    pp.query("bottomTemperature", TBottom);
  //    if (TBottom > -99)
  //    {
  //      // Assume porosity 1 at bottom
  //      HBottom = stefan*1 + tempTotheta(TBottom);
  //    }
  //  }

  plumeBounds.resize(2);
  if (pp.contains("plumeBounds"))
  {
    std::vector<Real>  temp = std::vector<Real>();
    pp.getarr("plumeBounds", temp, 0, 2);
    //    for (int dir=0; dir<SpaceDim; dir++)
    //    {
    //      plumeBounds[dir] = temp[dir];
    //    }
    plumeBounds[0] = temp[0]; plumeBounds[1] = temp[1];
  }


  pp.query("enthalpyPlume", HPlumeInflow);
  pp.query("bulkConcPlume", ThetaPlumeInflow);


  ParmParse ppBC("bc");

  // Define BC objects
  bcTypeVelLo.resize(SpaceDim, PhysBCUtil::SolidWall);
  bcTypeVelHi.resize(SpaceDim, PhysBCUtil::SolidWall);

  bcTypeBulkConcentrationLo.resize(SpaceDim, PhysBCUtil::Dirichlet);
  bcTypeBulkConcentrationHi.resize(SpaceDim, PhysBCUtil::Dirichlet);

  bcTypeEnthalpyLo.resize(SpaceDim, PhysBCUtil::Dirichlet);
  bcTypeEnthalpyHi.resize(SpaceDim, PhysBCUtil::Dirichlet);

  bcTypeLiquidConcentrationLo.resize(SpaceDim, PhysBCUtil::Dirichlet);
  bcTypeLiquidConcentrationHi.resize(SpaceDim, PhysBCUtil::Dirichlet);

  bcTypeTemperatureLo.resize(SpaceDim, PhysBCUtil::Dirichlet);
  bcTypeTemperatureHi.resize(SpaceDim, PhysBCUtil::Dirichlet);

  bcTypePorosityLo.resize(SpaceDim, PhysBCUtil::Dirichlet);
  bcTypePorosityHi.resize(SpaceDim, PhysBCUtil::Dirichlet);

  bcTypePermeabilityLo.resize(SpaceDim, PhysBCUtil::Dirichlet);
  bcTypePermeabilityHi.resize(SpaceDim, PhysBCUtil::Dirichlet);

  // Specify default values
  // dirichlet BCs in y-direction by default
  bcTypeScalarHi.resize(SpaceDim, PhysBCUtil::Dirichlet);
  bcTypeScalarLo.resize(SpaceDim, PhysBCUtil::Dirichlet);

  // Generally want neuman conditions in x-direction
  bcTypeScalarLo[0] = PhysBCUtil::Neumann;
  bcTypeScalarHi[0] = PhysBCUtil::Neumann;

  bool required = true;

  // Now get custom scalar BCs and apply
  parseBCs("scalarLo", &bcTypeScalarLo, required);
  parseBCs("scalarHi", &bcTypeScalarHi, required);
  parseBCs("velLo", &bcTypeVelLo, required);
  parseBCs("velHi", &bcTypeVelHi, required);

  // By default, set enthalpy and bulk conc to be the same as scalar bcs
  for (int dir=0; dir<SpaceDim; dir++)
  {
    bcTypeBulkConcentrationLo[dir] = bcTypeScalarLo[dir];
    bcTypeBulkConcentrationHi[dir] =  bcTypeScalarHi[dir];
    bcTypeEnthalpyLo[dir] = bcTypeScalarLo[dir];
    bcTypeEnthalpyHi[dir] =  bcTypeScalarHi[dir];

    bcTypeLiquidConcentrationLo[dir] = bcTypeScalarLo[dir];
    bcTypeLiquidConcentrationHi[dir] =  bcTypeScalarHi[dir];
    bcTypeTemperatureLo[dir] = bcTypeScalarLo[dir];
    bcTypeTemperatureHi[dir] =  bcTypeScalarHi[dir];
  }

  // Porosity is slightly different - by defualt, want fixed porosity at top and bottom (already set) and neumann conditions at sides
  bcTypePorosityLo[0] = PhysBCUtil::Neumann;
  bcTypePorosityHi[0] = PhysBCUtil::Neumann;
  bcTypePermeabilityLo[0] = PhysBCUtil::Neumann;
  bcTypePermeabilityHi[0] = PhysBCUtil::Neumann;


  // BC values. User must specify H, C and velocity as a minimum
  // If we have specified any options in our param files, overwrite defaults with them

  //
  // Optional extra BCs, if not specified, these are computed from the H and C already provided

  parseBCs("enthalpyLo", &bcTypeEnthalpyLo);
  parseBCs("enthalpyHi", &bcTypeEnthalpyHi);

  parseBCs("bulkConcentrationLo", &bcTypeBulkConcentrationLo);
  parseBCs("bulkConcentrationHi", &bcTypeBulkConcentrationHi);


    parseBCs("liquidConcentrationLo", &bcTypeLiquidConcentrationLo);
    parseBCs("liquidConcentrationHi", &bcTypeLiquidConcentrationHi);

    parseBCs("temperatureLo", &bcTypeTemperatureLo);
    parseBCs("temperatureHi", &bcTypeTemperatureHi);

    parseBCs("porosityLo", &bcTypePorosityLo);
    parseBCs("porosityHi", &bcTypePorosityHi);

    parseBCs("permeabilityLo", &bcTypePermeabilityLo);
    parseBCs("permeabilityHi", &bcTypePermeabilityHi);


  parseBCVals("velLoVal", bcValVelLo);
  parseBCVals("velHiVal", bcValVelHi);


  parseBCVals("enthalpyLoVal", bcValEnthalpyLo, required);
  parseBCVals("enthalpyHiVal", bcValEnthalpyHi, required);

  parseBCVals("bulkConcentrationLoVal", bcValBulkConcentrationLo, required);
  parseBCVals("bulkConcentrationHiVal", bcValBulkConcentrationHi, required);


  // For the plume
  ::computeBoundingEnergy(HPlumeInflow, ThetaPlumeInflow, HSolidusPlume, HLiquidusPlume, HEutecticPlume,
                          specificHeatRatio, stefan, compositionRatio, waterDistributionCoeff,
                          thetaEutectic, ThetaEutectic);
  ::computeEnthalpyVars(HPlumeInflow, ThetaPlumeInflow, porosityPlume, thetaPlumeInflow,
                        ThetaLPlumeInflow, ThetaSPlumeInflow,
                        HSolidusPlume, HLiquidusPlume, HEutecticPlume,
                        specificHeatRatio, stefan, compositionRatio, waterDistributionCoeff,
                        thetaEutectic, ThetaEutectic);
  permeabilityPlume = this->calculatePermeability(porosityPlume);


  // Now do defaults for the other boundary values
  for (int dir=0; dir<SpaceDim; dir++)
  {

    // Do this for every direction
//    Real solidus,liquidus,eutectic, ThetaS;

    ::computeBoundingEnergy(bcValEnthalpyHi[dir], bcValBulkConcentrationHi[dir], bcValSolidusHi[dir], bcValLiquidusHi[dir], bcValEutecticHi[dir],
                            specificHeatRatio, stefan, compositionRatio, waterDistributionCoeff,
                            thetaEutectic, ThetaEutectic);



    ::computeEnthalpyVars(bcValEnthalpyHi[dir], bcValBulkConcentrationHi[dir], bcValPorosityHi[dir], bcValTemperatureHi[dir], bcValLiquidConcentrationHi[dir], bcValSolidConcentrationHi[dir],
                          bcValSolidusHi[dir], bcValLiquidusHi[dir], bcValEutecticHi[dir],
                          specificHeatRatio, stefan, compositionRatio, waterDistributionCoeff,
                          thetaEutectic, ThetaEutectic);



    ::computeBoundingEnergy(bcValEnthalpyLo[dir], bcValBulkConcentrationLo[dir], bcValSolidusLo[dir], bcValLiquidusLo[dir], bcValEutecticLo[dir],
                            specificHeatRatio, stefan, compositionRatio, waterDistributionCoeff,
                            thetaEutectic, ThetaEutectic);

    ::computeEnthalpyVars(bcValEnthalpyLo[dir], bcValBulkConcentrationLo[dir], bcValPorosityLo[dir], bcValTemperatureLo[dir], bcValLiquidConcentrationLo[dir], bcValSolidConcentrationLo[dir],
                          bcValSolidusLo[dir], bcValLiquidusLo[dir], bcValEutecticLo[dir],
                          specificHeatRatio, stefan, compositionRatio, waterDistributionCoeff,
                          thetaEutectic, ThetaEutectic);


    bcValPermeabilityHi[dir] = calculatePermeability(bcValPorosityHi[dir]);
    bcValPermeabilityLo[dir] =  calculatePermeability(bcValPorosityLo[dir]);


    if (dir == SpaceDim - 1)
    {


      bool useTop = (bcValEnthalpyHi[dir] > bcValEnthalpyLo[dir]);

      if (useTop)
      {
        Hinitial = bcValEnthalpyHi[dir];
        //        ThetaLInitial = bcValLiquidConcentrationHi[dir];
        //        ThetaSInitial = bcValSolidConcentrationHi[dir];;
        //        thetaInitial =  bcValTemperatureHi[dir];;

      }
      else
      {
        Hinitial = bcValEnthalpyLo[dir];;
        //              ThetaLInitial = ThetaLBottom;
        //              ThetaSInitial = ThetaSBottom;
        //              thetaInitial = thetaBottom;
      }


    }
    //    else
    //    {
    //      // Neumann on the sides by default
    //      bcValEnthalpyHi[dir] = 0;
    //      bcValEnthalpyLo[dir]  = 0;
    //      bcValBulkConcentrationHi[dir]  = 0;
    //      bcValBulkConcentrationLo[dir]  = 0;
    //      bcValTemperatureHi[dir]  =0;
    //      bcValTemperatureLo[dir]  =0;
    //      bcValLiquidConcentrationLo[dir]  = 0;
    //      bcValLiquidConcentrationHi[dir]  = 0;
    //      bcValPorosityLo[dir]  = 0;
    //      bcValPorosityHi[dir]  = 0;
    //      bcValPermeabilityLo[dir]  = 0;
    //      bcValPermeabilityHi[dir]  = 0;
    //    }


    //    bcValVelLo[dir] = 0.0;
    //    bcValVelHi[dir] = 0.0;
  }




  // Now BC values


  parseBCVals("liquidConcentrationLoVal", bcValLiquidConcentrationLo);
  parseBCVals("liquidConcentrationHiVal", bcValLiquidConcentrationHi);

  parseBCVals("temperatureLoVal", bcValTemperatureLo);
  parseBCVals("temperatureHiVal", bcValTemperatureHi);

  parseBCVals("porosityLoVal", bcValPorosityLo);
  parseBCVals("porosityHiVal", bcValPorosityHi);

  parseBCVals("permeabilityLoVal", bcValPermeabilityLo);
  parseBCVals("permeabilityHiVal", bcValPermeabilityHi);

  if (physicalProblem == m_cornerFlow)
  {
    // Solid walls at x_lo and y_lo
    // Inflow at y_hi
    // Outflow at x_hi
    bcTypeVelHi[0] = PhysBCUtil::Outflow;
    bcTypeVelHi[1] = PhysBCUtil::Inflow;

    bcTypeVelLo[0]  =PhysBCUtil::noShear;
    bcTypeVelLo[1]  =PhysBCUtil::noShear;


  }
  else if (physicalProblem == m_mushyLayer
      //          && nonDimVel > 0 // I think we want outflow for all mushy layer simulations
  )
  {

    /*
    // Generally solid walls, with a few exceptions:
    if (darcy > 0 && reynolds > 0)
    {
      bcVelLo[1] = PhysBCUtil::Outflow;
    }
    else
    {
      // Darcy flow
      bcVelLo[1] =  PhysBCUtil::OutflowNormal;
    }

    // Scalar BCs - dirichlet at top, inflow/outflow at bottom
    bcBulkConcentrationLo[1] = PhysBCUtil::InflowOutflow;
    bcBulkConcentrationHi[1] = PhysBCUtil::Dirichlet;
    bcEnthalpyLo[1] = PhysBCUtil::Dirichlet;
    bcEnthalpyHi[1] = PhysBCUtil::Dirichlet;

    bcScalarLo[1] = PhysBCUtil::InflowOutflow;
     */


  }






  //Print out parameters
  bool printParams = false;
  if (printParams)
  {
    printParameters();
  }
}

void MushyLayerParams::parseBCs(string a_name, Vector<int>* a_bcHolder, bool required)
{
  std::vector<int>  temp = std::vector<int>();
  ParmParse ppBC("bc");

  if (ppBC.contains(a_name))
  {
    ppBC.getarr(a_name.c_str(), temp, 0, SpaceDim);

    for (int idir=0; idir<SpaceDim; idir++)
    {
      if (temp[idir] !=-1)
      {
        (*a_bcHolder)[idir] = temp[idir];
      }
    }
  }
  else
  {
    if (required)
    {
      pout() << "Can't find BC " << a_name << endl;
      MayDay::Error("Couldn't find BC");
    }
  }
}

void MushyLayerParams::parseBCVals(string a_name, RealVect& a_bcHolder, bool required)
{
  std::vector<Real>  temp = std::vector<Real>();
  ParmParse ppBC("bc");

  if (ppBC.contains(a_name))
  {
    ppBC.getarr(a_name.c_str(), temp, 0, SpaceDim);

    for (int idir=0; idir<SpaceDim; idir++)
    {
      //      if (temp[idir] !=-1)
      //      {
      a_bcHolder[idir] = temp[idir];
      //      }
    }
  }
  else
  {
    if (required)
    {
      pout() << "Can't find BC " << a_name << endl;
      MayDay::Error("Couldn't find BC");
    }
  }
}

void MushyLayerParams::setTime(Real a_time)
{
  m_time = a_time;

  // Update BCs if necessary
  if (m_timeDependentBC == m_sinusoid)

  {

    // dimensionless 0.1 mangitude is about 2 degrees C
    Real T = bcValTemperatureLo[1] - 0.15*sin(2*M_PI*a_time/365);

    // Ensure T isn't less than freezing
    T = max(T, 0.001);

    bcValTemperatureLo[1] = T;
    bcValEnthalpyLo[1] = stefan + T;

    //    pout() << "Set bottom temperature BC = " << T << endl;


    // m_BCtimescale = 1;
    //   m_time = -999;
    //   m_timeDependentBC = m_constant;)
  }


}


Real MushyLayerParams::calculatePermeability(Real liquidFraction)
{
  Real permeability = -1.0;
  Real solidFraction = 1-liquidFraction;
  //    Real referencePerm = params.referencePermeability;

  if (permeabilityFunction == m_pureFluid)
  {
    permeability = 1;
  }
  else if (permeabilityFunction == m_cubicPermeability)
  {
    permeability = pow(liquidFraction,3);
  }
  else if (permeabilityFunction == m_kozenyCarman)
  {
    permeability = pow(liquidFraction,3) / pow(solidFraction,2);
  }
  else if(permeabilityFunction == m_logPermeability)
  {
    permeability = - pow(liquidFraction,2) * log(solidFraction);
  }
  else if(permeabilityFunction == m_porosityPermeability)
  {
    permeability = liquidFraction;
  }
  else if(permeabilityFunction == m_permeabilityXSquared)
  {
    MayDay::Error("Can't calculate permeability for x^2 from just the porosity");
    //                  permeability = x*x;
    //
    //                  Real scale = 0.1;
    //
    //                  // Creates a channel of high permeability in the middle of the domain
    //                  permeability = exp(-pow(x-0.5,2)/scale);
    //
    //                  // Create a hole of low permeabilty in the middle of the domain
    //                  //                                      permeability = 1-exp(-pow(x-0.5,2)/scale)*exp(-pow(z-0.5,2)/scale);
    //
    //                  // Two permeability holes in the left and right of the domain
    //                  scale = 0.03;
    //                  Real zScale = 0.12;
    //                  permeability = exp(-pow(x-0.2,2)/scale)*exp(-pow(z-0.5,2)/zScale);
    //                  permeability = permeability + exp(-pow(x-0.8,2)/scale)*exp(-pow(z-0.5,2)/zScale);
    //                  permeability = 1-permeability;
    //
    //                  //Block flow at the boundaries
    //                  scale = 0.02;
    //                  permeability = exp(-pow(x-0.5,2)/scale);
  }
  else
  {
    permeability = -1;
    MayDay::Error("amrMushyLayer::calculatePermeability() - Unknown permeability function");
  }

  Real finalPermeability = permeability;
  if (heleShaw)
  {
    // Darcy number = K_0 / (h^2) where h is the domain height
    // Cell permeability = d^2/12 where d is the hele-shaw cell width
    // Want to take harmonic mean of cell permeability and this permeability

    //    Real nonDimCellPerm = d*d / (12 * darcy * height*height);
    Real nonDimCellPerm = 1/nonDimReluctance; // = d*d/(12*K_0)
    finalPermeability = 1 / (1/nonDimCellPerm + 1/permeability);
  }

  // Place a cap on the minimum permeability allowed to avoid dividing by 0
  //  Real minPermAllowed = 1e-10;
  //  finalPermeability = Max(finalPermeability, minPermAllowed);

  return finalPermeability;

}

Real MushyLayerParams::
directionalSolidificationMushyZ(Real theta, Real zEutectic)
{
  //Catch a few special situations
  if (stefan == 0)
  {
    //No mushy layer
    return zEutectic;
  }

  Real A = 0.5*(compositionRatio + thetaInf + stefan);
  Real B = sqrt(A*A - compositionRatio * thetaInf - stefan);
  Real alpha = A+B;
  Real beta = A-B;

  Real vel = nonDimVel/this->m_heatDiffusionCoeff;

  Real z =  zEutectic - (1/vel) * (  ((alpha - compositionRatio) / (alpha-beta)) * log((alpha)/(alpha - theta)) +
      ((compositionRatio - beta)  / (alpha-beta)) * log((beta) /(beta - theta))        );

  return z;
}

void MushyLayerParams::
printParameters()
{

  pout() << "Parameters: " << endl;

  SHOW(physicalProblem);
  SHOW(viscosity);
  SHOW(heatConductivityLiquid);
  SHOW(heatConductivitySolid);
  SHOW(specificHeatLiquid);
  SHOW(specificHeatSolid);
  SHOW(liquidDensity);
  SHOW(latentHeatDissolution);
  SHOW(thermalExpansivity);
  SHOW(solutalExpansivity);
  SHOW(eutecticTemp);
//  SHOW(bottomTemp);
//  SHOW(topTemp);
  SHOW(eutecticComposition);
  SHOW(initialComposition);
  SHOW(liquidusSlope);
  SHOW(waterDistributionCoeff);
  SHOW(heleShawCoolingCoeff);
  SHOW(liquidSoluteDiffusivity);
  SHOW(d);
  SHOW(height);
  SHOW(referencePermeability);
  SHOW(gravitationalAcceleration);
  SHOW(V);

  SHOW(deltaTemp);
  SHOW(stefan);
  SHOW(compositionRatio);
  SHOW(liquidHeatDiffusivity);
  SHOW(heatConductivityRatio);
  SHOW(specificHeatRatio);
  SHOW(lewis);
  SHOW(rayleighTemp);
  SHOW(rayleighComposition);
  SHOW(nonDimVel);
  SHOW(nonDimHeleShawCooling);
  SHOW(timescale);
  SHOW(nonDimReluctance);

  SHOW(thetaEutectic);
  SHOW(thetaInf);
  SHOW(thetaInitialLiquidus);
  SHOW(thetaInitial);
  SHOW(thetaInterface);
//  SHOW(thetaBottom);
//  SHOW(thetaTop);
  SHOW(ThetaEutectic);
  SHOW(ThetaInitial);
  SHOW(ThetaInf);
//  SHOW(ThetaTop);
//  SHOW(ThetaBottom);
//  SHOW(ThetaLBottom);
  SHOW(ThetaLInitial);
//  SHOW(ThetaLTop);
//  SHOW(porosityTop);
//  SHOW(porosityBottom);
//  SHOW(HBottom);
//  SHOW(HTop);
  SHOW(Hinitial);
}

Real MushyLayerParams::concToTheta(const Real C)
{
  Real Theta = (C - referenceSalinity) /
      deltaSalt;
  return Theta;
}

Real MushyLayerParams::tempTotheta(const Real T)
{
  Real theta = (T - referenceTemperature) / deltaTemp;
  return theta;
}



int MushyLayerParams::getVelBCType(int dir, Side::LoHiSide side)
{
  // First check to see if we've specified BCs in our param file
  int bcType;

  if (side == Side::Lo)
  {
    bcType = bcTypeVelLo[dir];
  }
  else
  {
    bcType = bcTypeVelHi[dir];
  }

  return bcType;
}


