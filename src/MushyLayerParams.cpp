/*
 * MushyLayerParams.cpp
 *
 *  Created on: 23 May 2016
 *      Author: parkinsonjl
 */

#include "MushyLayerParams.h"
#include "Logging.H"
#include "MushyLayerUtils.H"
#include "PhysBCUtil.H"
#include "ParmParse.H"
#include "phaseDiagram.H"
#include "CH_Timer.H"

// This is for printing out variable names, because c++ doesn't do reflection
#define SHOW(a) LOG(#a << ": " << (a));

/// Darcy timescale
string MushyLayerParams::s_DARCY_TIMESCALE = "Darcy";

/// Diffusive timescale
string MushyLayerParams::s_DIFFUSIVE_TIMESCALE = "Diffusive";

/// Buoyancy timescale
string MushyLayerParams::s_BUOYANCY_TIMESCALE = "Buoyancy";

/// Advective timescale
string MushyLayerParams::s_ADVECTIVE_TIMESCALE = "Advective";

/// Advective velocity scale
string MushyLayerParams::s_ADVECTIVE_VELOCITY_SCALE = "Advective";

/// Darcy velocity scale
string MushyLayerParams::s_DARCY_VELOCITY_SCALE = "Darcy";

void MushyLayerParams::printBCs(string bcName, Vector<int> bcTypeLo,
                                Vector<int> bcTypeHi, RealVect bcValLo,
                                RealVect bcValHi)
{
  char buffer[100];
  sprintf(buffer, "%-10s", bcName.c_str());
  LOG_NOEND(buffer);

  Vector<string> bcTypeNames;
  bcTypeNames.push_back("Fixed");
  bcTypeNames.push_back("No flux");

  for (int dir = 0; dir < SpaceDim; dir++)
  {
    sprintf(buffer, "|%-10s= %02.2f|%-10s= %02.2f",
            m_scalarBCTypes[bcTypeLo[dir]].c_str(), bcValLo[dir],
            m_scalarBCTypes[bcTypeHi[dir]].c_str(), bcValHi[dir]);
    pout() << buffer;
  }
  pout() << endl;
}

// Initialise in initialization list
MushyLayerParams::MushyLayerParams()
    : physicalProblem(PhysicalProblems::m_mushyLayer), viscosity(-999),
      heatConductivityLiquid(-999), heatConductivitySolid(-999),
      specificHeatLiquid(-999), specificHeatSolid(-999), liquidDensity(-999),
      latentHeatDissolution(-999), thermalExpansivity(-999),
      solutalExpansivity(-999), eutecticTemp(-999), eutecticComposition(-999),
      initialComposition(-999), liquidusSlope(-999),
      waterDistributionCoeff(-999), heleShawCoolingCoeff(-999),
      liquidSoluteDiffusivity(-999), d(-999), height(-999), width(-999),
      referencePermeability(-999), gravitationalAcceleration(-999), V(-999),
      deltaTemp(-999), deltaSalt(-999), stefan(-999), compositionRatio(-999),
      liquidHeatDiffusivity(-999), heatConductivityRatio(-999),
      specificHeatRatio(-999), lewis(-999), darcy(-999),
      nonDimReluctance(-999.0), heleShawPermeability(-999.0), reynolds(-999),
      prandtl(0.0), rayleighTemp(-999), rayleighComposition(-999),
      timescale(-999), m_nondimensionalisation(0), m_heatDiffusionCoeff(1),
      m_saltDiffusionCoeff(1), activeTracerDiffusionCoeff(0),
      passiveTracerDiffusionCoeff(0), activeTracerInitVal(0),
      passiveTracerInitVal(0), m_viscosityCoeff(1), m_buoyancyTCoeff(1.0),
      m_buoyancySCoeff(1), m_darcyCoeff(1), m_advectionCoeff(1.0),
      body_force(0.0), nonDimVel(-999), nonDimHeleShawCooling(-999),
      thetaEutectic(-999), thetaInf(-999), thetaInitialLiquidus(-999),
      thetaInitial(-999), thetaInterface(-999), ThetaEutectic(-999),
      ThetaInitial(-999), ThetaInf(-999), ThetaLInitial(-999),
      ThetaSInitial(-999), Hinitial(-999), thetaPlumeInflow(-999),
      HPlumeInflow(-999), ThetaPlumeInflow(-999), ThetaLPlumeInflow(-999),
      ThetaSPlumeInflow(-999), porosityPlume(-999), permeabilityPlume(-999),
      HLiquidusPlume(-999), HEutecticPlume(-999), HSolidusPlume(-999),
      referenceTemperature(-999), referenceSalinity(-999), inflowVelocity(-999),
      pressureHead(0), sinusoidal_temperature_bc_timescale(1),
      sinusoidal_temperature_bc_amplitude(0.5),
      sinusoidal_temperature_bc_av(0.5),
      sinusoidal_temperature_bc_phase_diff(0.0), max_bc_iter(1),
      bc_nonlinear_solve_method(NonlinearBCSolveMethods::picard),
      max_bc_residual(3), bc_relax_coeff(0.5), m_BCAccuracy(1),
      m_pressureBCAccuracy(1), m_time(-999), m_timeDependentBC(m_constant),
      m_BCamplitude(0), m_BCtimescale(1), fixedTempDirection(-999),
      permeabilityFunction(PermeabilityFunctions::m_kozenyCarman),
      heleShaw(false),
      m_porosityFunction(ParamsPorosityFunctions::m_porosityConstant),
      m_viscosityFunction(ViscosityFunction::uniformViscosity),
      max_viscosity(1.0)
{
  m_scalarBCTypes.push_back("fixed");
  m_scalarBCTypes.push_back("noflux");
  m_scalarBCTypes.push_back("open");
  m_scalarBCTypes.push_back("inflow");

  m_vectorBCTypes.push_back("noflow");
  m_vectorBCTypes.push_back("inflowVelocity");
  m_vectorBCTypes.push_back("open");
  m_vectorBCTypes.push_back("outflownormal");
  m_vectorBCTypes.push_back("inflowoutflow");
  m_vectorBCTypes.push_back("noshear");
  m_vectorBCTypes.push_back("symmetry");
  m_vectorBCTypes.push_back("inflowPlume");
  m_vectorBCTypes.push_back("outflowPressureGrad");
  m_vectorBCTypes.push_back("pressureHead");
}

MushyLayerParams::~MushyLayerParams() {}

void MushyLayerParams::computeDerivedBCs()
{
  pout() << "MushyLayerParams::computeDerivedBCs ThetaInitial=" << ThetaInitial
         << ", H top = " << bcValEnthalpyHi[SpaceDim - 1] << endl;
  // Now do defaults for the other boundary values
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    Real bulkCHi = bcValBulkConcentrationHi[dir];
    Real bulkCLo = bcValBulkConcentrationLo[dir];

    if (bcTypeBulkConcentrationHi[dir] == PhysBCUtil::Neumann)
    {
      bulkCHi = ThetaInitial;
    }

    if (bcTypeBulkConcentrationLo[dir] == PhysBCUtil::Neumann)
    {
      bulkCLo = ThetaInitial;
    }

    ::computeBoundingEnergy(
        bcValEnthalpyHi[dir], bulkCHi, bcValSolidusHi[dir],
        bcValLiquidusHi[dir], bcValEutecticHi[dir], specificHeatRatio, stefan,
        compositionRatio, waterDistributionCoeff, thetaEutectic, ThetaEutectic);
    ::computeEnthalpyVars(
        bcValEnthalpyHi[dir], bulkCHi, bcValPorosityHi[dir],
        bcValTemperatureHi[dir], bcValLiquidConcentrationHi[dir],
        bcValSolidConcentrationHi[dir], bcValSolidusHi[dir],
        bcValLiquidusHi[dir], bcValEutecticHi[dir], specificHeatRatio, stefan,
        compositionRatio, waterDistributionCoeff, thetaEutectic, ThetaEutectic);

    ::computeBoundingEnergy(
        bcValEnthalpyLo[dir], bulkCLo, bcValSolidusLo[dir],
        bcValLiquidusLo[dir], bcValEutecticLo[dir], specificHeatRatio, stefan,
        compositionRatio, waterDistributionCoeff, thetaEutectic, ThetaEutectic);
    ::computeEnthalpyVars(
        bcValEnthalpyLo[dir], bulkCLo, bcValPorosityLo[dir],
        bcValTemperatureLo[dir], bcValLiquidConcentrationLo[dir],
        bcValSolidConcentrationLo[dir], bcValSolidusLo[dir],
        bcValLiquidusLo[dir], bcValEutecticLo[dir], specificHeatRatio, stefan,
        compositionRatio, waterDistributionCoeff, thetaEutectic, ThetaEutectic);
    bcValPermeabilityHi[dir] = calculatePermeability(bcValPorosityHi[dir]);
    bcValPermeabilityLo[dir] = calculatePermeability(bcValPorosityLo[dir]);

    if (dir == SpaceDim - 1)
    {
      bool useTop = (bcValEnthalpyHi[dir] > bcValEnthalpyLo[dir]);
      if (useTop)
      {
        Hinitial = bcValEnthalpyHi[dir];
      }
      else
      {
        Hinitial = bcValEnthalpyLo[dir];
      }
    }
  }

  LOG("After computeDerivedBCs: ");
  LOG_NOEND("         ");
  Vector<string> direction;
  direction.push_back("x");
  direction.push_back("y");
  direction.push_back("z");
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    pout() << "| " << direction[dir] << "(lo)         | " << direction[dir]
           << "(hi)         ";
  }
  pout() << endl;

  //    printBCs("H", bcTypeEnthalpyLo, bcTypeEnthalpyHi, bcValEnthalpyLo,
  //    bcValEnthalpyHi); printBCs("C", bcTypeBulkConcentrationLo,
  //    bcTypeBulkConcentrationHi, bcValBulkConcentrationLo,
  //    bcValBulkConcentrationHi);
  printBCs("T", bcTypeTemperatureLo, bcTypeTemperatureHi, bcValTemperatureLo,
           bcValTemperatureHi);
  printBCs("Sl", bcTypeLiquidConcentrationLo, bcTypeLiquidConcentrationHi,
           bcValLiquidConcentrationLo, bcValLiquidConcentrationHi);
  printBCs("porosity", bcTypePorosityLo, bcTypePorosityHi, bcValPorosityLo,
           bcValPorosityHi);
}

void MushyLayerParams::getParameters()
{
  CH_TIME("MushyLayerParams::getParameters()");

  ParmParse ppParams("parameters");
  ParmParse ppMain("main");
  ParmParse ppBio("bio");
  ParmParse ppBC("bc");

  m_nondimensionalisation = 0;
  ppMain.query("nondimensionalisation", m_nondimensionalisation);

  int phys_problem = PhysicalProblems::m_mushyLayer;
  ppParams.query("problem_type", phys_problem);
  physicalProblem = PhysicalProblems(phys_problem);

  int perm_func = permeabilityFunction;
  ppParams.query("permeabilityFunction", perm_func);
  permeabilityFunction = PermeabilityFunctions(perm_func);

  int poros_func = m_porosityFunction;
  ppParams.query("porosityFunction", poros_func);
  m_porosityFunction = ParamsPorosityFunctions(poros_func);

  int viscous_func = m_viscosityFunction;
  ppParams.query("viscosity_function", viscous_func);
  m_viscosityFunction = ViscosityFunction(viscous_func);
  ppParams.query("max_viscosity", max_viscosity);

  ppParams.query("heleShaw", heleShaw);

  // I don't think I necessarily need these, if I'm specifying non dimensional
  // parameters
  ppParams.query("viscosity", viscosity);
  ppParams.query("heatConductivityLiquid", heatConductivityLiquid);
  ppParams.query("heatConductivitySolid", heatConductivitySolid);
  ppParams.query("specificHeatLiquid", specificHeatLiquid);
  ppParams.query("specificHeatSolid", specificHeatSolid);
  ppParams.query("liquidDensity", liquidDensity);
  ppParams.query("latentHeatDissolution", latentHeatDissolution);
  ppParams.query("thermalExpansivity", thermalExpansivity);
  ppParams.query("solutalExpansivity", solutalExpansivity);

  ppParams.query("heleShawCoolingCoeff", heleShawCoolingCoeff);
  ppParams.query("liquidSoluteDiffusivity", liquidSoluteDiffusivity);
  ppParams.query("d", d);
  ppParams.query("height", height);
  ppParams.query("referencePermeability", referencePermeability);
  ppParams.query("gravitationalAcceleration", gravitationalAcceleration);
  ppParams.query("V", V);

  // I do need these, to define the phase diagram
  // Defaults for sea ice:
  initialComposition = 30;
  eutecticComposition = 230;
  eutecticTemp = -23;
  liquidusSlope = -0.1;
  waterDistributionCoeff = 1e-5;

  ppParams.query("eutecticTemp", eutecticTemp);
  ppParams.query("eutecticComposition", eutecticComposition);
  ppParams.query("initialComposition", initialComposition);
  ppParams.query("liquidusSlope", liquidusSlope);
  ppParams.query("waterDistributionCoeff", waterDistributionCoeff);

  ppParams.query("fixedTempDirection", fixedTempDirection);
  ppParams.query("inflowVelocity", inflowVelocity);

  ppBC.query("sinusoidal_temperature_bc_timescale",
             sinusoidal_temperature_bc_timescale);
  ppBC.query("sinusoidal_temperature_bc_amplitude",
             sinusoidal_temperature_bc_amplitude);
  ppBC.query("sinusoidal_temperature_bc_av", sinusoidal_temperature_bc_av);
  ppBC.query("sinusoidal_temperature_bc_phase_diff",
             sinusoidal_temperature_bc_phase_diff);

  ppBC.query("timeDependent", m_timeDependentBC);

  // Derived parameters. Can enforce these if needed (e.g. for benchmarking with
  // a reduced model)
  if (ppParams.contains("deltaSalt"))
  {
    ppParams.get("deltaSalt", deltaSalt);
  }
  else
  {
    deltaSalt = (eutecticComposition - initialComposition);
  }

  if (ppParams.contains("deltaTemp"))
  {
    ppParams.get("deltaTemp", deltaTemp);
  }
  else
  {
    deltaTemp = -liquidusSlope * deltaSalt;
  }

  if (ppParams.contains("stefan"))
  {
    ppParams.get("stefan", stefan);
  }
  else
  {
    stefan = latentHeatDissolution / (specificHeatLiquid * deltaTemp);
  }

  if (ppParams.contains("compositionRatio"))
  {
    ppParams.get("compositionRatio", compositionRatio);
  }
  else
  {
    compositionRatio =
        (1 - waterDistributionCoeff) * eutecticComposition / deltaSalt;
  }

  if (ppParams.contains("liquidHeatDiffusivity"))
  {
    ppParams.get("liquidHeatDiffusivity", liquidHeatDiffusivity);
  }
  else
  {
    liquidHeatDiffusivity =
        heatConductivityLiquid / (liquidDensity * specificHeatLiquid);
  }

  if (ppParams.contains("heatConductivityRatio"))
  {
    ppParams.get("heatConductivityRatio", heatConductivityRatio);
  }
  else
  {
    heatConductivityRatio = heatConductivitySolid / heatConductivityLiquid;
  }

  if (ppParams.contains("specificHeatRatio"))
  {
    ppParams.get("specificHeatRatio", specificHeatRatio);
  }
  else
  {
    specificHeatRatio = specificHeatSolid / specificHeatLiquid;
  }

  if (ppParams.contains("lewis"))
  {
    ppParams.get("lewis", lewis);
  }
  else
  {
    lewis = liquidHeatDiffusivity / liquidSoluteDiffusivity;
  }

  if (ppParams.contains("reynolds"))
  {
    ppParams.get("reynolds", reynolds);
  }
  else
  {
    reynolds = liquidDensity * liquidHeatDiffusivity / viscosity;
  }

  if (ppParams.contains("prandtl"))
  {
    ppParams.get("prandtl", prandtl);
  }
  else
  {
    prandtl = viscosity / liquidDensity * liquidHeatDiffusivity;
  }

  if (ppParams.contains("darcy"))
  {
    ppParams.get("darcy", darcy);
  }
  else
  {
    darcy = referencePermeability / (height * height);
  }

  //  if (pp.contains("nonDimReluctance"))
  //  {
  //    pp.get("nonDimReluctance", nonDimReluctance);
  //  }
  //  else
  //  {
  //    nonDimReluctance = referencePermeability*12/(d*d);
  //  }

  if (ppParams.contains("heleShawPermeability"))
  {
    ppParams.get("heleShawPermeability", heleShawPermeability);
  }
  else if (ppParams.contains("nonDimReluctance"))
  {
    Real rel = 0.0;
    ppParams.get("nonDimReluctance", rel);
    heleShawPermeability = 1.0 / rel;
  }
  else
  {
    heleShawPermeability = d * d / (12 * referencePermeability);
  }

  if (ppParams.contains("nonDimHeleShawCooling"))
  {
    ppParams.get("nonDimHeleShawCooling", nonDimHeleShawCooling);
  }
  else
  {
    nonDimHeleShawCooling =
        heleShawCoolingCoeff * (height * height) / heatConductivityLiquid;
  }

  if (ppParams.contains("rayleighTemp"))
  {
    ppParams.get("rayleighTemp", rayleighTemp);
  }
  else
  {
    rayleighTemp = thermalExpansivity * liquidDensity *
                   gravitationalAcceleration * height * deltaTemp *
                   referencePermeability / (liquidHeatDiffusivity * viscosity);
  }

  if (ppParams.contains("rayleighComp"))
  {
    ppParams.get("rayleighComp", rayleighComposition);
  }
  else
  {
    rayleighComposition = solutalExpansivity * liquidDensity *
                          gravitationalAcceleration * height * deltaSalt *
                          referencePermeability /
                          (liquidHeatDiffusivity * viscosity);
  }

  timescale = height * height / liquidHeatDiffusivity;

  // Optional parameters
  //  bottomTemp = -999;
  //  topTemp = -999;
  //  pp.query("bottomTemp", bottomTemp);
  //  pp.query("topTemp", topTemp);

  // Some nondimensional parameters and boundary conditions
  if (ppParams.contains("nonDimVel"))
  {
    ppParams.get("nonDimVel", nonDimVel);
  }
  else
  {
    //    nonDimVel = (height/liquidHeatDiffusivity) * V;
    nonDimVel = 0.0;
  }

  ppParams.query("BCAccuracy", m_BCAccuracy);
  ppParams.query("pressureBCAccuracy", m_pressureBCAccuracy);

  body_force = 0.0;
  ppParams.query("body_force", body_force);

  // Default nondimensionalisation is about eutectic
  referenceTemperature = eutecticTemp;
  referenceSalinity = eutecticComposition;

  ppParams.query("referenceTemperature", referenceTemperature);
  ppParams.query("referenceSalinity", referenceSalinity);

  thetaEutectic = tempTotheta(eutecticTemp);
  ThetaEutectic = concToTheta(eutecticComposition);

  ThetaInitial = concToTheta(initialComposition);
  ThetaInf = ThetaInitial;

  ppParams.query("initialBulkConc", ThetaInitial);

  plumeBounds.resize(2);
  if (ppParams.contains("plumeBounds"))
  {
    std::vector<Real> temp = std::vector<Real>();
    ppParams.getarr("plumeBounds", temp, 0, 2);
    //    for (int dir=0; dir<SpaceDim; dir++)
    //    {
    //      plumeBounds[dir] = temp[dir];
    //    }
    plumeBounds[0] = temp[0];
    plumeBounds[1] = temp[1];
  }

  ppParams.query("enthalpyPlume", HPlumeInflow);
  ppParams.query("bulkConcPlume", ThetaPlumeInflow);
  
  
  
  
  //// jb - attempting to add multiple dirichlet BC values option for Enthalpy (maybe Salinity in future)
  ////    - currently keeping things in one spot, can redistribute after they're checked
  
  ////        here i'm introducing the boundaries for where the AltVal will occur
  diriBounds.resize(2);
  if (ppParams.contains("diriBounds"))
    {
    std::vector<Real>  temp = std::vector<Real>();
    ppParams.getarr("diriBounds", temp, 0, 2);
    diriBounds[0] = temp[0]; diriBounds[1] = temp[1];
    }

  ////        these are the AltVals
  parseBCVals("enthalpyAltLoVal", bcAltValEnthalpyLo, required);
  parseBCVals("enthalpyAltHiVal", bcAltValEnthalpyHi, required);

  ////        recompute bounding energy and enthalpy vars???
  
  
  

  //  ParmParse ppBC("bc");

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
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    bcTypeBulkConcentrationLo[dir] = bcTypeScalarLo[dir];
    bcTypeBulkConcentrationHi[dir] = bcTypeScalarHi[dir];
    bcTypeEnthalpyLo[dir] = bcTypeScalarLo[dir];
    bcTypeEnthalpyHi[dir] = bcTypeScalarHi[dir];

    bcTypeLiquidConcentrationLo[dir] = bcTypeScalarLo[dir];
    bcTypeLiquidConcentrationHi[dir] = bcTypeScalarHi[dir];
    bcTypeTemperatureLo[dir] = bcTypeScalarLo[dir];
    bcTypeTemperatureHi[dir] = bcTypeScalarHi[dir];
  }

  // Porosity is slightly different - by defualt, want fixed porosity at top and
  // bottom (already set) and neumann conditions at sides
  bcTypePorosityLo[0] = PhysBCUtil::Neumann;
  bcTypePorosityHi[0] = PhysBCUtil::Neumann;
  bcTypePermeabilityLo[0] = PhysBCUtil::Neumann;
  bcTypePermeabilityHi[0] = PhysBCUtil::Neumann;

  // BC values. User must specify H, C and velocity as a minimum
  // If we have specified any options in our param files, overwrite defaults
  // with them

  //
  // Optional extra BCs, if not specified, these are computed from the H and C
  // already provided

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

  for (int dir = 0; dir < SpaceDim; dir++)
  {
    bcValPressureHi[dir] = 0.0;
    bcValPressureLo[dir] = 0.0;
  }

  // Legacy option
  if (ppParams.contains("pressureHead"))
  {
    ppParams.query("pressureHead", pressureHead);
    bcValPressureHi[SpaceDim - 1] = pressureHead;
    bcValPressureLo[SpaceDim - 1] = 0.0;
  }

  parseBCVals("bcValPressureHi", bcValPressureHi);
  parseBCVals("bcValPressureLo", bcValPressureLo);

  // Now get extra parameters associated with boundary conditions
  m_bc_noFluxLimit = BCInfo(string("NoFluxLimit"), false);
  m_bc_bTref = BCInfo(string("TRef"), false);
  m_bc_b = BCInfo(string("b"), false);
  m_bc_a = BCInfo(string("a"), false);

  max_bc_iter = 2;
  ppMain.query("max_bc_iter", max_bc_iter);

  max_bc_residual = 1e-5;
  ppMain.query("max_bc_residual", max_bc_residual);

  bc_relax_coeff = 0.5;
  ppMain.query("bc_relax_coeff", bc_relax_coeff);

  bc_nonlinear_solve_method = NonlinearBCSolveMethods::picard;
  ppMain.query("nonlinear_bc_solve_method", bc_nonlinear_solve_method);

  // For the plume
  ::computeBoundingEnergy(HPlumeInflow, ThetaPlumeInflow, HSolidusPlume,
                          HLiquidusPlume, HEutecticPlume, specificHeatRatio,
                          stefan, compositionRatio, waterDistributionCoeff,
                          thetaEutectic, ThetaEutectic);
  ::computeEnthalpyVars(HPlumeInflow, ThetaPlumeInflow, porosityPlume,
                        thetaPlumeInflow, ThetaLPlumeInflow, ThetaSPlumeInflow,
                        HSolidusPlume, HLiquidusPlume, HEutecticPlume,
                        specificHeatRatio, stefan, compositionRatio,
                        waterDistributionCoeff, thetaEutectic, ThetaEutectic);
  permeabilityPlume = this->calculatePermeability(porosityPlume);

  // Now do defaults for the other boundary values
  computeDerivedBCs();

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

    bcTypeVelLo[0] = PhysBCUtil::noShear;
    bcTypeVelLo[1] = PhysBCUtil::noShear;
  }
  else if (physicalProblem == m_mushyLayer)
  {
  }

  // Now sort out nondimensionalisation
  if (!isDarcyBrinkman() &&
      m_nondimensionalisation != m_diffusiveTime_advectiveVel)
  {
    MayDay::Error("Need to choose diffusive timescale/advection velocity scale "
                  "for solving Darcy flow.");
  }

  if (m_nondimensionalisation == m_darcyTime_advectiveVel)
  {
    LOG("Darcy timescale, advective velocity scale");

    // To avoid dividing by 0 when Da = Pr = 0
    if (darcy == prandtl)
    {
      m_heatDiffusionCoeff = 1;
    }
    else
    {
      m_heatDiffusionCoeff = darcy / prandtl;
    }
    m_saltDiffusionCoeff = m_heatDiffusionCoeff / lewis;
    m_viscosityCoeff = darcy;

    m_buoyancyTCoeff = rayleighTemp * darcy * darcy * prandtl;
    m_buoyancySCoeff = rayleighComposition * darcy * darcy * prandtl;
    m_darcyCoeff = 1.0;
    m_advectionCoeff = 1.0;
  }
  else if (m_nondimensionalisation == m_diffusiveTime_advectiveVel)
  {
    LOG("Diffusive timescale, advective velocity scale");

    m_heatDiffusionCoeff = 1.0;
    m_saltDiffusionCoeff = 1 / lewis;
    m_viscosityCoeff = prandtl;
    m_buoyancyTCoeff = prandtl * rayleighTemp;
    m_buoyancySCoeff = prandtl * rayleighComposition;
    m_darcyCoeff = prandtl / darcy;
    m_advectionCoeff = 1.0;

    if (!isDarcyBrinkman())
    {
      m_buoyancyTCoeff = rayleighTemp;
      m_buoyancySCoeff = rayleighComposition;
    }
  }
  else if (m_nondimensionalisation == m_darcyTime_darcyVel)
  {

    LOG("Darcy timescale, darcy velocity scale");

    m_heatDiffusionCoeff = darcy / prandtl;
    m_saltDiffusionCoeff = m_heatDiffusionCoeff / lewis;
    m_viscosityCoeff = darcy;
    m_buoyancyTCoeff = 1.0; // rayleighTemp*darcy*darcy*prandtl;
    m_buoyancySCoeff = rayleighComposition /
                       rayleighTemp; // rayleighComposition*darcy*darcy*prandtl;
    m_darcyCoeff = 1.0;
    m_advectionCoeff = rayleighTemp * darcy * darcy / prandtl;
  }
  else if (m_nondimensionalisation == m_advectiveTime_darcyVel)
  {
    LOG("Advective timescale, darcy velocity scale");

    m_heatDiffusionCoeff = 1 / (darcy * rayleighTemp);
    m_saltDiffusionCoeff = m_heatDiffusionCoeff / lewis;
    m_viscosityCoeff = prandtl / (darcy * rayleighTemp);
    m_buoyancyTCoeff = prandtl / (darcy * rayleighTemp);
    m_buoyancySCoeff = m_buoyancyTCoeff * (rayleighComposition / rayleighTemp);
    m_darcyCoeff = prandtl / (darcy * darcy * rayleighTemp);
    m_advectionCoeff = 1.0;
  }
  else if (m_nondimensionalisation == m_buoyancyTime_advectiveVel)
  {
    LOG("Buoyancy timescale, advective velocity scale");

    Real R = rayleighComposition;
    if (R == 0)
    {
      if (rayleighTemp == 0)
      {
        MayDay::Error("MushyLayerParams - Can't nondimensionalise with "
                      "buoyancy timescale if RaT=RaC=0!!");
      }
      R = rayleighTemp;
    }

    m_heatDiffusionCoeff = 1 / sqrt(R * prandtl);
    m_saltDiffusionCoeff = m_heatDiffusionCoeff / lewis;

    m_viscosityCoeff = sqrt(prandtl / R);

    m_darcyCoeff = (1 / darcy) * sqrt(prandtl / R);
    m_advectionCoeff = 1.0;

    // Avoid dividing by zero
    if (rayleighComposition != 0)
    {
      m_buoyancySCoeff = 1;
      m_buoyancyTCoeff = rayleighTemp / rayleighComposition;
    }
    else
    {
      m_buoyancySCoeff = 0;
      m_buoyancyTCoeff = 1;
    }
  }
  else
  {
    MayDay::Error("Unknown non dimensionalisation");
  }

  // Finally, option to manually set certain terms if we want (for testing
  // purposes)
  ppMain.query("heatDiffusionCoeff", m_heatDiffusionCoeff);
  ppMain.query("saltDiffusionCoeff", m_saltDiffusionCoeff);
  ppMain.query("viscosityCoeff", m_viscosityCoeff);
  ppMain.query("buoyancyTCoeff", m_buoyancyTCoeff);
  ppMain.query("buoyancySCoeff", m_buoyancySCoeff);
  ppMain.query("darcyCoeff", m_darcyCoeff);
  ppMain.query("advectionCoeff", m_advectionCoeff);

  ppBio.query("activeTracerDiffusionCoeff", activeTracerDiffusionCoeff);
  ppBio.query("passiveTracerDiffusionCoeff", passiveTracerDiffusionCoeff);

  ppBio.query("activeTracerInitVal", activeTracerInitVal);
  ppBio.query("passiveTracerInitVal", passiveTracerInitVal);

  // uncomment to print the parameters
  //  printParameters();

  LOG("BCs used: ");
  LOG_NOEND("         ");
  Vector<string> direction;
  direction.push_back("x");
  direction.push_back("y");
  direction.push_back("z");
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    pout() << "| " << direction[dir] << "(lo)         | " << direction[dir]
           << "(hi)         ";
  }
  pout() << endl;

  printBCs("H", bcTypeEnthalpyLo, bcTypeEnthalpyHi, bcValEnthalpyLo,
           bcValEnthalpyHi);
  printBCs("C", bcTypeBulkConcentrationLo, bcTypeBulkConcentrationHi,
           bcValBulkConcentrationLo, bcValBulkConcentrationHi);
  printBCs("T", bcTypeTemperatureLo, bcTypeTemperatureHi, bcValTemperatureLo,
           bcValTemperatureHi);
  printBCs("Sl", bcTypeLiquidConcentrationLo, bcTypeLiquidConcentrationHi,
           bcValLiquidConcentrationLo, bcValLiquidConcentrationHi);
  printBCs("porosity", bcTypePorosityLo, bcTypePorosityHi, bcValPorosityLo,
           bcValPorosityHi);
}

void MushyLayerParams::parseBCs(string a_name, Vector<int> *a_bcHolder,
                                bool required)
{
  std::vector<int> temp = std::vector<int>();
  std::vector<string> temp_str = std::vector<string>();
  ParmParse ppBC("bc");

  if (ppBC.contains(a_name))
  {
    string val;
    ppBC.get(a_name.c_str(), val);

    // determine how we should read the bcs in
    if (is_integer(val))
    {
      ppBC.getarr(a_name.c_str(), temp, 0, SpaceDim);
    }
    else
    {
      ppBC.getarr(a_name.c_str(), temp_str, 0, SpaceDim);

      temp.resize(temp_str.size());

      // Convert strings to numbers
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (temp_str[idir] == "noflux")
        {
          temp[idir] = 1;
        }
        else if (temp_str[idir] == "fixed")
        {
          temp[idir] = 0;
        }
        else if (temp_str[idir] == "open")
        {
          temp[idir] = 2;
        }
        else if (temp_str[idir] == "inflow")
        {
          temp[idir] = 3;
        }
        else if (temp_str[idir] == "mixed")
        {
          temp[idir] = 9;
        }
        /**
         * SolidWall,
                Inflow,
                Outflow,
                OutflowNormal, // only a normal velocity
                VelInflowOutflow, // both inflow and outflow possible
                noShear,
                Symmetry,
                VelInflowPlume,
                OutflowPressureGrad,
                PressureHead,
         */
        else if (temp_str[idir] == "solidwall" || temp_str[idir] == "noflow")
        {
          temp[idir] = 0;
        }
        else if (temp_str[idir] == "inflowVelocity")
        {
          temp[idir] = 1;
        }
        else if (temp_str[idir] == "outflow" || temp_str[idir] == "open")
        {
          temp[idir] = 2;
        }
        else if (temp_str[idir] == "outflownormal")
        {
          temp[idir] = 3;
        }
        else if (temp_str[idir] == "inflowoutflow")
        {
          temp[idir] = 4;
        }
        else if (temp_str[idir] == "noshear")
        {
          temp[idir] = 5;
        }
        else if (temp_str[idir] == "symmetry")
        {
          temp[idir] = 6;
        }
        else if (temp_str[idir] == "inflowPlume")
        {
          temp[idir] = 7;
        }
        else if (temp_str[idir] == "outflowPressureGrad")
        {
          temp[idir] = 8;
        }
        else if (temp_str[idir] == "pressureHead")
        {
          temp[idir] = 9;
        }

        else
        {
          LOG("Unknown BC " << temp_str[idir]);
          temp[idir] = 0;
        }
      }
    }

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      // We're parsing bc types here, so the values should be >= 0
      if (temp[idir] != -1)
      {
        (*a_bcHolder)[idir] = temp[idir];
      }
    }
  }
  else
  {
    if (required)
    {
      LOG("Can't find BC " << a_name);
      MayDay::Error("Couldn't find BC");
    }
  }
}

void MushyLayerParams::parseBCVals(string a_name, RealVect &a_bcHolder,
                                   bool required)
{
  std::vector<Real> temp = std::vector<Real>();
  ParmParse ppBC("bc");

  if (ppBC.contains(a_name))
  {
    ppBC.getarr(a_name.c_str(), temp, 0, SpaceDim);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_bcHolder[idir] = temp[idir];
    }
  }
  else
  {
    if (required)
    {
      LOG("Can't find BC " << a_name);
      MayDay::Error("Couldn't find BC");
    }
  }
}

bool MushyLayerParams::isViscous() { return (m_viscosityCoeff > 0); }

bool MushyLayerParams::isDarcyBrinkman()
{
  if (darcy > 1e-10)
  {
    return true;
  }
  else
  {
    return false;
  }
}

void MushyLayerParams::setTime(Real a_time) { m_time = a_time; }

Real MushyLayerParams::calculatePermeability(Real liquidFraction)
{
  Real permeability = -1.0;
  Real solidFraction = 1 - liquidFraction;
  //    Real referencePerm = params.referencePermeability;

  if (permeabilityFunction == PermeabilityFunctions::m_pureFluid)
  {
    permeability = 1;
  }
  else if (permeabilityFunction ==
           PermeabilityFunctions::m_cubicPermeability)
  {
    permeability = pow(liquidFraction, 3);
  }
  else if (permeabilityFunction == PermeabilityFunctions::m_kozenyCarman)
  {
    permeability = pow(liquidFraction, 3) / pow(solidFraction, 2);
  }
  else if (permeabilityFunction == PermeabilityFunctions::m_logPermeability)
  {
    permeability = -pow(liquidFraction, 2) * log(solidFraction);
  }
  else if (permeabilityFunction ==
           PermeabilityFunctions::m_porosityPermeability)
  {
    permeability = liquidFraction;
  }
  else if (permeabilityFunction ==
           PermeabilityFunctions::m_permeabilityXSquared)
  {
    MayDay::Error(
        "Can't calculate permeability for x^2 from just the porosity");
  }
  else
  {
    permeability = -1;
    MayDay::Error("amrMushyLayer::calculatePermeability() - Unknown "
                  "permeability function");
  }

  Real finalPermeability = permeability;
  if (heleShaw)
  {
    // Darcy number = K_0 / (h^2) where h is the domain height
    // Cell permeability = d^2/12 where d is the hele-shaw cell width
    // Want to take harmonic mean of cell permeability and this permeability

    //    Real nonDimCellPerm = d*d / (12 * darcy * height*height);
    //    Real nonDimCellPerm = 1/nonDimReluctance; // = d*d/(12*K_0)

    finalPermeability = 1 / (1 / heleShawPermeability + 1 / permeability);
  }

  // Place a cap on the minimum permeability allowed to avoid dividing by 0
  //  Real minPermAllowed = 1e-10;
  //  finalPermeability = Max(finalPermeability, minPermAllowed);

  return finalPermeability;
}

Real MushyLayerParams::directionalSolidificationMushyZ(Real theta,
                                                       Real zEutectic)
{
  // Catch a few special situations
  if (stefan == 0)
  {
    // No mushy layer
    return zEutectic;
  }

  Real A = 0.5 * (compositionRatio + thetaInf + stefan);
  Real B = sqrt(A * A - compositionRatio * thetaInf - stefan);
  Real alpha = A + B;
  Real beta = A - B;

  Real vel = nonDimVel / this->m_heatDiffusionCoeff;

  Real z =
      zEutectic - (1 / vel) * (((alpha - compositionRatio) / (alpha - beta)) *
                                   log((alpha) / (alpha - theta)) +
                               ((compositionRatio - beta) / (alpha - beta)) *
                                   log((beta) / (beta - theta)));

  return z;
}

void MushyLayerParams::printParameters()
{

  LOG("Parameters: ");

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
  //  SHOW(nonDimReluctance);
  SHOW(heleShawPermeability);

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
  Real Theta = (C - referenceSalinity) / deltaSalt;
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

void MushyLayerParams::writeToHDF5(HDF5HeaderData &a_header) const
{
  a_header.m_real["C"] = compositionRatio;
  a_header.m_real["k"] = heatConductivityRatio;
  a_header.m_real["cp"] = specificHeatRatio;
  a_header.m_real["pc"] = waterDistributionCoeff;
  a_header.m_real["S"] = stefan;
  a_header.m_real["L"] = lewis;
  a_header.m_real["D"] = darcy;
  a_header.m_string["timescale"] = getTimescale();
  a_header.m_string["velocity_scale"] = getVelocityScale();
}

string MushyLayerParams::getTimescale() const
{
  switch (m_nondimensionalisation)
  {
  case m_advectiveTime_darcyVel:
    return s_ADVECTIVE_TIMESCALE;
    break;

  case m_diffusiveTime_advectiveVel:
    return s_DIFFUSIVE_TIMESCALE;
    break;

  case m_darcyTime_advectiveVel:
    return s_DARCY_TIMESCALE;
    break;

  case m_darcyTime_darcyVel:
    return s_DARCY_TIMESCALE;
    break;

  case m_buoyancyTime_advectiveVel:
    return s_BUOYANCY_TIMESCALE;
    break;

  default:
    MayDay::Error("Unknown nondimensionalisation");
    return "";
    break;
  }
}

// Utility function
Real MushyLayerParams::computePorosity(Real H, Real C)
{
  CH_TIME("MushyLayerParams::computePorosity");

  Real porosity = ::computePorosity(
      H, C, compositionRatio, specificHeatRatio, stefan, waterDistributionCoeff,
      specificHeatRatio, thetaEutectic, ThetaEutectic);

  return porosity;
}

Real MushyLayerParams::compute_dHdT(Real H, Real C)
{
  Real H_e, H_s, H_l, dHdT;
  H_e = std::nan("1");
  H_s = std::nan("1");
  H_l = std::nan("1");

  ::computeBoundingEnergy(H_e, C, H_s, H_l, H_e, specificHeatRatio, stefan,
                          compositionRatio, waterDistributionCoeff,
                          thetaEutectic, ThetaEutectic);

  if (H <= H_s)
  {
    dHdT = 1 / specificHeatRatio;
  }
  else if (H > H_s && H <= H_e)
  {
    dHdT = 0;
  }
  else if (H > H_e && H < H_l)
  {
    Real porosity =
        computePorosityMushyLayer(H, C, compositionRatio, specificHeatRatio,
                                  stefan, waterDistributionCoeff);
    //      theta = - (C + compositionRatio*(1-porosity)) / (porosity +
    //      waterDistributionCoeff*(1-porosity));
    Real A = compositionRatio * (specificHeatRatio - 1) - stefan;
    Real B = H + compositionRatio * (1 - 2 * specificHeatRatio) -
             C * (specificHeatRatio - 1);
    Real Cc = specificHeatRatio * (compositionRatio + C);

    dHdT = pow(porosity, 2) * (C + compositionRatio) * (-2 * A) /
           (1 + B / sqrt(pow(B, 2) - 4 * A * Cc));
  }
  else
  {
    dHdT = 1;
  }

  return dHdT;
}

Real MushyLayerParams::computeTemperature(Real H, Real C)
{
  CH_TIME("MushyLayerParams::computeTemperature");

  Real temperature = ::computeTemperature(
      H, C, compositionRatio, specificHeatRatio, stefan, waterDistributionCoeff,
      specificHeatRatio, thetaEutectic, ThetaEutectic);

  return temperature;
}

void MushyLayerParams::computeDiagnosticVars(Real H, Real C, Real T,
                                             Real porosity, Real Cl, Real Cs)
{
  Real H_s, H_l, H_e;
  ::computeBoundingEnergy(H, C, H_s, H_l, H_e, specificHeatRatio, stefan,
                          compositionRatio, waterDistributionCoeff,
                          thetaEutectic, ThetaEutectic);
  ::computeEnthalpyVars(H, C, porosity, T, Cl, Cs, H_s, H_l, H_e,
                        specificHeatRatio, stefan, compositionRatio,
                        waterDistributionCoeff, thetaEutectic, ThetaEutectic);
}

string MushyLayerParams::getVelocityScale() const
{
  switch (m_nondimensionalisation)
  {
  case m_advectiveTime_darcyVel:
    return s_DARCY_VELOCITY_SCALE;
    break;

  case m_diffusiveTime_advectiveVel:
    return s_ADVECTIVE_VELOCITY_SCALE;
    break;

  case m_darcyTime_advectiveVel:
    return s_ADVECTIVE_VELOCITY_SCALE;
    break;

  case m_darcyTime_darcyVel:
    return s_DARCY_VELOCITY_SCALE;
    break;

  case m_buoyancyTime_advectiveVel:
    return s_ADVECTIVE_VELOCITY_SCALE;
    break;

  default:
    MayDay::Error("Unknown nondimensionalisation");
    return "";
    break;
  }
}
