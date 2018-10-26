/*
 * MushyLayerParams.h
 *
 *  Created on: 23 May 2016
 *      Author: parkinsonjl
 */

#ifndef SRC_MUSHYLAYERPARAMS_H_
#define SRC_MUSHYLAYERPARAMS_H_

#include <cmath>
#include "ParmParse.H"
#include "LoHiSide.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "RealVect.H"



/// Class to handle the physical parameters of a mushy layer simulation
/**
 * This includes material properties, boundary conditions, initial conditions and more.
 * Properties are read in and then non-dimensionalised here.
 */
class MushyLayerParams {
public:
	MushyLayerParams();
	virtual ~MushyLayerParams();

	// So many parameters

	int
	/// Benchmark problem to solve
        physicalProblem;

	Real
	/// Viscosity, \f$ \eta \f$, used in rayleigh number calculations
	viscosity,

	/// Liquid heat conductivity, \f$ k_l \f$.
	heatConductivityLiquid,

	/// Solid heat conductivity, \f$ k_s \f$.
	heatConductivitySolid,

	/// Specific heat capacity of the liquid phase, \f$ c_l \f$.
	specificHeatLiquid,

	/// Specific heat capacity of the solid phase, \f$c_s\f$.
	specificHeatSolid,

	/// Density of the liquid
	liquidDensity,

	/// Latent heat for liquid->solid phase change
	latentHeatDissolution,

	/// Thermal expansivity (\f$\mbox{K}^{-1} \f$)
	thermalExpansivity,

	/// Solutal expansitivity ( \f$\mbox{(g/kg)}^{-1}\f$ )
	solutalExpansivity,

	/// Temperature at the eutectic point
	eutecticTemp,

	/// Temperature at the bottom of the domain
//	bottomTemp,

	/// Temperature at the top of the domain
//	topTemp,

	/// Bulk concentration at the eutectic point
	eutecticComposition,

	/// Initial bulk concentration
	initialComposition,

	/// Slope of the linearised liquidus \f$\mbox{ K (g/kg)}^{-1} \f$
	liquidusSlope,

	/// Water distribution coefficient, \f$ p_c \f$
	/**
	 * ratio of solid to liquid concentration, \f$ S_s/S_l \f$.
	 */
	waterDistributionCoeff,

	/// Hele-Shaw cooling coefficient - not currently in use
	heleShawCoolingCoeff,

	/// Diffusivity of salt in the liquid phase, \f$ D_l \f$
	liquidSoluteDiffusivity,

	/// Hele-Shaw cell width - not currently in use
	d,

	/// Domain height
	height,

	/// Domain width - not currently in use
	width,

	/// Reference permeability \f$ K_0 \f$
	referencePermeability,

	/// Gravitational acceleration, \f$ g \f$
	gravitationalAcceleration,

	/// Frame advection velocity \f$ V \f$
	V;

	//Derived parameters

	/// Temperature differenece used for nondimensionalisation
	/**
	 * \f$ \Delta T = - \Gamma \Delta S \f$
	 */
	Real deltaTemp,

	/// Salt difference used for nondimensionalisation
	/**
	 * \f$ \Delta S = S_e - S_i \f$
	 */
	deltaSalt,

	/// Stefan number
	/**
	 * \f[  = \frac{L}{c_{p,l} (T_i - T_e)} \f]
	 */
	stefan,

	/// Composition ratio = \f$ \frac{S_e (p_c - 1)}{S_i -S_e} \f$
	compositionRatio,

	/// Heat diffusivity in the liquid phase
	liquidHeatDiffusivity,

	/// Thermal conductivity ratio \f$ \bar{k} = k_s/k_l \f$
	heatConductivityRatio,

	/// Specific heat ratio \f$ \bar{c} = c_s/c_l \f$
	specificHeatRatio,

	/// Lewis number \f$ Le = \kappa_l/D_l \f$
	lewis,

	/// Darcy number \f$ K_0/h^2 \f$
	/**
	 * Set to zero to turn off viscosity.
	 */
	darcy,

	/// Non dimensional reluctance of the hele-shaw cell
	/**
	 * Inverse of permeability. Equivalent to the Darcy number if Katz & Worster (2008)
	 * Reluctance = Da * 12 * (h/d)^2
	 */
	nonDimReluctance,



	/// Reynolds number \f$ \rho_0 \kappa_l / \eta \f$. Not used any more.
	reynolds,

	/// Prandtl number \f$ \eta / \rho_0 \kappa_l \f$
	prandtl,

	/// Rayleigh number for temperature contributions to buoyancy
	/**
	 * \f[ Ra_T = \frac{\rho_0 g h K_0 \alpha (T_i - T_e)}{\kappa_l \eta} \f]
	 */
	rayleighTemp,

	/// Rayleigh number for salinity contributions to buoyancy
	/**
	         * \f[ Ra_T = \frac{\rho_0 g h K_0 \beta (S_i - S_e)}{\kappa_l \eta} \f]
	         */
	rayleighComposition,

	/// Timescale for nondimensionalisation - not currently used
	/**
	 * \f$ \tau = h^2 / \kappa_l \f$
	 */
	timescale;

	// For defining the equations:
//	  bool m_darcyTimescale;
	  int m_nondimensionalisation;

	  enum nondimensionalisations
	  {
	    m_diffusiveTime_advectiveVel,
	    m_darcyTime_advectiveVel,
	    m_darcyTime_darcyVel,
	    m_advectiveTime_darcyVel,
	    m_buoyancyTime_advectiveVel,

	    m_num_nondimensionalisations
	  };

	  Real m_heatDiffusionCoeff,
	  m_saltDiffusionCoeff,
	  m_viscosityCoeff,
	  m_buoyancyTCoeff,
	  m_buoyancySCoeff,
	  m_darcyCoeff,
	  m_advectionCoeff;

	//Nondimensional parameters and boundary conditions

	/// Dimensionless frame advection velocity
	Real nonDimVel,

	/// Dimensionless hele-shaw cooling coefficient
	nonDimHeleShawCooling,

	/// Dimensionless temperature at the eutectic, \f$ \theta_e \f$
	thetaEutectic,

	/// Dimensionless temperature at infinity, \f$ \theta_\infty \f$
	thetaInf,

	/// Dimensionless initial temperature calculated from liquidus (not currently used)
	thetaInitialLiquidus,

	/// Dimensionless initial temperature \f$ \theta_i \f$
	thetaInitial,

	/// Dimensionless temperature at mush-liquid interface
	/**
	 * Used particularly for solidification without flow benchmark problem
	 */
	thetaInterface,

	/// Dimensionless temperature at the bottom of the domain \f$ \theta_b \f$
//	thetaBottom,

	/// Dimensionless temperature at the top of the domain \f$ \theta_t \f$
//	thetaTop,

	/// Dimensionless bulk concentration at the eutectic \f$ \Theta_e \f$
	ThetaEutectic,

	/// Dimensionless initial bulk concentration  \f$ \Theta_i \f$
	ThetaInitial,

	/// Dimensionless far field bulk concentration \f$ \Theta_\infty \f$
	ThetaInf,

	/// Dimensionless bulk concentration at the top of the domain \f$ \Theta_t \f$
//	ThetaTop,

	/// Dimensionless bulk concentration at the bottom of the domain \f$ \Theta_b \f$
//	ThetaBottom,

	/// Dimensionless liquid concentration at the bottom of the domain \f$ \Theta_l \f$
//	ThetaLBottom,

	/// Dimensionless initial liquid concentration \f$ \Theta_{l,i} \f$
	ThetaLInitial,

	/// Dimensionless liquid concentration at the top of the domain \f$ \Theta_{l,t} \f$
//	ThetaLTop,

	/// Dimensionless solid concentration at the top of the domain \f$ \Theta_{s,t} \f$
//	ThetaSTop,

	/// Dimensionless solid concentration at the bottom of the domain \f$ \Theta_{s,b} \f$
//	ThetaSBottom,

	/// Dimensionless initial solid concentration \f$ \Theta_{s,i} \f$
	ThetaSInitial,

	/// Porosity at the top of the domain
//	porosityTop,

	/// Porosity at the bottom of the domain
//	porosityBottom,

	/// Permeability at the top of the domain
//	permeabilityTop,

	/// Permeability at the bottom of the domain
//	permeabilityBottom,

	/// Dimensionless enthalpy at the bottom of the domain
//	HBottom,

	/// Dimensionless enthalpy at the top of the domain
//	HTop,

	/// Dimensionless initial enthalpy
	Hinitial,

	/// Dimensionless enthalpy solidus at the bottom of the domain
//	HSolidusBottom,

	/// Dimensionless enthalpy solidus at the top of the domain
//	HSolidusTop,

	/// Dimensionless eutectic enthalpy at the bottom of the domain
//	HEutecticBottom,

	/// Dimensionless eutectic enthalpy at the top of the domain
//	HEutecticTop,

	/// Dimensionless enthalpy liquidus at the bottom of the domain
//	HLiquidusBottom,

	/// Dimensionless enthalpy liquidus at the top of the domain
//	HLiquidusTop,

        /// For plumes
        thetaPlumeInflow,
        HPlumeInflow,
        ThetaPlumeInflow,
        ThetaLPlumeInflow,
        ThetaSPlumeInflow,
        porosityPlume,
        permeabilityPlume,
        HLiquidusPlume,
        HEutecticPlume,
        HSolidusPlume;

	// To determine how we do nondimensionalisation
	Real referenceTemperature, referenceSalinity;

	/// Inflow velocity, when required
	Real inflowVelocity;

	/// Specify start and end of an inflow plume
	Vector<Real> plumeBounds;

	/// Velocity boundary conditions (lo side, for each spatial direction)
	Vector<int> bcTypeVelLo,

	/// Velocity boundary conditions (hi side, for each spatial direction)
	bcTypeVelHi,

	/// Enthalpy boundary conditions (lo side, for each spatial direction)
	bcTypeEnthalpyLo,

	/// Enthalpy boundary conditions (hi side, for each spatial direction)
	bcTypeEnthalpyHi,

	/// Bulk concentration boundary conditions (lo side, for each spatial direction)
	bcTypeBulkConcentrationLo,

	/// Bulk concentration boundary conditions (hi side, for each spatial direction)
	bcTypeBulkConcentrationHi,

	bcTypeTemperatureLo, bcTypeTemperatureHi,
	bcTypeLiquidConcentrationLo, bcTypeLiquidConcentrationHi,
	bcTypePorosityLo, bcTypePorosityHi,
	bcTypePermeabilityLo, bcTypePermeabilityHi;

	/// Unified BCs for all scalars (lo side, for each spatial direction)
	Vector<int> bcTypeScalarLo,

	/// Unified BCs for all scalars (hi side, for each spatial direction)
	bcTypeScalarHi;

	/// Vector containing boundary values
	RealVect bcValEnthalpyHi,
	bcValEnthalpyLo,
	bcValBulkConcentrationHi,
	bcValBulkConcentrationLo,
	bcValTemperatureHi,
	bcValTemperatureLo,
	bcValLiquidConcentrationLo,
	bcValLiquidConcentrationHi,
	bcValPorosityLo, bcValPorosityHi,
	bcValPermeabilityLo, bcValPermeabilityHi,
	bcValVelHi,
	bcValVelLo,
	bcValSolidusHi, bcValSolidusLo,
	bcValLiquidusHi, bcValLiquidusLo,
	bcValEutecticHi, bcValEutecticLo,
	bcValSolidConcentrationLo, bcValSolidConcentrationHi;

	/// Time, in case BCs are time-dependent
	Real m_time;

	/// To specify which of timeDependentBCtypes we want to use
	int m_timeDependentBC;

	/// Possible time dependences for BCs
	enum timeDependentBCtypes{
	  m_constant,
	  m_sinusoid, // val = constant val + m_BCamplitude * sin(2 pi t/m_BCtimescale)
	  m_custom
	};

	Real m_BCamplitude;
	Real m_BCtimescale;




	/// 0 means sidewall heating, 1 means vertical heating (sea ice)
	int fixedTempDirection;

	/// Which permeability function should we use?
	int permeabilityFunction;

	/// Are we running the experiment in a Hele-Shaw cell?
	bool heleShaw;

	/// Different possible permeability functions
	enum permeabilityFunctions {
		m_pureFluid,
		m_cubicPermeability,
		m_kozenyCarman,
		m_logPermeability,
		m_permeabilityXSquared,
		m_porosityPermeability
	};

	/// For cases where want to impose a porosity, e.g. for benchmarking
	int m_porosityFunction;
	enum porosityFunctions {
	  m_porosityConstant,
	  m_porosityLinear,
	  m_porosityGaussian,
	  m_porosityEdge,
	  m_porosityTimeDependent,
	};

	/// Different physical problems
	enum physicalProblems {
			m_mushyLayer,
			m_burgersSin,
			m_burgersPeriodic,
			m_poiseuilleFlow,
			m_diffusion,
			m_solidificationNoFlow,
			m_cornerFlow,
			m_sidewallHeating,
			m_HRL,
			m_rayleighBenard,
			m_soluteFluxTest,
			m_refluxTest,
			m_zeroPorosityTest,
			m_meltingIceBlock,
			m_convectionMixedPorous,
			m_vortexPair,
		};

	/// Get all parameters from the inputs file
	void getParameters();

	/// Write out parameters to the command line (pout)
	void printParameters();

	/// Calculate \f$ z(\theta) \f$ for the directional solidification without flow benchmark
	Real directionalSolidificationMushyZ(Real theta, Real zEutectic=1.0);

	/// Read in BCs from the inputs file
	void parseBCs(string a_name, Vector<int>* a_bcHolder, bool required = false);

	/// Read in BC vals from inputs file
	void parseBCVals(string a_name, RealVect& a_bcHolder, bool required = false);

	/// Convert dimensional to non-dimensional temperature
	Real tempTotheta(const Real T);

	/// Convert dimensional to non-dimensional composition
	Real concToTheta(const Real C);

	/// Get the velocity boundary condition
	int getVelBCType(int dir, Side::LoHiSide side);

	/// Compute permeability from porosity
	Real calculatePermeability(Real liquidFraction);

	/// Set the time, might be used for BCs
	void setTime(Real a_time);


};

#endif /* SRC_MUSHYLAYERPARAMS_H_ */
