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
#include "CH_HDF5.H"
#include "BCInfo.h"

// Comment below ensures that enums doxygen generates documentation for enums declared in this file
/*!
@file
Contains immunity flag definitions
*/

/// Different options for enforced porosity
enum ParamsPorosityFunctions {
  /// Sets porosity to it's boundary value on the bottom of the domain, i.e. MushyLayerParams::bcValPorosityLo[0]
  m_porosityConstant,

  /// Porosity varies linearly from the middle to the boundary.
  /**
   * Porosity = 1.0 in the middle, and takes the value given by MushyLayerParams::bcValPorosityLo[0] on all sides.
   * There is a box in the middle of width MushyLayerOptions::fixedPorosityFractionalInnerRadius*domain width which
   * is all at porosity = 1.0.
   */
  m_porosityLinear,

  /// Gaussian profile.
  /**
   * \f$ \chi =  \chi_b + \chi_{max} \exp\left[-((x-x_c)^2 + (y-y_c)^2)/(s L) \right]\f$
   *
   * where \f$ \chi_b \f$ is the porosity enforced at the bottom boundary,
   *
   * \f$ \chi_{max} \f$ is given by MushyLayerOptions::fixedPorosityMaxChi,
   *
   * \f$ (x_c, y_c)\f$ is the centre of the domain.
   *
   * \f$ s \f$, the standard deviation for the gaussian, is given by MushyLayerOptions::FixedPorositySTD
   *
   * \f$ L \f$, the domain width, is given by MushyLayerOptions::domainWidth
   *
   *
   *    */
  m_porosityGaussian,

  /// Gaussian profile but with highest porosity at the edge of the domain, rather than in the middle
  /**
   * \f$ \chi =  \chi_b \left\{ 1 - \chi_{max} \exp\left[-((x-x_c)^2 + (y-y_c)^2)/(s L) \right] \right\} \f$
   *
   * where all the variables are defined as the same way as in m_porosityGaussian
   */
  m_porosityEdge,

  /// Porosity which varies with time
  /**
   * \f$ \chi = 1- 0.3 \alpha (y/H)^2  \sin [\pi x / (2 L)] \f$
   *
   * where \f$ \alpha = \min(4.0, t / \tau) \f$, with the timescale \f$ \tau \f$ defined by MushyLayerOptions::porosityTimescale
   *
   * \f$ H \f$ is the domain height,
   *
   * \f$ L \f$ is the domain width
   */
  m_porosityTimeDependent,
};

/// How does the fluid viscosity depend on the solute concentration (if at all)?
enum ViscosityFunction {
  /// Uniform viscosity, \f$ \mu = 1.0 \f$ everywhere
  uniformViscosity,

  /// Viscosity varies linearly with liquid concentration, from 1 to \f$ \mu_{max} \f$
  /**
   * \f$ \mu = 1 + (\mu_{max} - 1)(\Theta_l + \mathcal{C})/\mathcal{C} \f$
   *
   * where \f$\mu_{max} \f$ is defined by MushyLayerParams::max_viscosity
   */
  linearViscosity,
};


/// Different possible permeability functions
enum PermeabilityFunctions {
  /// Permeability = 1.0
  m_pureFluid,

  /// Cubic permeability function: \f$ \pi = \chi^3 \f$
  m_cubicPermeability,

  /// Kozeny-Carman function: \f$ \Pi = \frac{\chi^3}{(1-\chi)^2} \f$
  m_kozenyCarman,

  /// Logarithmic function: \f$ \Pi = - \chi^2 \log (1-\chi) \f$
  m_logPermeability,

  /// Permeability varies in space. Currently \f$ \Pi = -exp(-(x-0.5)^2/\alpha) \f$ where \f$\alpha=0.02\f$ is some scale.
  m_permeabilityXSquared,

  /// Permeability = porosity, i.e. \f$ \pi = \chi \f$
  m_porosityPermeability
};

/// Different physical problems. Most of these, labelled (Old) should be considered fragile.
enum PhysicalProblems {
  /// Mushy layer - coupled fluid flow, heat, and salt transfer
  m_mushyLayer,

  /// (Old) Burgers equation
  m_burgersSin,

  /// (Old) Burgers equation
  m_burgersPeriodic,

  /// (Old) Poiseuille flow
  m_poiseuilleFlow,

  /// (Old) diffusion (no flow flow)
  m_diffusion,

  /// Solidification with out flow - useful benchmark problem
  m_solidificationNoFlow,

  /// (Old) corner flow
  m_cornerFlow,

  /// (Old) Heating at the sidewalls, with a fixed porosity
  m_sidewallHeating,

  /// (Old) Horton-Rogers-Lapwood problem: heating at the top and bottom in a fixed porous medium
  m_HRL,

  /// (Old) Rayleigh-Benard convection (pure fluid - no porosity effects)
  m_rayleighBenard,

  /// (Old)
  m_soluteFluxTest,

  /// (Old)
  m_refluxTest,

  /// (Old)
  m_zeroPorosityTest,

  /// (Old) Initialise with a block of ice and melt it
  m_meltingIceBlock,

  /// Convection in a porous material with fixed porosity
  m_convectionMixedPorous,

  /// (Old) Vortex merger benchmark problem
  m_vortexPair,
};

/// Class to handle the physical parameters of a mushy layer simulation
/**
 * This includes material properties, boundary conditions, initial conditions and more.
 * Properties are read in and then non-dimensionalised here.
 */
class MushyLayerParams {
public:
        /// Default constructor
	MushyLayerParams();

	/// Default destructor
	virtual ~MushyLayerParams();


	/// Benchmark problem to solve
	PhysicalProblems physicalProblem;

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

	// Temperature at the bottom of the domain
//	bottomTemp,

	// Temperature at the top of the domain
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

	/// How to nondimensionalise the governing equations
	  int m_nondimensionalisation;

	  /// Different options for nondimensionalisation
	  enum nondimensionalisations
	  {
	    /// Diffusive timescale, advective velocity scale
	    m_diffusiveTime_advectiveVel,

	    /// Darcy timescales, advective velocity scale
	    m_darcyTime_advectiveVel,

	    /// Darcy timescale, darcy velocity scale
	    m_darcyTime_darcyVel,

	    /// Advective timescale, darcy velocity scale
	    m_advectiveTime_darcyVel,

	    /// Buoyancy timescale, advective velocity scale
	    m_buoyancyTime_advectiveVel,

	    /// Number of nondimensional schemes
	    m_num_nondimensionalisations
	  };

	  static string s_DARCY_TIMESCALE;
	  static string s_DIFFUSIVE_TIMESCALE;
	  static string s_BUOYANCY_TIMESCALE;
	  static string s_ADVECTIVE_TIMESCALE;

	  static string s_ADVECTIVE_VELOCITY_SCALE;
	  static string s_DARCY_VELOCITY_SCALE;

	  /// Heat diffusion coefficient
	  Real m_heatDiffusionCoeff,

	  /// Salt diffusion coefficient
	  m_saltDiffusionCoeff,

	  /// Viscosity
	  m_viscosityCoeff,

	  /// Buoyancy due to temperature
	  m_buoyancyTCoeff,

	  /// Buoyancy due to liquid concentration
	  m_buoyancySCoeff,

	  /// Darcy coefficient
	  m_darcyCoeff,

	  /// Coefficient for advection terms
	  m_advectionCoeff;


	/// Constant body force, which we add to the buoyancy
	  /**
	   * Buoyancy \f$ = \chi  (Ra_T  \theta - Ra_c*\Theta_l) \mathbf{k} + \mathbf{F}_B \f$
	   *
	   * where \f$ \mathbf{F}_B \f$ is the buoyancy force
	   */
	Real body_force;

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


	/// Dimensionless bulk concentration at the eutectic \f$ \Theta_e \f$
	ThetaEutectic,

	/// Dimensionless initial bulk concentration  \f$ \Theta_i \f$
	ThetaInitial,

	/// Dimensionless far field bulk concentration \f$ \Theta_\infty \f$
	ThetaInf,



	/// Dimensionless initial liquid concentration \f$ \Theta_{l,i} \f$
	ThetaLInitial,



	/// Dimensionless initial solid concentration \f$ \Theta_{s,i} \f$
	ThetaSInitial,



	/// Dimensionless initial enthalpy
	Hinitial,


        /// For plumes


	/// Temperature inflow value in plume
        thetaPlumeInflow,

        /// Enthalpy inflow value in plume
        HPlumeInflow,

        /// Bulk concentration inflow value in plume
        ThetaPlumeInflow,

        /// Liquid concentration inflow value in plume
        ThetaLPlumeInflow,

        /// Solid concentration inflow value in plume
        ThetaSPlumeInflow,

        /// Porosity inflow value in plume
        porosityPlume,

        /// Permeability inflow value in plume
        permeabilityPlume,

        /// Liquidus inflow value in plume
        HLiquidusPlume,

        /// Eutectic inflow value in plume
        HEutecticPlume,

        /// Solidus inflow value in plume
        HSolidusPlume;

	/// Reference temperature for nondimensionalisation
	Real referenceTemperature,

	/// Reference salinity for nondimensionalisation
	referenceSalinity;

	/// Inflow velocity, when required
	Real inflowVelocity;

	/// Pressure difference between top and bottom boundaries
	Real pressureHead;

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

	///Temperature boundary conditions (low side, for each spatial direction)
	bcTypeTemperatureLo,

	/// Temperature boundary conditions (hi side, for each spatial direction)
	bcTypeTemperatureHi,

	/// Liquid concentration boundary conditions (low side, for each spatial direction)
	bcTypeLiquidConcentrationLo,

	/// Liquid concentration boundary conditions (hi side, for each spatial direction)
	bcTypeLiquidConcentrationHi,

	/// Porosity boundary conditions (low side, for each spatial direction)
	bcTypePorosityLo,

	/// Porosity boundary conditions (hi side, for each spatial direction)
	bcTypePorosityHi,

	/// Permeability boundary conditions (low side, for each spatial direction)
	bcTypePermeabilityLo,

	/// Permeability boundary conditions (hi side, for each spatial direction)
	bcTypePermeabilityHi;

	/// Unified BCs for all scalars (lo side, for each spatial direction)
	Vector<int> bcTypeScalarLo,

	/// Unified BCs for all scalars (hi side, for each spatial direction)
	bcTypeScalarHi;

	/// Vector containing boundary values

	/// Enthalpy BCs on hi-side boundaries
	RealVect bcValEnthalpyHi,

	/// Enthalpy BCs on lo-side boundaries
	bcValEnthalpyLo,

	/// Bulk Concentration BCs on hi-side boundaries
	bcValBulkConcentrationHi,

	/// Bulk Concentration BCs on lo-side boundaries
	bcValBulkConcentrationLo,

	/// Temperature BCs on hi-side boundaries
	bcValTemperatureHi,

	/// Temperature BCs on lo-side boundaries
	bcValTemperatureLo,

	/// Liquid concentration BCs on lo-side boundaries
	bcValLiquidConcentrationLo,

	/// Liquid concentration BCs on hi-side boundaries
	bcValLiquidConcentrationHi,

	/// Porosity BCs on lo-side boundaries
	bcValPorosityLo,

	/// Porosity BCs on hi-side boundaries
	bcValPorosityHi,

	/// Permeability BCs on lo-side boundaries
	bcValPermeabilityLo,

	/// Permeability BCs on hi-side boundaries
	bcValPermeabilityHi,

	/// Velocity BCs on hi-side boundaries
	bcValVelHi,

	/// Velocity BCs on lo-side boundaries
	bcValVelLo,

	/// Solidus BCs on hi-side boundaries
	bcValSolidusHi,

	/// Solidus BCs on lo-side boundaries
	bcValSolidusLo,

	/// Liquidus BCs on hi-side boundaries
	bcValLiquidusHi,

	/// Liquidus BCs on lo-side boundaries
	bcValLiquidusLo,

	/// Eutectic BCs on hi-side boundaries
	bcValEutecticHi,

	/// Eutectic BCs on lo-side boundaries
	bcValEutecticLo,

	/// Solid concentration BCs on lo-side boundaries
	bcValSolidConcentrationLo,

	/// Solid concentration BCs on hi-side boundaries
	bcValSolidConcentrationHi,

	/// Pressure BCs on hi/lo sides.
	/**
	 * Note that the type of pressure bc is determine from the
	 * type of velocity bc specified, as the two are coupled.
	 * The values specified here are only used for enforcing a pressure head
	 */
	bcValPressureHi,

	/// Pressure BC on the low side of the domain in each direction (x, y, z)
	bcValPressureLo;


	BCInfo
	/// Position along each face, below which to enforce no heat flux (above this position, enforce some other specified heat flux)
	m_bc_noFluxLimit,

	/// Reference temperature for radiation boundary conditions
	m_bc_bTref,

	/// Some generic variable which can be used in various boundary conditions if required
	m_bc_b;

	int max_bc_iter;

	/// First or second order BCs
	int m_BCAccuracy;

	/// Accuracy for pressure bcs
	int m_pressureBCAccuracy;

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


	/// Amplitude of time periodic boundary conditions
	Real m_BCamplitude;

	/// Time scale (period) for time periodic boundary conditions
	Real m_BCtimescale;




	/// 0 means sidewall heating, 1 means vertical heating (sea ice)
	int fixedTempDirection;

	/// Which permeability function should we use?
	PermeabilityFunctions permeabilityFunction;

	/// Are we running the experiment in a Hele-Shaw cell?
	bool heleShaw;

	/// For cases where want to impose a porosity, e.g. for benchmarking
	ParamsPorosityFunctions m_porosityFunction;

	/// For when we want the fluid viscosity to depend on the solute concentration
	ViscosityFunction m_viscosityFunction;

	/// Maximum dimensionless viscosity the fluid can have
	/**
	 * Only used when m_viscosityFunction is not = 0
	 */
	Real max_viscosity;

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

	/// Write out the key parameters to hdf5 file, so that solution can be reconstructed
	void writeToHDF5(HDF5HeaderData& a_header) const;

	/// Returns a string indicating the velocity scale used for nondimensionalisation
	string getVelocityScale() const;

	/// Utility function to compute porosity given a single enthalpy/bulk concentration value
	Real computePorosity(Real H, Real C);

	Real compute_dHdT(Real H, Real C);

	/// Utility function to compute temperature given a single enthalpy/bulk concentration value
	Real computeTemperature(Real H, Real C);

	/// Compute everything from enthalpy and bulk concentration
	//todo - implement
	void computeDiagnosticVars(Real H, Real C, Real T, Real porosity, Real Cl, Real Cs);

	/// Returns a string indicating the time scale used for nondimensionalisation
	string getTimescale() const;

	/// Returns whether or not we're solving the Darcy-Brinkman equation
	bool isDarcyBrinkman();

	/// Returns whether or not the problem we're solving contains viscous terms
	bool isViscous();

	/// Compute boundary conditions for fields which can be found from the enthalpy and salinity via the phase diagram
	void computeDerivedBCs ();
};

#endif /* SRC_MUSHYLAYERPARAMS_H_ */
