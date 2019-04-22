#ifndef PHYSBCUTIL_H
#define PHYSBCUTIL_H

#include "Tuple.H"
#include "BCFunc.H"
#include "MushyLayerParams.h"
#include "PhysIBC.H"
#include "AdvectIBC.H"



/**
 * Function interface for applying enthalpy boundary conditions
 * which satisfy some condition based on the temperature field.
 * Computes some residual, which is characterises the extent
 * to which we have not satisfied the temperature constraint, which we attempt to minimise.
 */
//typedef void(*EnthalpyBCResidual)(Real&   a_residual,
//                      const Real  a_dx,
//                      const Real a_ghost_
//                      Real                 a_dx,
//                      bool                 a_homogeneous);

/// Class for implementing temperature bcs on the enthalpy field
/**
 * This should be an abstract base class - different applications should inherit from this
 */
class NonlinearTemperatureBC
{
public:
  NonlinearTemperatureBC();

  // destructor
  virtual ~NonlinearTemperatureBC()
  {

  }

  NonlinearTemperatureBC(const MushyLayerParams a_params,
                         const Real a_dx)
  {
    CH_TIME("NonlinearTemperatureBC::NonlinearTemperatureBC");

    m_params = a_params;
    m_dx = a_dx;
  }

  /// Compute the residual given the current estimate of the enthalpy on the boundary
  virtual void computeResidual(Real& a_residual, const Real a_ghost_enthalpy, const Real a_ghost_bulk_concentration, const Real a_interior_temperature)
  {

  }

protected:

  MushyLayerParams m_params;
  Real m_dx;

};

/// Robin temperature bcs, for applying to the enthalpy field
/**
 * Compute the residual of
 *   \f$ a \frac{dT}{dz} - F_0 - b(T - T_{ref}) \f$
 *
 *  Specifically, the residual \f$ r \f$ is
 *
 *  \f$ r = \Delta x \left[F_0 - b( T_b - T_{ref}) \right] - (T_g - T_i) \f$
 *
 *  where \f$ T_b = (T_g + T_i)/2 \f$ is the temperature on the boundary,
 *  \f$ T_g \f$ is the temperature in the ghost cell and \f$ T_i \f$ is the temperature
 *  on the interior cell.
 *
 */
class NonlinearTemperatureBCRobin : public NonlinearTemperatureBC
{
public:

  virtual ~NonlinearTemperatureBCRobin()
  {

  }

  NonlinearTemperatureBCRobin();

  NonlinearTemperatureBCRobin(const MushyLayerParams a_params,
                              const Real a_dx,
                              const Real a,
                              const Real b,
                              const Real F,
                              const Real T_ref) : NonlinearTemperatureBC(a_params, a_dx)
  {
    m_a = a;
    m_b = b;
    m_F = F;
    m_T_ref = T_ref;
  }

  virtual void computeResidual(Real& a_residual, const Real a_ghost_enthalpy, const Real a_ghost_bulk_concentration, const Real a_interior_temperature)
  {
    Real a_ghost_T = m_params.computeTemperature(a_ghost_enthalpy, a_ghost_bulk_concentration);

    Real Tflux = (a_ghost_T - a_interior_temperature);
    Real Tboundary = (a_ghost_T + a_interior_temperature)/2;

    a_residual = m_dx * (m_F - m_b*(Tboundary - m_T_ref)) - Tflux;
  }


protected:
  Real m_a, m_b, m_F, m_T_ref;
};

/// Big class to encapsulate physical boundary conditions
class PhysBCUtil
{

public:
	/** this will access the ParmParse database using the
      physBC prefix and use it to define the physical BC types (hi
      and lo) in each direction.  The BC type "custom" will imply
      that the BC is a somewhat nonstandard one, but that the
      derived class will know what to do with it. */
	PhysBCUtil(MushyLayerParams a_params, Real a_dx);

	PhysBCUtil();

	///
	virtual ~PhysBCUtil();
	/// Different types of velocity boundary conditions
	/**
	 * SolidWall - no normal flow, no transverse flow is viscous
	 * Inflow -
	 * Outflow -
	 * OutflowNormal - normal component has zero gradient at boundary, transverse component = 0
	 * VelInflowOutflow -
	 * noShear -
	 * Symmetry -
	 */
	enum VelBCType
	{
		///
		bogusBC = -1,
		///
		SolidWall,
		Inflow,
		Outflow,
		OutflowNormal, // only a normal velocity
		VelInflowOutflow, // both inflow and outflow possible
		noShear,
		Symmetry,
		VelInflowPlume,
		OutflowPressureGrad,
		PressureHead,
		NUM_PHYS_BC_TYPES
	};

	/// Different options for scalar boundary conditions
	enum scalarBCType
	{
	  /// Fix value at boundary
	  Dirichlet,
	  /// Fix gradient at boundary
	  Neumann,
	  /// Separate conditions for inflow and outflow regions
	  /**
	   * Where there's inflow - dirichlet,
	   * where there's outflow - zero normal gradient
	   */
	  InflowOutflow,
	  /// Allow inflow but no outflow
	  /**
	   * For e.g. a plume forced into the domain.
	   * Enforce dirichlet BCs everywhere, but with a different value for
	   * inflow vs outflow regions
	   */
	  OnlyInflow,
	  /// Robin
	  /**
	   * Not implemented yet
	   */
	  Robin,
	  /// Variable flux
	  /**
	   * Flux which varies along the boundary
	   */
	  VariableFlux,
	  /// Fixed temperature
	  /**
	   * This is a special case for the heat equation, where we apply BCs to enthalpy
	   * and need to compute the boundary from the interior porosity and specified temperature
	   */
	  FixedTemperature,
	  /// Temperature Flux
	  /**
	   * BCs for enthalpy. Assume bulk concentration \f$ C_g \f$ in the ghost cell is equal to
	   * it's value in the interior of the domain (i.e. assume a no normal flux condition) then
	   * iteratively estimate the enthalpy \f$ H_g \f$ such that
	   * \f$ \frac{T_g(H_g,C_g) - T_{interior} }{\Delta x}= F \f$,
	   * in order to satisfy the temperature flux condition. Finally, set the enthalpy ghost value to \f$ H_g \f$.
	   *
	   * Only do this for \f$z>0.75\f$, otherwise, no flux.
	   *
	   */
	  TemperatureFlux,

	  /// Temperature flux + radiation
	  TemperatureFluxRadiation
	};

	/// Physical parameters for the system
	MushyLayerParams m_params;

	/// "virtual" constructor
	virtual PhysBCUtil* newPhysBCUtil() const;

	/// returns velocity BC for tracing (inviscid)
	virtual Tuple<BCHolder, SpaceDim> tracingVelFuncBC() const;

	/// returns edge-centered vel BC for edge-ctrd vel (single component)
	virtual Tuple<BCHolder, SpaceDim> advectionVelFuncBC(bool a_isViscous) const;

	/// Currently unused
	virtual Tuple<BCHolder, SpaceDim> fluxExtrapBC() const;

	/// BC holder for \f$ \mathbf{u}^* \f$
	virtual Tuple<BCHolder, SpaceDim> uStarFuncBC(bool a_isViscous = false) const;


	/// BC holder for \f$ \mathbf{u} \f$ (edge centred)
	virtual Tuple<BCHolder, SpaceDim> edgeVelFuncBC(bool a_isViscous = false,
	                                                LevelData<FluxBox>* velocityBCVals = NULL) const;


	/// BC holder for porosity \f$ \chi \f$, (edge centred)
	virtual Tuple<BCHolder, SpaceDim> porosityFaceBC() const;

	/// BCs for permeability \f$ \Pi \f$ (edge centred)
	virtual Tuple<BCHolder, SpaceDim> permeabilityFaceBC() const;

	/// BC holder for velocity \f$ \mathbf{u} \f$ - extrapolation
	virtual Tuple<BCHolder, SpaceDim> velExtrapBC(Interval ghostInterval=Interval(0,0)) const;

	/// returns single-component BC for viscous refluxing solves
	virtual BCHolder viscousRefluxBC(int a_dir, bool a_isViscous=false) const;

	/// returns a BC object compatible with AMRGodunov advection infrastructure
	virtual PhysIBC* advectionVelIBC() const;

	/// Is object defined?
	bool isDefined();

	/// sets cell spacing
	void Dx(const Real a_dx);

	/// Convert to AdvectIBC BCs
	void convertBCType(const int a_implicitBC,  const Real a_implicitVal,
                           int& a_explicitBC, Real& a_explicitVal) const;

	/// Handy to have a function of this form, too
	int convertBCType(const int a_implicitBC) const;

	/// returns cell spacing
	Real Dx() const;

	/// sets time
	void Time(const Real a_time);

	/// returns current time
	Real Time() const;

	/// Set pointer to advection velocity
	void setAdvVel(LevelData<FluxBox>* a_advVel);

	/// Get pointer to advection velocity
	LevelData<FluxBox>* getAdvVel();

	/// BC holder for pressure \f$ P \f$
	virtual BCHolder LevelPressureFuncBC() const;

	/// BC holder for synchronisation correction
	/**
	 * Just use same BC as \f$ \nabla P \f$
	 */
	virtual BCHolder  gradESyncFuncBC() const;

	/// BC holder for \f$ \nabla \Pi \f$
	/*
	 * Just use same BC as \f$ \nabla P \f$
	 */
	virtual BCHolder  gradPiFuncBC() const;

	/// Extrap BCs on all domain boundaries
	virtual BCHolder extrapFuncBC() const;

	 // BC on stream function (used when computing stream func. from vorticity)
//	  virtual BCFunc streamBC() const;

	/// BC holder for freestream correction
	/**
	 * Normal gradient  = 0
	 */
	virtual BCHolder FreestreamCorrFuncBC() const;


	/// BC holder for synchronisation projection
	/**
	 * Same BCs as for pressure
	 */
	virtual BCHolder SyncProjFuncBC() const;

	/// BC holder for \f$ \nabla P \f$ in the MAC projection
	/**
	 * Same as for \f$ \nabla P \f$ in general
	 */
	virtual BCHolder gradMacPressureFuncBC() const;

	/// BC holder for freestream preservation correction \f$ \nabla e_\lambda \f$
	        /**
	         * Same as for \f$ \nabla P \f$ in general
	         */
	virtual BCHolder gradELambdaFuncBC() const;

	/// returns single-component BC for viscous solves
	/**
	 * Applies BCs to the component in the Interval(0,0)
	 */
	virtual BCHolder velFuncBC(int a_dir, bool a_viscous, Interval interval=Interval(0,0) ) const;

	/// returns single-component BC
	/**
	 * Applies BCs to component in the Interval(a_dir, a_dir)
	 */
//	virtual BCHolder velFuncMultiCompBC(int a_dir) const;

	/// returns multi-component BC for applying to fields which have been calculated
	/**
	 * Not for use when computing fields - use velFuncBC() for solves.
	 * 0th order extrapolation, also fills interior ghost cells
	 */
	virtual BCHolder extrapolationFuncBC(int a_order=0) const;

	/// BC holder for \f$ \nabla P \f$
	/**
	 * Extrapolation BCs on all boundaries
	 */
	virtual BCHolder BasicGradPressureFuncBC() const;

	/// BCs for streamfunction
	virtual BCHolder streamFunctionBC() const;

	/// BC holder for pressure
	virtual BCHolder BasicPressureFuncBC(bool a_isHomogeneous) const;

//	virtual BCHolder BasicPressureFuncBCNonsubcycle(bool a_isHomogeneous,
//			Vector<RefCountedPtr<LevelData<FArrayBox> > > a_theta =  Vector<RefCountedPtr<LevelData<FArrayBox> > >(0)) const;

	/// BC holder for liquid concentration, \f$ \Theta_l \f$
	virtual BCHolder ThetaLFuncBC(bool a_homogeneous = false, LevelData<FluxBox>* a_advVel= NULL) const;


	/// BC holder for temperature, \f$ \theta \f$
	virtual BCHolder BasicthetaFuncBC(bool a_homogeneous = false, LevelData<FluxBox>* a_advVel = NULL) const;

//	virtual BCHolder porosityFuncBC(bool a_homogeneous = false) const;

	/// BC holder for porosity, \f$ \chi \f$ (cell centred)
	virtual BCHolder BasicPorosityFuncBC(bool a_homogeneous = false) const;

	/// BC holder for porosity, \f$ \chi \f$ (edge centred)
	virtual BCHolder BasicPorosityFaceFuncBC(int a_comp, bool a_homogeneous = false) const;

	/// BC holder for permeability, \f$ \Pi \f$
	virtual BCHolder BasicPermeabilityFuncBC(bool a_homogeneous = false) const;

	/// BC holder for permeability, \f$ \Pi \f$ (edge centred)
	virtual BCHolder BasicPermeabilityFaceFuncBC(int a_comp, bool a_homogeneous = false) const;

	/// BC holder for enthalpy
	virtual BCHolder BasicEnthalpyFuncBC(bool a_homogeneous = false,
	                                     LevelData<FluxBox>* a_advVel = NULL) const;

	/// Boundary conditions for calculating the liquidus
	virtual BCHolder BasicLiquidusBC(bool a_homogeneous = false,
	                                             LevelData<FluxBox>* a_advVel = NULL) const;

	/// Boundary conditions for calculating the eutectic
	virtual BCHolder BasicEutecticBC(bool a_homogeneous = false,
	                                             LevelData<FluxBox>* a_advVel = NULL) const;

	/// Boundary conditions for calculating the solidus
	virtual BCHolder BasicSolidusBC(bool a_homogeneous = false,
	                                             LevelData<FluxBox>* a_advVel = NULL) const;

	/// BC holder for some generic scalar
	virtual BCHolder BasicScalarFuncBC(Real a_topVal, Real a_bottomVal,
	                                   Real a_plumeVal,
	                                   bool a_homogeneous=false) const;

	/// BC holder for bulk concentration \f$ \Theta \f$
	virtual BCHolder ThetaFuncBC(bool a_homogeneous = false, LevelData<FluxBox>* a_advVel = NULL) const;

	/// BC holder for \f$ \lambda \f$ (for finding freestream preservation errors)
	virtual BCHolder lambdaBC(bool a_homogeneous = false, bool a_scaleWithPorosity=false)const;

	/// BCs for no flux boundaries
	virtual BCHolder noFluxBC() const;

	/// Initial and Boundary Conditions for advecting scalars
	virtual PhysIBC* scalarTraceIBC(RealVect bcValLo, RealVect bcValHi,
	                                IntVect bcTypeLo, IntVect bcTypeHi,
	                                Real a_bcPlumeVal = 0.0, Vector<Real> a_plumeBounds = Vector<Real>(0)) const;

	/// For multi component enthalpy-salinity solve
	virtual PhysIBC* scalarTraceHC_IBC() const;

	/// Boundary conditions for Enthalpy during advection solves
	virtual PhysIBC* scalarTraceH_IBC() const;

	/// Boundary conditions for Bulk Concentration during advection solves
	virtual PhysIBC* scalarTraceC_IBC() const;

	/// Generic boundary conditions for Bulk Concentration during advection solves
	virtual PhysIBC*  scalarTrace_IBC(Vector<int> a_bcTypeHi, Vector<int> a_bcTypeLo,
	                                  RealVect a_bcValHi, RealVect a_bcValLo,
	                                  Real a_plumeInflow) const;

	/// Apply frame advection BCs if appropriate
	void applyFrameAdvectionBC(IntVect& bcTypeHi, IntVect& bcTypeLo) const;

	/// For multi component enthalpy-salinity solve
	virtual PhysIBC* scalarTraceTSl_IBC() const;


	/// BC for object with two components: temperature and liquid salinity
	virtual BCHolder temperatureLiquidSalinityBC(bool a_homogeneous = false) const;

	/// BC for object with two components: enthalpy and bulk salinity
	virtual BCHolder enthalpySalinityBC(bool a_homogeneous = false) const;

	/// Apply extrapolation BCs to velocity
	        virtual BCHolder extrapVelFuncBC( int  a_comp,
	                        Interval ghostInterval = Interval(0,0)) const;

protected:


	/// BCHolder for edge centred velocity
	virtual BCHolder basicECVelFuncBC(bool a_isHomogeneous,
			bool a_isViscous,
			int  a_comp,
			const Interval& a_interval,
			LevelData<FluxBox>* a_velocityBCVals = NULL) const;



	/// returns a BCHolder for cell-centered velocity
	virtual BCHolder basicCCVelFuncBC(bool a_isHomogeneous,
			bool a_isViscous,
			int a_velComponent,
			const Interval& a_interval) const;

//	virtual BCHolder uStarCCVelFuncBC(bool a_isHomogeneous,
//			bool a_isViscous,
//			int  a_comp,
//			const Interval& a_interval) const;



	/// Apply extrapolation BCs to velocity (interior ghost cells)
	virtual BCHolder extrapVelFuncInterior(bool a_isHomogeneous,
			int  a_comp,
			const Interval& a_interval,
			int n_ghost) const;

	///  interact with ParmParse to set physical BC types
	virtual void setBCs();

	/// contains the enumerated values for the lo-end BC
	Tuple<int, SpaceDim> m_loBC;

	/// contains the enumerated values for the hi-end BC
	Tuple<int, SpaceDim> m_hiBC;

private:
	/// cell spacing
	Real m_dx;

	/// current time
	Real m_time;

	/// Is this object defined?
	bool m_defined;

	/// Pointer to advection velocity
	LevelData<FluxBox>* m_advVel;

	/// copy constructor
	PhysBCUtil(const PhysBCUtil& rhs);

	///
	virtual PhysBCUtil& operator=(const PhysBCUtil& rhs);

};

#endif
