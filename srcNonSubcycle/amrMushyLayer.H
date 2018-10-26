/* _______              __
  / ___/ /  ___  __ _  / /  ___
 / /__/ _ \/ _ \/  ' \/ _ \/ _ \
 \___/_//_/\___/_/_/_/_.__/\___/
 */
//
// This software is copyright (C) by the Lawrence Berkeley
// National Laboratory.  Permission is granted to reproduce
// this software for non-commercial purposes provided that
// this notice is left intact.
//
// It is acknowledged that the U.S. Government has rights to
// this software under Contract DE-AC03-765F00098 between
// the U.S.  Department of Energy and the University of
// California.
//
// This software is provided as a professional and academic
// contribution for joint exchange. Thus it is experimental,
// is provided ``as is'', with no warranties of any kind
// whatsoever, no support, no promise of updates, or printed
// documentation. By using this software, you acknowledge
// that the Lawrence Berkeley National Laboratory and
// Regents of the University of California shall have no
// liability with respect to the infringement of other
// copyrights by any part of this software.
//

#ifndef AMR_MUSHYLAYER_H
#define AMR_MUSHYLAYER_H

#include "LevelData.H"
#include "FArrayBox.H"
#include <assert.h>

#include "../src/MushyLayerUtils.H"
#include "RealVect.H"
#include "BackwardEuler.H"
#include "AMRPoissonOp.H"
#include "LevelAdvect.H"
#include "BiCGStabSolver.H"
#include "LevelTGA.H"
#include "VCAMRPoissonOp2.H"
#include "CCProjectorComp.H"
//#include "PhysBCUtil.H"
#include "LevelFluxRegisterEdge.H"


#include "MushyLayerParams.h"

/// class to manage non-subcycled AMR computation of mushy layer dynamics
/**
 */
class amrMushyLayer
{

public:
	/// Default constructor
	/** At the moment, there is only one constructor, which defines itself
      based on an inputs file.  (This may change, of course)
	 */
	amrMushyLayer();

	/// destructor
	/** destructor
	 */
	~amrMushyLayer();

	/// set default values before initialization
	void setDefaults();

	/// initializes object based on inputs data passed in through parmParse
	/** initializes new object based on data in ParmParse database.  This
      may change to a more explicit initialization if it's thought to be
      a good idea.
	 */
	void initialize();


	/// advance solution until either max_time or max_step are reached
	/**
	 */
	void run(Real a_max_time, int a_max_step);

	/// is this object defined and initialized?
	/**
	 */
	bool isDefined() const;

#ifdef CH_USE_HDF5

	/// write hdf5 plotfile to the standard location
	/** writes plotfile to the standard location; filename is
      <plot_root>.<step>.<DIM>d.hdf5
	 */
	void writePlotFile();

	//Write plotfile at a specific iteration
	void writePlotFile(int iteration);

	void writePlotHeader(HDF5Handle& a_handle) const;


	/// write checkpoint file out for later restarting
	/** writes checkpoint data to the standard location; filename is
      <check_root>.<step>.<DIM>d.hdf5
	 */
	void writeCheckpointFile() const;


	/// read checkpoint file for restart
	/** read checkpoint data from file pointed to by a_handle
	 */
	void readCheckpointFile(HDF5Handle& a_handle);

	/// set up for restart
	void restart(string& a_restart_file);

#endif

protected:
	/// compute one timestep
	void timeStep();

	/// do regridding
	void regrid();

	/// compute tags for regridding
	void tagCells(Vector<IntVectSet>& a_tags);

	/// compute tags for the level a_level
	void tagCellsLevel(IntVectSet& a_tags, int a_level);

	/// compute tags at initial time
	void tagCellsInit(Vector<IntVectSet>& a_tags);

	/// Initialize all global variables on this level
	void initVars(const int lev);

	/// initialize grids at initial time
	void initGrids(int a_finest_level);

	/// set up grids from grids file
	void setupFixedGrids(const std::string& a_gridFile);

	/// initialize data on hierarchy
	void initData();

	/// Fill variable with data from file
	void loadFromFile(const int a_var, const string filename, const int num_ghost, const int comp=0);

	/// Calculate time dependent porosity (for some benchmarks)
	Real timeDepPorosity(Real x, Real z, Real t);

	/// compute timestep
	void computeDt();

	/// compute timestep at initial time
	void computeInitialDt();

	/// Check if we've reached steady state
	void convergedToSteadyState();

	/// Initialise elliptic operators for \f$ \theta, \Theta_l \f$.
	void setupPoissonSolvers(
			Vector<LevelData<FArrayBox>* >& a_thetaNew,
			Vector<LevelData<FArrayBox>* >& thetaSource,
			Vector<LevelData<FArrayBox>* >& a_ThetaLNew,
			Vector<LevelData<FArrayBox>* >& ThetaSource,
			Vector<DisjointBoxLayout>& activeGrids);

	/// Setup multigrid solvers and timestepping (backward euler, TGA etc.)
	void setupMultigrid(
			Vector<LevelData<FArrayBox>* >& a_thetaNew,
			Vector<LevelData<FArrayBox>* >& thetaSource,
			Vector<LevelData<FArrayBox>* >& a_ThetaLNew,
			Vector<LevelData<FArrayBox>* >& ThetaLSource,
			Vector<DisjointBoxLayout>& activeGrids);

	/// Initialise Godunov operators for frame advection calculations (enthalpy and bulk concentration)
	void setupAdvectionSolvers();

	/// Initialise various initial and boundary conditions.
	void setupIBCS();

	/// Get grids and scalar fields for the levels that we're currently working on. \n
	/// May well be less than the maximum number of levels.
	void getActiveGridsVars(Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > >& a_scalarNew,
			Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > >& a_scalarOld,
			Vector<DisjointBoxLayout>&a_amrGrids);

	/// Get all the variables needed to calculate advection on this level
	void getLevelAdvectionsVars(const int a_var,
			const int lev,
			int& refRatio,
			const bool hasCoarser,
			const bool hasFiner,
			const bool use_limiting,
			RefCountedPtr<LevelData<FArrayBox> >& coarserDataOldPtrEnthalpy,
			RefCountedPtr<LevelData<FArrayBox> >& coarserDataNewPtrEnthalpy,
			RefCountedPtr<LevelFluxRegister>& coarserFRPtrEnthalpy,
			RefCountedPtr<LevelFluxRegister>& finerFRPtrEnthalpy,
			DisjointBoxLayout& coarseGrid,
			LevelAdvect& levAdvect,
			Real& tCoarserNew,
			Real& tCoarserOld);

	/// Populate m_frameAdv with the frame advection velocity \f$ \mathbf{V} = - |V| \mathbf{e}_z \f$
	/// which is constant in space and time
	void getFrameAdvection(const int lev);

	/// Convert the FArrayBox of fluid velocities to a fluxbox
	void getFluidAdvection(const int lev);

	/// Get coefficients for the \f$ \theta \f$ VCAMRPoissonOp2 object.
	/**
	  aCoef = \f$ \chi + (1-\chi) c_p \f$ \n
	 bCoef = \f$ \chi + (1-\chi) k \f$
	 */
	void thetaVCOpCoeffs(Real& alpha, Real& beta,
			Vector<RefCountedPtr<LevelData<FArrayBox> > >& aCoef,
			Vector<RefCountedPtr<LevelData<FluxBox> > >& bCoef);

	/// Advect the scalar field by one timestep.
	/// The source term is used to calculate intermediate fluxes, but not included in the final calculation.
	void advectScalar(int a_var,
			Vector<LevelData<FArrayBox>* > a_source,
			Vector<LevelData<FluxBox>* > a_advVel);

	/// reflux
	void doReflux(int a_var, int refluxLevel);

	/// Reflux for when we have diffusion as well. Ripped from AMRGodunov example
	void doImplicitReflux(int a_var, int lev);

	/// Advect scalar and put U.grad(var) into a_advTerm
	void makeAdvectionTerm(const int a_var,
			const Vector<LevelData<FArrayBox>* > a_source,
			const Vector<LevelData<FluxBox>* > a_advVel,
			Vector<RefCountedPtr<LevelData<FArrayBox> > > a_advTerm);

	/// Copy across all the active levels of a Vector<LevelData> to another
	void activeLevelCopy(Vector<RefCountedPtr<LevelData<FArrayBox> > > from,
			Vector<RefCountedPtr<LevelData<FArrayBox> > > to);

	///Calculate the source term for the \f$ \theta \f$ equation, \n
	/// \f$ B \theta_\infty + V \frac{\partial H^{n+1/2}}{\partial z} +
	/// \left(\theta^{n+1/2} (c_p - 1) - S \right) \frac{\chi^* - \chi^n}{\Delta t}\f$
	void calculatethetaSourceTerm(Vector<LevelData<FArrayBox>* > thetaSource,
			Vector<LevelData<FArrayBox>* > V_gradH);

	///Calculate \f$ H = \chi S + (\chi + (1-\chi) c_p)\theta \f$
	void calculateEnthalpy();

	///Calculate \f$ \Theta = \chi \Theta_l + (1-\chi) \Theta_s\f$
	void calculateBulkConcentration();

	///Calculate porosity using the phase diagram
	void calculatePorosity();

	//	void correctThetaL();

	///Calculate \f$ \nabla^2 \phi \f$ on this level if we have an operator for that
	void makeLaplacian(const int a_var, const int lev, LevelData<FArrayBox>& a_laplacian);

	///Update enthalpy method variables
	void updateEnthalpyVariables();

	//	Real calculateHeatEqnResidual();

	///Calculate \f$ H^{k} - H^{k+1} \f$ where \f$ k \f$ is the iteration number
	Real calculateHeatEqnResidual2(const Vector<LevelData<FArrayBox>* >a_thetaNew, const Vector<LevelData<FArrayBox>* > a_thetaNewPrev);

	//	Real calculateConcEqnResidual();

	/// Calculate the analytic solution for this particular problem if possible
	void calculateAnalyticSolution();

	/// Analytic solution for benchmark 1 - directional solidification without flow (Worster 1991)
	void calculateSolnBm1();

	/// Analytic solution for simple temperature diffusion \f$ \partial \theta / \partial t = \nabla^2 \theta \f$
	void calculateSolnDiffusion();

	/// Analytic solution for benchmark 2 - flow in a fixed porous medium
	void calculateSolnBm2();

	/// Fill the fields used in calculations with the analytic values, if possible
	void enforceAnalyticSolution();

	/// Write the primary variables of the current solution to a text file
	void writeSolutionToTextFile(int lev);

	/// replace coarse data with the average of finer data where appropriate
	void averageCoarseToFineSolutions();

	/// Fill fine levels with data from coarser levels
	void interpToFineSolutions(int lev,
			const std::vector<int> scalarVars = std::vector<int>(),
			const std::vector<int> vectorVars = std::vector<int>());

	/// Check a sequence of residuals isn't diverging
	bool residualsDiverging(const Vector<Real> residuals);

	/// Check that a sequence of residuals isn't equal (relative to our tolerance level)
	bool residualsConstant(const Vector<Real> residuals);

	/// Set up flux registers
	void setupFluxRegisters();

	/// A test to see what CF interp does
	void doCFInterpTest();

	/// Calculate \f$ V d \phi / dt \f$ using a one sided, second order, finite difference
	void calculateFrameAdvection(int a_var, Vector<LevelData<FArrayBox>* >& a_grad, Real frameAdvVel);

	/// Calculate \f$ V d \phi / dt \f$ using a Godunov (flux conservative) method. Think this is first order in space.
	void calculateFrameAdvectionGodunov(const int a_var);

	/// Calculate the source term for the \f$ \Theta_l \f$ equation,
	/// \f$ V \frac{\partial \Theta^{n+1/2}}{\partial z} - \frac{(1-\chi^*)\Theta_s^* - (1-\chi^n)\Theta_s^n}{\Delta t} - \Theta_l^{n+1/2} \frac{\chi^* - \chi^n}{\Delta t} \f$
	/// where superscripts denote time centering
	void calculateThetaLSourceTerm(Vector<LevelData<FArrayBox>* > V_dThetadz_n,
			Vector<LevelData<FArrayBox>* > ThetaDiffusion_n);

	/// Calculate \f$ Le^{-1} \nabla \chi \nabla \Theta_l \f$
	void calculateThetaDiffusion(Vector<LevelData<FArrayBox>* > a_ThetaDiffusion);

	/// Calculate \f$ \nabla^2 \theta \f$
	void calculatethetaDiffusion(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_ThetaDiffusion,
			Vector<LevelData<FArrayBox>* >& a_thetaOld);

	///calculate the concentration equations residual Theta^{k} - Theta^{k-1} where k is iteration number.
	/// Returned value can be positive or negative!
	Real calculateResidual2(const Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phiOld,
			const Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phiNew);

	/// Calculate the Darcy velocity (mass flux) based on the current temperature and porosity fields
	/// Alpha determines which timestep to use for the scalar fields.
	/// alpha = 0 - old time, alpha = 1 - new time, alpha = 1/2 - mean of both fields
	void updateVelocity(Vector<DisjointBoxLayout> activeGrids,  Real alpha=0);

	/// Not used anymore
	void updateVelocityComp(Vector<DisjointBoxLayout> activeGrids,  Real alpha=0);

	/// Do implicity velocity reflux on level
	void doImplicitVelocityReflux(int m_level);

	/// Do explicit velocity reflux
	void doExplicitVelocityReflux();

	/// Initialise all scalar fields on this level, so the data structures can be used
	void initScalarVars(const int lev);

	/// Initialise all vector fields on this level, so the data structures can be used
	void initVectorVars(const int lev);

	/// Calculate the average of m_scalarNew[a_var] and m_scalarOld[a_var], i.e.
	/// \f$ (\phi^n + \phi^{n-1}) / 2 \f$
	Vector<LevelData<FArrayBox>*> timeCenteredScalar(const int a_var);

	//This doesn't seem to work
	//	void initialiseVarToZero(Vector<LevelData<FArrayBox>* > a_phi);

	/// For each value x in a_phi, if |x| > |a_cap| set x=a_cap (with the relevant sign)
	void enforceAbsCap(Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phi, Real a_cap);

	///Remove nan values and replace them with 0
	void removeNanValues(FArrayBox &fab);

	/// Do one iteration of the nonlinear solver
	bool nonlinearSolver(Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > a_scalarNew,
			Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > >a_scalarOld,
			Vector<DisjointBoxLayout> a_amrGrids,
			Vector<LevelData<FArrayBox>* >  a_thetaNew,
			Vector<LevelData<FArrayBox>* > a_thetaOld,
			Vector<LevelData<FArrayBox>* > a_thetaDiff,
			Vector<LevelData<FArrayBox>* > a_thetaOldTemp,
			Vector<LevelData<FArrayBox>* > a_ThetaLNew,
			Vector<LevelData<FArrayBox>* > a_ThetaLOld,
			Vector<LevelData<FArrayBox>* > a_thetaNewPrev,
			Vector<LevelData<FArrayBox>* > zeroSource,
			Vector<LevelData<FArrayBox>* > V_gradH_n,
			Vector<LevelData<FArrayBox>* > V_dThetadz_n,
			Vector<LevelData<FArrayBox>* > thetaSource,
			Vector<LevelData<FArrayBox>* > ThetaLSource,
			Vector<RefCountedPtr<LevelData<FArrayBox> > > a_ThetaPrevIteration,
			Vector<RefCountedPtr<LevelData<FArrayBox> > > a_HPrevIteration,
			Vector<RefCountedPtr<LevelData<FArrayBox> > > a_ThetaTwoPrev,
			Vector<RefCountedPtr<LevelData<FArrayBox> > > a_HTwoPrev,
			Real residualEnthalpy,
			Real residualConc,
			Vector<Real> &heatResiduals, Vector<Real> &concResiduals,
			int &iteration);




	/// Get the coefficients for the \f$ \Theta_l \f$ elliptic operator, VCAMRPoissonOp2.
	///  (alpha * aCoef(x) * I - beta * Div(bCoef(x) . Grad)) phi = rho \n
	/// alpha, beta are set later by the BackwardEuler/TGA solver. \n
	/// aCoef(x) = \f$ (\chi^{n+1} + \chi^{n})/2 \f$. \n
	/// bCoef(x) = \f$ Le^{-1} \chi^{n+1} \f$ \n
	void getThetaLOpCoeffs(Real& alpha, Real& beta,
			Vector<RefCountedPtr<LevelData<FArrayBox> > >& aCoef,
			Vector<RefCountedPtr<LevelData<FluxBox> > >& bCoef);

	/// helper function - if priority >= s_verbosity, write out message
	void logMessage(int priority, string message);

	/// Take a pointer and return a RefCountedPtr to the same data
	void  ptrToRefCountedPtr(Vector<LevelData<FArrayBox>* >& ptr,
			Vector<RefCountedPtr<LevelData<FArrayBox> > >& refPtr);

	/// Take a RefCountedPtr and return a pointer to the same data
	void refCountedPtrToPtr(Vector<RefCountedPtr<LevelData<FArrayBox> > >& refPtr,
			Vector<LevelData<FArrayBox>* >& ptr);

	/// Get the actual location (x,y,z) of an IntVect
	void getLocation(const IntVect iv, const int lev, RealVect &loc, Real ccOffsetX=0.5, Real ccOffsetY=0.5);

	/// Get a vector containing all var names (scalar and vector components)
	void getVarNames(Vector<string>& varNames);

	/// Get the coefficients for the VCAMRPoissonOp2 objects which we use to calculate the pressure. \n
	/// aCoef = \f$ 0\f$, \n
	/// bCoef = \f$ k(\chi)\f$ (permeability)
	void getPressureVCcoefs(Vector<RefCountedPtr<LevelData<FArrayBox> > >& aCoef,
			Vector<RefCountedPtr<LevelData<FluxBox> > >& bCoef);

	/// apply boundary conditions to a specific variable (if possible) and on a specific level (if specified)
	void applyBCs(const int a_var, const int lev);

	/// apply boundary conditions to vector field
	void applyVectorBCs(const int a_var, const int lev);

	/// apply boundary conditions to a specific variable (if possible) on all levels
	void applyBCs(int a_var);

	/// Subtract one field from the other and find the maximum difference (for normal pointers)
	Real calculateMaxDiff(Vector<LevelData<FArrayBox> * > a_phiNew, Vector<LevelData<FArrayBox> * > a_phiOld);

	/// Subtract one field from the other and find the maximum difference (for RefCounted pointers)
	Real calculateMaxDiff(Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phiNew, Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phiOld);

	/// Calculate the difference between two fields at a specific z-index,
	/// assuming there is no variation in the x-direction.
	Real calculateDiffAtPoint(Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phiNew,
			Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phiOld,
			int z_i);

	/// Gives a Vector<LevelData<FArrayBox>* > on the current AMR structure,
	/// with every value initialised to zero
	Vector<LevelData<FArrayBox>* >	newEmptyPtr();

	/// Gives a Vector<RefCountedPtr<LevelData<FArrayBox> > > on the current AMR structure,
	/// with every value initialised to zero
	Vector<RefCountedPtr<LevelData<FArrayBox> > >	newEmptyRefPtr();

	/// Define projection operators on current levels
	void setupProjectionOperators();


	/// max number of levels
	int m_max_level;

	/// current finest level
	int m_finest_level;

	/// blocking factor
	int m_block_factor;

	/// grid efficiency
	Real m_fill_ratio;

	/// proper nesting radius
	int m_nesting_radius;

	/// max box size
	int m_max_box_size;

	/// Regrid interval.
	/*
	 * Set to 1 to regrid at every timestep. Less than 1 means no regridding at all.
	 * */
	int m_regrid_interval;

	/// tagging value (undivided gradient threshold for regridding)
	Real m_tagging_val;

	/// Scalar field to use for refinement. If > -1 use this.
	int m_tagging_scalar_var;

	/// Vector field to use for refinement. Only use if m_tagging_scalar_var = -1.
	int m_tagging_vector_var;

	/// amount to buffer tags used in regridding
	int m_tags_grow;


	/// refinement ratio between this level and the next finer
	Vector<int>  m_refinement_ratios;

	/// Get refinement ratio between coarse level and next finer
	int getRefinementRatio(int coarseLevel);

	/// cell spacing at each level
	Vector<Real> m_amrDx;

	/// problem domains at each level
	Vector<ProblemDomain> m_amrDomains;

	/// problem domain size
	RealVect m_domainSize;

	/// How many grid cells we have on the coarsest level
	Vector<int> m_ncells;

	/// current grids
	Vector<DisjointBoxLayout> m_amrGrids;

	/// keeps track of which levels are completely covered
	Vector<int> m_covered_level;

	/// book-keeping; keeps strack of number of cells per level
	Vector<int> m_num_cells;

	/// current time
	Real m_time;

	/// most recent timestep
	Real m_dt;

	/// Cap dt, just incase something odd is happening with fluid velocities
	Real m_max_dt;

	/// In case we want to enforce a minimum number of timesteps
	int m_minTimestep;

	/// timestep scaling
	Real m_cfl;

	/// cfl number for initial timestep (useful if initial data needs small cfl)
	Real m_initial_cfl;

	/// maximum amount cfl number may grow in one timestep
	Real m_max_dt_grow;

	/// current step
	int m_cur_step;

	/// have we reached steady state?
	bool m_steady_state;

	/// Parameters used in the model
	MushyLayerParams m_parameters;

	/// Try to fix velocity
	bool m_debugging_fixVelocity;

	/// Resets temperature field
	bool m_debugging_resetTemperatureField;



	/**************
	 * New stuff
	 */
	//Vector<Vector<LevelData<FArrayBox>*> > m_scalarNew, m_scalarOld, m_dScalar;
	/// Scalar fields at new timestep
	Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > m_scalarNew;

	/// Scalar fields at old timestep
	Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > m_scalarOld;

	/// Difference between some scalar fields, source term etc - used in different ways
	/// to pass data between functions
	Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > m_dScalar;

	/// Vector fields at new timestep
	Vector<Vector<LevelData<FArrayBox>*> > m_vectorNew;

	/// Vector fields at old timestep
	Vector<Vector<LevelData<FArrayBox>*> > m_vectorOld;

	/// Difference between vectors field, or a vector field source term - used in
	/// different ways at different points in the code
	Vector<Vector<LevelData<FArrayBox>*> > m_dVector;

	/// List of vector fields
	enum vectorVars {
		m_fluidVel,
		m_fluidVelPred,
		m_gradPressure,
		m_fluidVelAnalytic,
		m_fluidVelErr,
		m_gradPressureErr
	};

	/// Number of vector fields. Must be set manually.
	int m_numVectorVars;

	/// List of scalar fields
	enum vars {
		m_enthalpy,
		m_enthalpySolid,
		m_enthalpyEutectic,
		m_enthalpyLiquid,
		m_composition,
		m_theta,
		m_thetaLiquidus,
		m_thetaSolidus,
		m_porosity,
		m_porosityEutectic,
		m_compositionLiquid,
		m_compositionSolid,
		m_thetaForcing,
		m_liquidCompositionGrad,
		m_solidFraction,
		m_steadyStateImbalance,
		m_Tanalytic,
		m_ThetaLAnalytic,
		m_solidFractionTrue,
		m_thetaTrue,
		m_ThetaAnalytic,
		m_enthalpyAnalytic,
		m_enthalpyAdvection,
		m_thetaLaplacian,
		m_streamFunction,
		m_resid,
		m_concSource,
		m_ThetaDiffusion,
		m_ThetaDiffusionN,
		m_thetaBcoef,
		m_ThetaBCoef,
		m_permeability,
		m_pressure,
		m_divU,
		m_ThetaFrameAdvection,
		m_ThetaSSource,
		m_ThetaPorositySource,
		m_pressureError,
		m_pressureAnalytic,
		m_divUerror,
		m_thetaForcingError,
		m_thetaForcingAnalytic,
		m_CFInterpTest,
		m_divUstar,
		m_divUstarErr,
		m_lambda,
		m_lambdaPostCorr,
		m_buoyancyFilter
	};


	/// Number of scalar fields. Must be set manually.
	int m_numVars;

	/// Frame advection velocity, constant
	Vector<LevelData<FluxBox>* > m_frameAdv;

	/// Darcy velocity (mass flux)
	Vector<LevelData<FluxBox>* > m_fluidAdv;

	/// To map a scalar variable to it's name
	Vector<string> m_varNames;

	/// To map a vector variable to it's name
	Vector<string> m_vectorVarNames;

	/// Which problem to solve?
	/**
	 0 - full problem \n
	 1 - benchmark 1, directional solidification without flow \n
	 2 - benchmark 2, \n
	 3 - heat diffusion in a liquid
	 */
	enum problemType {
		m_fullProblem,
		m_bmDiffusiveSolidification,
		m_bmHRL,
		m_bmDiffusion,
		m_bmAdvection,
		m_enforcedPorosity
	};

	/// How should we calculate \f$ V d\phi / dt \f$ terms?
	int m_frameAdvectionMethod;

	/// Different ways of calculating the frame advection term
	enum frameAdvectionMethods{
		m_godunov,
		m_finiteDifference
	};

	/// 2 or 4. How accurately to calculate divergences and gradients
	/// when finding U
	int m_fluidVelDiffOrder;

	/// 0,1,2. How accurately to interpolate from cell centred to face centred values
	/// when finding U
	int m_fluidVelInterpOrder;

	/// Whether to calculate velocity using my hand written code or the code copied from CCProjector
	bool m_ownVelCalc;

	/// calculate the permeability \f$ k=k(\chi) \f$, using the method prescribed by m_permeabilityFunction
	void calculatePermeability(bool oldTime = false);

	/// Different phases we could be in
	enum domainRegions {
		m_solidRegion,
		m_eutecticRegion,
		m_mushyRegion,
		m_liquidRegion
	};

	/**************************************
	 * Boundary conditions
	 **************************************/


	/**************************************
	 * Solver variables
	 *************************************/

	/// Should we do 1st or 2nd order integration in time? \n
	/// 1-backward euler
	/// 2-TGA
	int m_timeIntegrationOrder;

	/// Advection operator for temperature
	Vector<RefCountedPtr<AdvectPhysics> > m_advPhystheta,

	/// Advection operator for liquid concentration
	m_advPhysConcLiquid,

	/// Advection operator for enthalpy
	m_advPhysEnthalpy,

	/// Advection operator for bulk concentration
	m_advPhysConc,

	/// Advection operator for lambda
	m_advPhysLambda;

	/// LevelAdvect operator for temperature
	Vector<RefCountedPtr<LevelAdvect> > m_godunovtheta,

	/// LevelAdvect operator for liquid concentration
	m_godunovConcLiquid,

	/// LevelAdvect operator for bulk concentration
	m_godunovConc,

	/// LevelAdvect operator for enthalpy
	m_godunovEnthalpy;


	/// Flux registers for scalar vars
	/**
	 * Flux register lives on the coarser level
	 */
	Vector<Vector<RefCountedPtr<LevelFluxRegister> > > m_fluxRegister;

	/// Flux registers for vector vars
	/**
	 * Flux register lives on the coarser level
	 */
	Vector<Vector<RefCountedPtr<LevelFluxRegister> > > m_vectorFluxRegister;

	/// Flux registers for velocity
//	Vector<RefCountedPtr<LevelFluxRegisterEdge> > m_velFluxReg;
	//	LevelFluxRegister m_fluxReg;

	/// Operator factory for \f$ \theta \f$ (using a simple poisson solver)
	RefCountedPtr<AMRPoissonOpFactory> thetaPoissonOpFact;

	/// Operator factory for theta, using a variable coefficient poisson solver
	RefCountedPtr<VCAMRPoissonOp2Factory> thetaVCOpFact,

	/// Operator factory for Theta_l, using a variable coefficient poisson solver
	ThetaLVCOpFact;

	/// Op factory for solving for reflux corrections
	RefCountedPtr<AMRPoissonOpFactory> m_refluxOpFact;

	/// Diffusion/multigrid bottom solver
	BiCGStabSolver<LevelData<FArrayBox> > s_bottomSolver,

	/// Bottom solver for variable coefficient computations
	s_bottomSolverVC;

	/// Backward Euler timestepping (temperature)
	RefCountedPtr<BackwardEuler> m_diffuseBEtheta,

	/// Backward Euler timestepping (liquid salinity)
	m_diffuseBEThetaL;

	///  TGA timestepping for temperature
	RefCountedPtr<AMRTGA<LevelData<FArrayBox> > > m_diffuseTGAtheta,

	///  TGA timestepping for liquid salinity
	m_diffuseTGAThetaL;

	///Multigrid solver for \f$ \Theta_l \f$
	RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > > m_amrSolverVCThetaL;

	///Multigrid solver for \f$ \theta \f$
	RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > > m_amrSolverVCtheta;

	///Multigrid solver for reflux corrections
	RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > > m_refluxAMRMG;

	/// Pointer to physics boundary condition object
	PhysBCUtil* m_physBCPtr;

	/**********************************/


	/// how verbose should we be?
	static int s_verbosity;

	/// is this object initialized?
	bool m_is_defined;

	///
//	bool m_do_restart;

	/// if starting from a restart, timestep of restart
	int m_restart_step;

	/// initially refine entire domain (rather than generating tags)
	/** If true, the entire domain is refined to the finest possible
      resolution at the initial step, rather than generating tags, etc.
      This can be useful for random initial conditions 
	 */
	bool m_refine_initial_domain;

	/// Folder (relative to current directory) into which to write output file
	string m_output_folder;

	/// Prefix for plot files
	string m_plot_prefix;

	/// Prefix for checkpoint files
	string m_check_prefix;

	/// Determines after how many timesteps we write a plot file (1 = every timestep, -1 = never)
	int m_plot_interval;

	/// Determines after how many iterations (within a timestep) we write a plot file (1 = every timestep, -1 = never)
	int m_iteration_plot_interval;

	/// Determines after how many timesteps we write a checkpoint file (1 = every timestep, -1 = never)
	int m_check_interval;

	/// Whether or not the model should just output the analytic solution and stop
	bool m_printAnalyticSoln;

	/// Number of ghost cells
	int m_num_ghost;

	/// Is domain periodic in each direction?
	bool m_is_periodic[SpaceDim];

	/// Ghostvector of length m_num_ghost
	IntVect m_ghostVect;

	/// Whether or not to initialise from the analytic solution (if possible)
	bool m_enforceAnalyticSoln;

	/// Keep a count of how many iterations we've done in total
	int m_total_num_iterations;


};  



#endif


