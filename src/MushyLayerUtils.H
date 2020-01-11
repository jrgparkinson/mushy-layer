
#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MUSHYLAYERUTILS_H_
#define _MUSHYLAYERUTILS_H_

#include <cmath>

#include "AMR.H"
#include "memtrack.H"
#include "parstream.H"
#include "ParmParse.H"
#include "MushyLayerParams.h"

#include "NamespaceHeader.H"



/// Compute solute flux
Real computeSoluteFlux(const LevelData<FArrayBox>& S, const LevelData<FArrayBox>& Sl, const LevelData<FArrayBox>& porosity,
                       const LevelData<FluxBox>& vel, Real m_dx, MushyLayerParams m_parameters,
                       int a_dir = 0,
                       Side::LoHiSide side = Side::Lo);

/// Compute nusselt number
Real computeNusselt(const LevelData<FArrayBox>& T, const LevelData<FArrayBox>& vel,
                    Real a_dx, MushyLayerParams a_params, Real a_domWidth, Real a_domHeight);

/// Compute the largest value and retain it's sign
/**
 * Sometimes we want to know the sign of the largest value, not just it's magnitude.
  *e.g. to tell is residuals are oscillating or decreasing monotonically
*/
Real computeNormWithSign(const Vector<LevelData<FArrayBox>* >& a_phi,
		const Vector<int>&                    a_nRefFine,
		const Interval                        a_comps = Interval(0,0),
		const int                             a_lBase = 0);



//extern Real
//dtgNorm(const Vector< LevelData<FArrayBox>* >& a_src,
//		const Vector< DisjointBoxLayout >&     a_grids,
//		const Vector<int>&                     a_refRatio,
//		const ProblemDomain&                   a_coarseDom,
//		const int& a_comp,
//		const int& a_pval);

/// Define the ProblemDomain object from the inputs files
/**
 */
extern void
getProblemDomain(ProblemDomain& a_domain);

/// Setup AMR for new run
/**
 * Works out if we're restarting, using fixed grids or starting a new AMR run
 */
extern void
setupAMRForAMRRun(AMR& a_amr, ProblemDomain prob_domain);

/// Get variable names for I/O
void getVarNames(Vector<string>& varNames,
		const int numScalarVars, const int numVectorVars,
		const Vector<string> scalarVarNames, const Vector<string> vectorVarNames);

/// Get variable names for I/O
void getScalarVarNames(Vector<string>& varNames,
                 const int numScalarVars,
                 const Vector<string> scalarVarNames);

/// Get variable names for I/O
void getChkVectorVarNames(Vector<string>& chkVectorVarNames,
                          const Vector<int> chkVectorVars,
                          const Vector<string> vectorVarNames
                          );

/// Get variable names for I/O
void getVectorVarNames(Vector<string>& varNames,
                       const int numVectorVars,
                       const Vector<string> vectorVarNames,
                       const int startPos = 0);

/// Get variable names for I/O
void getVarNames(Vector<string>& varNames,
		const Vector<int> scalarVars, const Vector<int> vectorVars,
		const Vector<string> scalarVarNames, const Vector<string> vectorVarNames);

/// Print out which revision of the repository we're using
void printRepoVersion();

/// Calculate the permeability
void
calculatePermeability(FArrayBox& permeabilityFAB, FArrayBox& solidFractionFAB,
		MushyLayerParams& params, Real a_dx);

/// Calculate the permeability - moved to MushyLayerParams
//Real calculatePermeability(Real liquidFraction,
//		MushyLayerParams& params);

/// Get the physical location in the domain
/**
 * default offset is for cell centred points.
 * note that ccOffset should be 0, 0.5 or 1 in each direction
 */
void getLocation(const IntVect iv, RealVect& loc, const Real a_dx,
		const RealVect ccOffset = 0.5*RealVect::Unit);

/// Get location in domain
/**
 * default offset is for cell centred points.
 * note that ccOffset should be 0, 0.5 or 1 in each direction
 */
void getLocation(const IntVect iv, const Real a_dx, Real& x, Real& y, Real& z,
                 const RealVect ccOffset = 0.5*RealVect::Unit);

/// Initial data for burgers problem benchmark
Real burgersPeriodicInit(Real x, Real y, int dir, MushyLayerParams params);

/// Initial data for burgers problem benchmark
Real burgersSinInit(Real x, Real y, int dir, const MushyLayerParams& params);

/// Initial data for Stokes-Darcy channel flow benchmark
Real stokesDarcyInit(Real x, Real y, int dir, const MushyLayerParams& params);

/// Read in fixed grids
extern bool
getFixedGrids(Vector<Vector<Box> >& amrGrids,  ProblemDomain prob_domain, string gridfileParam = "gridfile");

/// Check if string is also an integer
bool is_integer(const std::string& s);

/// Create data structures for holding phase boundaries
void createPhaseBoundaryStructures(const int lev, const Vector<DisjointBoxLayout>& grids, const IntVect ivGhost,
                                   Vector<RefCountedPtr<LevelData<FArrayBox> > >& enthalpySolidus,
                                   Vector<RefCountedPtr<LevelData<FArrayBox> > >& enthalpyEutectic,
                                   Vector<RefCountedPtr<LevelData<FArrayBox> > >& enthalpyLiquidus);


#include "NamespaceFooter.H"

#endif
