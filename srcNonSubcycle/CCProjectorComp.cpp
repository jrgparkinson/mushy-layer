#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CCProjectorComp.H"
#include "Gradient.H"
#include "GradientF_F.H"
#include "Divergence.H"
#include "AMRPoissonOp.H"
#include "RelaxSolver.H"
#include "CoarseAverageEdge.H"
#include "CornerCopier.H"
#include "AMRIO.H"
#include "computeSum.H"
#include "SetValLevel.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "BiCGStabSolver.H"
#include "VelBCHolder.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"
#include "PiecewiseLinearFillPatchFace.H"
#include "VCAMRPoissonOp2.H"
//#include "NodeQuadCFInterp2.H"

#ifdef CH_USE_HDF5
#include "CH_HDF5.H"
#endif

/// define static variables here
bool CCProjectorComp::m_doSyncProjection = true;
bool CCProjectorComp::m_applySyncCorrection = true;
bool CCProjectorComp::s_applyVDCorrection = true;
bool CCProjectorComp::m_doQuadInterp = true;
Real CCProjectorComp::m_etaLambda = 0.9;

#if defined(CH_USE_DOUBLE)
Real CCProjectorComp::s_solver_tolerance = 1.0e-10; // pmc, 18 oct 2007:  was -10
#elif defined(CH_USE_FLOAT)
Real CCProjectorComp::s_solver_tolerance = 1.0e-5; // pmc, 18 oct 2007:  was -10
#endif

int  CCProjectorComp::s_num_smooth_up = 4;
int  CCProjectorComp::s_num_smooth_down = 4;
int  CCProjectorComp::s_num_precond_smooth = 4;
bool CCProjectorComp::s_constantLambdaScaling = false;
Real CCProjectorComp::s_lambda_timestep = 0.0;
bool CCProjectorComp::pp_init = false;
int  CCProjectorComp::s_verbosity = 1;

/// first define quick-n-easy accessfunctions

// ------------------------------------------------------
int CCProjectorComp::getLevel() const
{
	return m_level;
}



// ------------------------------------------------------
Real CCProjectorComp::etaLambda() const
{
	return m_etaLambda;
}

// ------------------------------------------------------

// ------------------------------------------------------
//LevelData<FArrayBox>& CCProjectorComp::eSync()
//{
//	return m_eSync;
//}
//
//// ------------------------------------------------------
//const LevelData<FArrayBox>& CCProjectorComp::eSync() const
//{
//	return m_eSync;
//}

// ------------------------------------------------------
void CCProjectorComp::getPressure(Vector<RefCountedPtr<LevelData<FArrayBox> > >& pi)
{
	for (int lev = 0; lev <= finest_level; lev++)
	{
		m_pressure[lev]->copyTo(*pi[lev]);
	}
}

void CCProjectorComp::getGradPressure(Vector<RefCountedPtr<LevelData<FArrayBox> > >& gradP)
{
	for (int lev = 0; lev <= finest_level; lev++)
	{
		m_gradPressure[lev]->copyTo(*gradP[lev]);
	}
}

void CCProjectorComp::getGradPressure(Vector<LevelData<FArrayBox> * >& gradP)
{
	for (int lev = 0; lev <= finest_level; lev++)
	{
		m_gradPressure[lev]->copyTo(*gradP[lev]);
	}
}

void CCProjectorComp::getGradPressureEdge(Vector<LevelData<FluxBox> * >& gradPedge)
{
	for (int lev = 0; lev <= finest_level; lev++)
	{
		m_gradPressureEdge[lev]->copyTo(*gradPedge[lev]);
	}
}

// ------------------------------------------------------
//LevelData<FArrayBox>& CCProjectorComp::eLambda()
//{
//	return m_eLambda;
//}
//
//// ------------------------------------------------------
//const LevelData<FArrayBox>& CCProjectorComp::eLambda() const
//{
//	return m_eLambda;
//}

// ------------------------------------------------------
//LevelData<FluxBox>& CCProjectorComp::grad_eLambda()
//{
//	return m_grad_eLambda;
//}
//
//// ------------------------------------------------------
//const LevelData<FluxBox>& CCProjectorComp::grad_eLambda() const
//{
//	return m_grad_eLambda;
//}

// ------------------------------------------------------
bool CCProjectorComp::doSyncProjection() const
{
	return m_doSyncProjection;
}

// ------------------------------------------------------
bool CCProjectorComp::doMacSync() const
{
	return (m_etaLambda > 0.0);
}

// ------------------------------------------------------
bool CCProjectorComp::isInitialized() const
{
	return m_isInitialized;
}

// ------------------------------------------------------
bool CCProjectorComp::doQuadInterp() const
{
	return m_doQuadInterp;
}

// ------------------------------------------------------
QuadCFInterp& CCProjectorComp::quadCFInterpolator()
{
	return m_cfInterp;
}

/// returns predefined intvectset which shows coverage
const LayoutData<IntVectSet>& CCProjectorComp::getGridIVS()
{
	return m_gradIVS;
}

// ------------------------------------------------------
bool CCProjectorComp::isFinestLevel() const
{
	return m_finest_level;
}

// ------------------------------------------------------
void CCProjectorComp::isFinestLevel(bool a_finest_level)
{
	m_finest_level = a_finest_level;
}

// ------------------------------------------------------
void CCProjectorComp::verbosity(int a_verbosity)
{
	s_verbosity = a_verbosity;
}

// ------------------------------------------------------
int CCProjectorComp::verbosity() const
{
	return s_verbosity;
}

// ------------------------------------------------------
void CCProjectorComp::setPhysBC(PhysBCUtil& a_bc)
{
	//Don't do this - we don't own it!
	//  if (m_physBCPtr != nullptr)
	//    {
	//      delete m_physBCPtr;
	//      m_physBCPtr = nullptr;
	//    }
	//  m_physBCPtr = a_bc.newPhysBCUtil();
	m_physBCPtr = (&a_bc);
}

// ------------------------------------------------------
PhysBCUtil* CCProjectorComp::getPhysBCPtr() const
{
	return m_physBCPtr;
}

// ------------------------------------------------------
// default constructor
CCProjectorComp::CCProjectorComp()
{
	finest_level = -1;
	m_refinement_ratios = Vector<int>(1,1);
	m_isInitialized = false;
	// have to initialize this to _something_!
	m_finest_level = false;
	finest_level = -1;
	m_physBCPtr = nullptr;
	m_bottomSolverLevel = nullptr;
	m_limitSolverCoarsening = false;
	m_amrDx = Vector<Real>(1,1);
}

// ------------------------------------------------------

void CCProjectorComp::undefine()
{
	for (int lev=0; lev <=finest_level; lev++)
	{
		if(m_pressure[lev] != nullptr)
		{
			delete m_pressure[lev];
			//				m_pressure[lev] = nullptr;
		}

		if(m_gradPressure[lev] != nullptr)
		{
			delete m_gradPressure[lev];
			//				m_gradPressure[lev] = nullptr;
		}
		if(m_phi[lev] != nullptr)
		{
			delete m_phi[lev];
			//				m_phi[lev] = nullptr;
		}
		if(m_divUstar[lev] != nullptr)
		{
			delete m_divUstar[lev];
			//				m_divUstar[lev] = nullptr;
		}

		if(m_amrGradIVS[lev] != nullptr)
		{
			delete m_amrGradIVS[lev];
			//				m_amrGradIVS[lev] = nullptr;
		}
		if(m_eLambda[lev] != nullptr)
		{
			delete m_eLambda[lev];
		}
		if(m_grad_eLambda[lev] != nullptr)
		{
			delete m_grad_eLambda[lev];
		}
		if(m_eSync[lev] != nullptr)
		{
			delete m_eSync[lev];
		}
		if(m_gradPressureEdge[lev] != nullptr)
		{
			delete m_gradPressureEdge[lev];
		}

	}

	if (m_bottomSolver!=nullptr)
	{
		delete m_bottomSolver;
	}

	m_isInitialized = false;

}
// destructor
CCProjectorComp::~CCProjectorComp()
{
	//We don't own this! It belongs to amrMushyLayer
	//  if (m_physBCPtr != nullptr)
	//    {
	//      delete m_physBCPtr;
	//      m_physBCPtr = nullptr;
	//    }

	//	if(m_bottomSolverLevel != nullptr)
	//	{
	//		delete m_bottomSolverLevel;
	//		m_bottomSolverLevel = nullptr;
	//	}


	// everything else should be automatic here
}

// ------------------------------------------------------
void
CCProjectorComp::define(const Vector<DisjointBoxLayout>& a_amrGrids,
		const Vector<ProblemDomain>& a_amrDomains,
		const Vector<Real>&  a_amrDx,
		const Vector<int>& a_refinementRatios,
		const int a_finest_level,
		PhysBCUtil& a_physBC,
		Vector<RefCountedPtr<LevelData<FArrayBox> > > a_theta,
		Vector<RefCountedPtr<LevelData<FArrayBox> > > a_permeability,
		int a_numGhost)
{
	// set physical boundary condition object
	setPhysBC(a_physBC);

	if (!pp_init)
	{
		variableSetUp();
	}

	CH_assert (a_finest_level > -1);
	m_amrDx = a_amrDx;
	m_amrDomains = a_amrDomains;

	finest_level = a_finest_level;
	m_amrGrids = a_amrGrids;
	m_refinement_ratios = a_refinementRatios;


	num_ghost = a_numGhost;
	IntVect ghostVect(num_ghost*IntVect::Unit);
	IntVect smallGhostVect((num_ghost-1)*IntVect::Unit); // Grad Pressure should have less ghost cells than pressure
	int numLevels = a_finest_level + 1;

	m_amrCFInterps.resize(numLevels, nullptr);
	m_amrFillPatchFace.resize(numLevels, nullptr);

	m_permeability.resize(numLevels);

	m_amrGradIVS.resize(numLevels, nullptr);
	m_pressure.resize(numLevels, nullptr);
	m_gradPressure.resize(numLevels, nullptr);
	m_gradPressureEdge.resize(numLevels, nullptr);
	m_grad_eLambda.resize(numLevels, nullptr);
	m_eLambda.resize(numLevels, nullptr);
	m_eSync.resize(numLevels, nullptr);
	m_phi.resize(numLevels, nullptr);
	m_divUstar.resize(numLevels, nullptr);

	m_theta = a_theta;

	for (int lev = 0; lev <= finest_level; lev++)
	{
		if (lev > 0)
		{
			m_amrCFInterps[lev] = new QuadCFInterp(m_amrGrids[lev], &m_amrGrids[lev-1], m_amrDx[lev], m_refinement_ratios[lev-1],
					1,m_amrDomains[lev]);


			m_amrFillPatchFace[lev] = new PiecewiseLinearFillPatchFace(m_amrGrids[lev], m_amrGrids[lev-1], 1, m_amrDomains[0], m_refinement_ratios[lev-1], 2);
		}

		// We generally store permeability as cell centered, but need it on the face
		// for the VCAMRPoissonOp
		m_permeability[lev] = RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>(m_amrGrids[lev], 1, ghostVect));
		CellToEdge(*a_permeability[lev], *m_permeability[lev]);

//		for (DataIterator dit = m_permeability[lev]->dataIterator(); dit.ok(); ++dit)
//		{
//			FluxBox& perm = (*m_permeability[lev])[dit];
////			perm.setVal(1.0);
//			int temp=0;
//		}

		m_pressure[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, ghostVect);
		m_eLambda[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, ghostVect);
		m_eSync[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, ghostVect);
		m_gradPressure[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], SpaceDim, ghostVect);
		m_gradPressureEdge[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, ghostVect); //smallGhostVect*IntVect::Unit
		m_grad_eLambda[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, ghostVect);
		m_phi[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, ghostVect);
		m_divUstar[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);

		m_amrGradIVS[lev] = new LayoutData<IntVectSet>();
		//		Gradient::createGridIVS(*m_amrGradIVS[lev], m_amrGrids[lev], num_ghost);
		Gradient::createGridIVS(*m_amrGradIVS[lev], m_pressure[lev]->disjointBoxLayout(), 1); // trying to calculate ghost cells aswell

		// Initialize to a bogus value
		setValLevel(*m_pressure[lev], 0);
		setValLevel(*m_gradPressure[lev], 0);
		setValLevel(*m_phi[lev], 0);
		setValLevel(*m_divUstar[lev], 0);
		setValLevel(*m_eLambda[lev], 0);
		setValLevel(*m_grad_eLambda[lev], 0);
		setValLevel(*m_gradPressureEdge[lev], 0);
		setValLevel(*m_eSync[lev], 0);
	}



	// define persistent storage


	// these need to start as 0 for boundary conditions
	// for first intermediate solve if more than one level
	// of refinement
	//	setValLevel(m_eLambda, 0.0);
	//	setValLevel(m_eSync, 0.0);

	//	defineSolverMGlevel(a_grids, a_crseGridsPtr);
	m_bottomSolver = new BiCGStabSolver<LevelData<FArrayBox> >;
	m_bottomSolver->m_verbosity = s_verbosity-2;
	defineSolverMG();

	m_isInitialized = false;
}

// ------------------------------------------------------
CCProjectorComp::CCProjectorComp(const Vector<DisjointBoxLayout>& a_amrGrids,
		const Vector<ProblemDomain>& a_amrDomains,
		const Vector<Real>&  a_amrDx,
		const Vector<int>& a_refinementRatios,
		const int a_finest_level,
		PhysBCUtil& a_physBC,
		Vector<RefCountedPtr<LevelData<FArrayBox> > > a_theta,
		Vector<RefCountedPtr<LevelData<FArrayBox> > > a_permeability,
		int a_numGhost)
{
	define(a_amrGrids, a_amrDomains, a_amrDx, a_refinementRatios, a_finest_level, a_physBC,
			a_theta, a_permeability, a_numGhost);
}

// -------------------------------------------------------
void CCProjectorComp::init(const CCProjectorComp& a_oldProj)
{
	// this function performs a partial initialization
	// of this object using the oldProj...

	// at the current state of this algorithm,
	// this doesn't need to do anything (everything will
	// need to be re-initialized anyway)
}

void CCProjectorComp::variableSetUp()
{
	// first set up parm parse object
	ParmParse pp("projection");

	int tempBool;

	tempBool = (int) m_doSyncProjection;
	pp.query("doSyncProjection",tempBool);
	m_doSyncProjection = (tempBool == 1);

	tempBool = (int) m_applySyncCorrection;
	pp.query("applySyncCorrection", tempBool);
	m_applySyncCorrection = (tempBool == 1);

	tempBool = (int) s_applyVDCorrection;
	pp.query("applyFreestreamCorrection", tempBool);
	s_applyVDCorrection = (tempBool == 1);

	tempBool = (int) m_doQuadInterp;
	pp.query("doQuadInterp", tempBool);
	m_doQuadInterp = (tempBool == 1);

	pp.query("eta", m_etaLambda);

	pp.query("solverTol", s_solver_tolerance);
	pp.query("numSmoothUp", s_num_smooth_up);
	pp.query("numSmoothDown", s_num_smooth_down);
	pp.query("numPrecondSmooth", s_num_precond_smooth);

	tempBool = (int) s_constantLambdaScaling;
	pp.query("constantLambdaScaling", tempBool);
	s_constantLambdaScaling = (tempBool == 1);

	// now dump this all back out
	if (s_verbosity > 2)
	{
		pout () << "Projection inputs: " << endl;

		pout() << "  doSyncProjection = "  << m_doSyncProjection << endl;

		pout() << "  applySyncCorrection = "  << m_applySyncCorrection << endl;

		pout() << "  applyFreestreamCorrection = "  << s_applyVDCorrection
				<< endl;

		pout() << "  doQuadInterp = "  << m_doQuadInterp << endl;

		pout() << "  eta = "  <<  m_etaLambda << endl;

		pout() << "  constantLambdaScaling = "  << s_constantLambdaScaling
				<< endl;
	}

	// set flag so we don't have to call this function again.
	pp_init = true;

}



void CCProjectorComp::limitSolverCoarsening(bool a_limitSolverCoarsening)
{
	m_limitSolverCoarsening = a_limitSolverCoarsening;
}

// --------------------------------------------------------------
// note that this function assumes that all BC's on phi have been set
void CCProjectorComp::gradPhi(LevelData<FArrayBox>& a_gradPhi, int a_dir) const
{
	MayDay::Error("gradPhi - not implemented");
	// first compute edge-centering of gradPhi
	//	DataIterator dit = a_gradPhi.dataIterator();
	//	dit.reset();
	//
	//	int edgeDir = -1;
	//	const Box& edgeBox = a_gradPhi[dit].box();
	//	for (int dir=0; dir<SpaceDim; dir++)
	//	{
	//		if (edgeBox.type(dir) == IndexType::NODE)
	//		{
	//			CH_assert (edgeDir == -1);
	//			edgeDir = dir;
	//		}
	//	}
	//	CH_assert(edgeDir != -1);
	//
	//	const DisjointBoxLayout& grids = getBoxes();
	//
	//	for (dit.reset(); dit.ok(); ++dit)
	//	{
	//		// call gradient subroutine directly, since this is one-directional
	//
	//		// I don't think this function is called by anybody ever
	//		// I think this is correct, but just in case, flag it
	//		// (dfmartin@lbl.gov)
	//		CH_assert (false);
	//
	//		Gradient::singleBoxMacGrad(a_gradPhi[dit],
	//				m_phi[dit],
	//				0, 0, 1, grids[dit], m_dx,
	//				a_dir, edgeDir, m_gradIVS[dit]);
	//	} // end loop over grids
}

// --------------------------------------------------------------
// note that this function assumes that all BC's on phi have been set
//void CCProjectorComp::gradPhi(LevelData<FluxBox>& a_gradPhi) const
//{
//	Gradient::levelGradientMAC(a_gradPhi, m_phi, m_dx);
//}

///
void
CCProjectorComp::gradPhiComp(Vector<RefCountedPtr<LevelData<FluxBox> > >& a_gradPhi) const
{
	//	Gradient::levelGradientMAC(a_gradPhi, m_phi, m_dx);
	MayDay::Error("gradPhiComp - not implemented");
	//	LevelData<FArrayBox>* a_phiCrsePtr = nullptr;
	//	LevelData<FArrayBox>* a_phiFinePtr = nullptr;
	//
	//	for (int lev=0; lev <= finest_level; lev++)
	//	{
	//		if (lev > 0)
	//		{
	//			a_phiCrsePtr = &(*m_phi[lev-1]);
	//
	//		}
	//		if (lev < finest_level)
	//		{
	//			a_phiFinePtr = &(*m_phi[lev+1]);
	//		}
	//
	//		Gradient::compGradientMAC(*a_gradPhi[lev], *m_phi[lev], a_phiCrsePtr, a_phiFinePtr ,m_amrDx[lev],
	//				m_refinementRatio[lev-1], m_refinementRatio[lev], m_amrDomains[lev]);
	//	}
}


// --------------------------------------------------------------
// note that this function assumes that all BC's on Pi have been set
void CCProjectorComp::gradPi(LevelData<FArrayBox>& a_gradPi, int a_dir) const
{
	MayDay::Error("gradPi - not implemented");
	//	Gradient::compGradientCC();
	//
	//	for (int lev=0; lev <= finest_level; lev++)
	//	{
	//		DataIterator dit = a_gradPi.dataIterator();
	//
	//		for (dit.reset(); dit.ok(); ++dit)
	//		{
	//			// call FORTRAN subroutine directly...
	//			FORT_GRADCC(CHF_FRA1(a_gradPi[dit],0),
	//					CHF_CONST_FRA1(m_Pi[dit],0),
	//					CHF_BOX(getBoxes()[dit]),
	//					CHF_CONST_REAL(m_dx),
	//					CHF_INT(a_dir));
	//		}
	//	}
}

// --------------------------------------------------------------
// note that this function assumes that all BC's on Pi have been set
void CCProjectorComp::gradPi(LevelData<FArrayBox>& a_gradPi) const
{
	//todo - implement this
}

// -------------------------------------------------------------
// note that this function assumes all BC's on eSync have been set
// also only returns valid grad(eSync) in valid regions of grids
void CCProjectorComp::grad_eSync(LevelData<FArrayBox>& a_grad_eSync) const
{
	MayDay::Error("grad_eSync not implemented");
	//	if (m_fineProjPtr != nullptr)
	//	{
	//		const LevelData<FArrayBox>& fine_eSync = m_fineProjPtr->eSync();
	//		int nRefFine = m_fineProjPtr->nRefCrse();
	//
	//		Gradient::compGradientCC(a_grad_eSync, m_eSync,  &fine_eSync,
	//				m_dx, nRefFine, m_domain);
	//	}
	//	else
	//	{
	//		// no finer level -- can do level-operator gradient
	//		Gradient::levelGradientCC(a_grad_eSync, m_eSync, m_dx);
	//	}
}

// --------------------------------------------------------------
// note that this function assumes all BC's on eSync have been set
// also only returns valid grad(eSync) in valid regions of grids
void CCProjectorComp::grad_eSync(LevelData<FArrayBox>& a_grad_eSync, int a_dir) const
{
	// this is not yet implemented (not sure if i'll really need it!)
	MayDay::Warning("directional grad_eSync not yet implemented!");
}

// --------------------------------------------------------------
void CCProjectorComp::levelMacProject(LevelData<FluxBox>& a_uEdge,
		Real a_oldTime, Real a_dt)
{
	//todo - write this

	//	// MAC rhs should have no ghost values
	//	LevelData<FArrayBox> MacRHS(getBoxes(),1);
	//
	//	Divergence::levelDivergenceMAC(MacRHS, a_uEdge, m_dx);
	//
	//	LevelData<FArrayBox>* crseBCPtr = nullptr;
	//	Real CFscale = 0.5*a_dt;
	//
	//	if (m_crseProjPtr != nullptr)
	//	{
	//		// coarse-fine BC is 0.5*dt*(coarse Pi)
	//		const DisjointBoxLayout crseGrids = m_crseProjPtr->getBoxes();
	//		crseBCPtr = new LevelData<FArrayBox>(crseGrids,1);
	//		LevelData<FArrayBox>& crseBC = *crseBCPtr;
	//
	//		Interval comps(0,0);
	//		const LevelData<FArrayBox>& crsePi = m_crseProjPtr->Pi();
	//		crsePi.copyTo(comps,crseBC,comps);
	//
	//		DataIterator dit = crseBC.dataIterator();
	//		for (dit.reset(); dit.ok(); ++dit)
	//		{
	//			crseBC[dit] *= CFscale;
	//		}
	//	}
	//
	//	// report sum(rhs)
	//	if (s_verbosity >= 2)
	//	{
	//		DisjointBoxLayout* finerGridsPtr = nullptr;
	//		int nRefFine = -1;
	//		Real sumRHS = computeSum(MacRHS, finerGridsPtr,
	//				nRefFine, m_dx, MacRHS.interval());
	//		pout() << "MAC projection -- sum(RHS) = " << sumRHS << endl;
	//	}
	//
	//	// now solve for phi
	//	solveMGlevel(m_phi, crseBCPtr, MacRHS);
	//
	//	// now correct
	//
	//	// first re-apply physical and copy BC's
	//	Interval phiComps(0,0);
	//	m_phi.exchange(phiComps);
	//	BCHolder bcHolder = m_physBCPtr->gradMacPressureFuncBC();
	//	const DisjointBoxLayout& levelGrids = getBoxes();
	//
	//	DataIterator dit = m_phi.dataIterator();
	//	for (dit.reset(); dit.ok(); ++dit)
	//	{
	//		bcHolder.operator()(m_phi[dit],
	//				levelGrids[dit],
	//				m_domain,
	//				m_dx,
	//				false); // not homogeneous
	//	}
	//
	//#define CORNERCOPIER
	//#ifdef CORNERCOPIER
	//CornerCopier cornerExchangeCopier(m_phi.getBoxes(),
	//		m_phi.getBoxes(),
	//		m_domain,
	//		m_phi.ghostVect(),
	//		true);
	//
	//m_phi.exchange(m_phi.interval(), cornerExchangeCopier);
	//#endif
	//
	//applyMacCorrection(a_uEdge, CFscale);
	//
	//if (crseBCPtr != nullptr)
	//{
	//	delete crseBCPtr;
	//	crseBCPtr = nullptr;
	//}
	// that should be it!
}

// --------------------------------------------------------------
void CCProjectorComp::applyMacCorrection(LevelData<FluxBox>& a_uEdge,
		Real CFscale,
		int lev)
{

	// this function computes a_uEdge - grad(phi)
	// (assumes all BC's already set)

	//		const DisjointBoxLayout& levelGrids = a_uEdge.getBoxes();
	//		LevelData<FluxBox> gradPhi(levelGrids,1);
	//
	//		LevelData<FArrayBox>* crseBCDataPtr=nullptr;
	//		// should C/F BC's already be set? assume not and reset them here
	////		if (lev > 0)
	////		{
	////			const DisjointBoxLayout crseGrids = m_amrGrids[lev-1];
	////			crseBCDataPtr = new LevelData<FArrayBox>(crseGrids,1);
	////			LevelData<FArrayBox>& crseBCData = *crseBCDataPtr;
	////
	////			Interval comps(0,0);
	////			m_pressure[lev-1]->copyTo(comps,crseBCData,comps);
	////
	////			// now scale the coarse-level BC
	////			DataIterator ditCrse = crseBCData.dataIterator();
	////
	////			for (ditCrse.reset(); ditCrse.ok(); ++ditCrse)
	////			{
	////				crseBCData[ditCrse] *= CFscale;
	////			}
	////		}
	//
	//		// compute gradient
	//		Gradient::levelGradientMAC(*m_gradPressureEdge[lev],*m_pressure[lev], m_pressure[lev-1], m_amrDx[lev],
	//				*m_amrGradIVS[lev],
	//				*m_amrCFInterps[lev]);
	//
	//		// now rescale and subtract from uEdge
	//		DataIterator dit = a_uEdge.dataIterator();
	//
	//		// assumes that uEdge and phi can use same dataIterator
	//		for (dit.reset(); dit.ok(); ++dit)
	//		{
	//			FluxBox& thisGrad = gradPhi[dit];
	//			FluxBox& thisEdgeVel = a_uEdge[dit];
	//
	//			for (int dir=0; dir<SpaceDim; dir++)
	//			{
	//				FArrayBox& thisGradDir = thisGrad[dir];
	//				FArrayBox& thisVelDir = thisEdgeVel[dir];
	//
	//				thisVelDir -= thisGradDir;
	//			}
	//		}
	//
	////		if (crseBCDataPtr!=nullptr)
	////		{
	////			delete crseBCDataPtr;
	////			crseBCDataPtr = nullptr;
	////		}
}

// ---------------------------------------------------------------
void CCProjectorComp::LevelProject(LevelData<FArrayBox>& a_velocity,
		LevelData<FArrayBox>* a_crseVelPtr,
		const Real a_newTime, const Real a_dt)
{
	//todo - write this

	//	LevelData<FArrayBox>* crseDataPtr=nullptr;
	//	LevelData<FArrayBox> levelProjRHS(getBoxes(),1);
	//
	//	// just to be safe.  proper place for this may be outside this function
	//	Interval velComps(0,SpaceDim-1);
	//	a_velocity.exchange(velComps);
	//
	//	if (m_level > 0)
	//	{
	//		if (doQuadInterp())
	//		{
	//			CH_assert(a_crseVelPtr != nullptr);
	//			CH_assert(m_crseProjPtr != nullptr);
	//
	//			// will need to add dt*gradPi from crse vel.
	//			const DisjointBoxLayout& crseLevelGrids = m_crseProjPtr->getBoxes();
	//			crseDataPtr= new LevelData<FArrayBox>(crseLevelGrids, SpaceDim);
	//			LevelData<FArrayBox>& crseData = *crseDataPtr;
	//			m_crseProjPtr->gradPi(crseData);
	//
	//			// now multiply by dt and add crseVel
	//			DataIterator ditCrse = crseData.dataIterator();
	//			LevelData<FArrayBox>& crseVel = *a_crseVelPtr;
	//
	//			// inherent assumption that crseVel and crseGradPi share
	//			// the same layout here
	//			for (ditCrse.reset(); ditCrse.ok(); ++ditCrse)
	//			{
	//				FArrayBox& thisCrseVel = crseVel[ditCrse];
	//				FArrayBox& thisCrseGradPi = crseData[ditCrse];
	//
	//				// note that we multiply by dt on _this_ level
	//				Real scale = a_dt;
	//				scale = 0; //think this is actually correct
	//				thisCrseGradPi *= scale;
	//
	//				// now add crseVel
	//				thisCrseGradPi.plus(thisCrseVel);
	//			}
	//		}
	//		else
	//		{
	//			// alternative is not implemented at this point
	//			MayDay::Error("non-quadratic interpolation BC's not implemented");
	//		}
	//	} // end if coarser level exists
	//
	//	QuadCFInterp velCFInterp;
	//	if (doQuadInterp() && m_crseProjPtr != nullptr)
	//	{
	//		// define two-component CF-interp object for velocities
	//		velCFInterp.define(a_velocity.getBoxes(),
	//				&(crseDataPtr->getBoxes()),
	//				m_dx, m_nRefCrse, SpaceDim, m_domain);
	//	}
	//
	//	// now compute RHS for projection
	//	Divergence::levelDivergenceCC(levelProjRHS, a_velocity, crseDataPtr,
	//			m_dx, doQuadInterp(), velCFInterp);
	//
	//	// for proper scaling of pi, divide this by dt
	//	DataIterator dit = levelProjRHS.dataIterator();
	//	Real dtScale = 1.0/a_dt;
	//	//jparkinson - don't do this scaling, it doesn't work with BCs
	//	dtScale = 1.0;
	//	for (dit.reset(); dit.ok(); ++dit)
	//	{
	//		levelProjRHS[dit] *= dtScale;
	//		//levelProjRHS[dit].setVal(1.4); //Hack! testing if divergence is calculated wrong
	//	}
	//
	//	// set up coarse BC's for solve, then solve
	//	const LevelData<FArrayBox>* crsePiPtr = nullptr;
	//	if (m_crseProjPtr != nullptr) crsePiPtr = &(m_crseProjPtr->Pi());
	//
	//	// report sum(rhs)
	//	if (s_verbosity >= 2)
	//	{
	//		DisjointBoxLayout* finerGridsPtr = nullptr;
	//		int nRefFine = -1;
	//		Real sumRHS = computeSum(levelProjRHS, finerGridsPtr,
	//				nRefFine, m_dx, levelProjRHS.interval());
	//		pout() << "Level projection -- sum(RHS) = " << sumRHS << endl;
	//	}
	//
	//	// now solve for pi
	//	solveMGlevel(m_Pi, crsePiPtr, levelProjRHS);
	//
	//	// apply appropriate physical BC's
	//	BCHolder bcHolder = m_physBCPtr->gradPiFuncBC();
	//	const DisjointBoxLayout& levelGrids = getBoxes();
	//
	//	// loop over grids...
	//	for (dit.reset(); dit.ok(); ++dit)
	//	{
	//		bcHolder.operator()(m_Pi[dit],
	//				levelGrids[dit],
	//				m_domain,
	//				m_dx,
	//				false); // not homogeneous
	//	}
	//
	//	if (m_crseProjPtr != nullptr)
	//	{
	//		// reapply coarse-fine BC's here if necessary
	//		m_cfInterp.coarseFineInterp(m_Pi, *crsePiPtr);
	//		// also clean up after ourselves!
	//		delete crseDataPtr;
	//		crseDataPtr = nullptr;
	//
	//	}
	//
	//	dtScale = a_dt;
	//	dtScale = 1.0; //Don't do this scaling, it doesn't work with BCs
	//	correctCCVelocities(a_velocity, dtScale);

	// check resulting velocity field
	// to do this, need to reset physical BC's
}

// ---------------------------------------------------------------
void CCProjectorComp::correctCCVelocities(LevelData<FArrayBox>& a_velocity,
		const Real scale) const
{
	MayDay::Error("correctCCVelocities - not implemented");
	//	const DisjointBoxLayout& grids = getBoxes();
	//	LevelData<FArrayBox> gradPi(grids, SpaceDim);
	//
	//	// assume that all relevant BC's have already been set
	//	Gradient::levelGradientCC(gradPi, m_Pi, m_dx);
	//
	//	// vel = vel - scale*gradPi
	//	DataIterator dit = gradPi.dataIterator();
	//
	//	for (dit.reset(); dit.ok(); ++dit)
	//	{
	//		FArrayBox& thisGradPi = gradPi[dit];
	//		FArrayBox& thisVel = a_velocity[dit];
	//
	//		thisGradPi *= scale;
	//		thisVel -= thisGradPi;
	//
	//		FArrayBox gradPz(gradPi[dit].box(), 1);
	//		gradPz.copy(gradPi[dit], 1, 0);
	//		//      int temp=1;
	//	}
}

// ---------------------------------------------------------------
void CCProjectorComp::doSyncOperations(Vector<LevelData<FArrayBox>* >& a_velocity,
		//		Vector<LevelData<FArrayBox>* >& a_lambda,
		const Real a_newTime, const Real a_dtLevel)
{
	RelaxSolver<LevelData<FArrayBox> >* bottomSolver = new RelaxSolver<LevelData<FArrayBox> >;
	bottomSolver-> m_verbosity = s_verbosity;

	AMRMultiGrid<LevelData<FArrayBox> >* dspSolver = new
			AMRMultiGrid<LevelData<FArrayBox> >;

	defineMultiGrid(*dspSolver, bottomSolver, false);

	// now call sync projection
	doSyncProjection(a_velocity, a_newTime, a_dtLevel, *dspSolver);

	// Need to delete dspSolver->m_bottomSolver from AMRMultiGrid.
	delete dspSolver;

	//	AMRMultiGrid<LevelData<FArrayBox> >* cvdcSolver = new
	//			AMRMultiGrid<LevelData<FArrayBox> >;
	//	defineMultiGrid(*cvdcSolver, bottomSolver, a_velocity, false);
	//
	//	// now do freestream preservation solve
	//	computeVDCorrection(a_lambda, a_newTime, a_dtLevel, *cvdcSolver);
	//
	//	// Need to delete cvdcSolver->m_bottomSolver from AMRMultiGrid.
	//	delete cvdcSolver;

	// Now delete the bottom solver since AMRMultiGrid won't do this
	delete bottomSolver;
}

// ---------------------------------------------------------------
void CCProjectorComp::doPostRegridOps(Vector<LevelData<FArrayBox>* >& a_velocity,
		Vector<LevelData<FArrayBox>* >& a_lambda,
		const Real a_dt, const Real a_time)
{
	MayDay::Error("doPostRegridOps not implemented");
	//	//
	//	RelaxSolver<LevelData<FArrayBox> >* bottomSolver = new RelaxSolver<LevelData<FArrayBox> >;
	//
	//	AMRMultiGrid<LevelData<FArrayBox> >* bigSolverPtr = new
	//			AMRMultiGrid<LevelData<FArrayBox> >;
	//
	//	defineMultiGrid(*bigSolverPtr, bottomSolver, false);
	//
	//	// for inviscid flow, only do this
	//	// now do freestream preservation solve
	//	computeVDCorrection(a_lambda, a_time, a_dt, *bigSolverPtr);
	//
	//	// Need to delete bigSolverPtr->m_bottomSolver from AMRMultiGrid.
	//	delete bigSolverPtr;
	//
	//	//Now delete the bottom solver since AMRMultiGrid won't do this
	//	delete bottomSolver;
}

// ---------------------------------------------------------------
void CCProjectorComp::doSyncProjection(Vector<LevelData<FArrayBox>* >& a_velocity,
		const Real a_newTime, const Real a_dtSync,
		AMRMultiGrid<LevelData<FArrayBox> >& a_solver
)
{
	m_level = 0;
	//
	//		if (doSyncProjection())
	//		{
	//			// set up temp storage
	//			int vectorSize = a_velocity.size();
	//
	//			Vector<LevelData<FArrayBox>* > syncRHS(vectorSize,nullptr);
	//			Vector<LevelData<FArrayBox>* > syncCorr(vectorSize,nullptr);
	//
	//			int finestLevel = finest_level;
	//
	//			// this is a bogus leveldata pointer for the composite divergences
	//			LevelData<FArrayBox>* crseVelPtr=nullptr;
	//			LevelData<FArrayBox>* fineVelPtr = nullptr;
	//
	//			// loop over levels to allocate storage and compute RHS
	//			for (int lev = m_level; lev<= finestLevel; lev++)
	//			{
	//				const DisjointBoxLayout& levelGrids = a_velocity[lev]->getBoxes();
	//
	//				syncRHS[lev]= new LevelData<FArrayBox>(levelGrids,1);
	//				syncCorr[lev] = &(levelProjPtr->eSync());
	//
	//				Real dx = levelProjPtr->dx();
	//				const ProblemDomain& levelDomain = levelProjPtr->dProblem();
	//				int nRefCrse = levelProjPtr->nRefCrse();
	//				int nRefFine = -1;
	//
	//				if (lev > 0) crseVelPtr = a_velocity[lev-1];
	//				if (lev < finestLevel)
	//				{
	//					fineVelPtr = a_velocity[lev+1];
	//					nRefFine = levelProjPtr->fineProjPtr()->nRefCrse();
	//				}
	//				else
	//				{
	//					fineVelPtr = nullptr;
	//				}
	//
	//				// do exchange here
	//				a_velocity[lev]->exchange(a_velocity[lev]->interval());
	//
	//				// now compute divergence
	//				Divergence::compDivergenceCC(*syncRHS[lev],*a_velocity[lev],
	//						crseVelPtr, fineVelPtr, dx,
	//						nRefCrse, nRefFine, levelDomain,
	//						doQuadInterp());
	//
	//				// increment projection level
	//				levelProjPtr = levelProjPtr->fineProjPtr();
	//
	//			} // end loop over levels for computation of RHS
	//
	//			// take care of C/F Boundary condition for lBase
	//			if (m_level > 0)
	//			{
	//				syncCorr[m_level-1] = &(m_crseProjPtr->eSync());
	//				// need to rescale crse-level BC
	//				LevelData<FArrayBox>& crseSyncCorr = *(syncCorr[m_level-1]);
	//				DataIterator ditCrse = crseSyncCorr.dataIterator();
	//				for (ditCrse.reset(); ditCrse.ok(); ++ditCrse)
	//				{
	//					crseSyncCorr[ditCrse] *= a_dtSync;
	//				}
	//			}
	//
	//			// debugging tool -- want sum(rhs) = 0
	//			Vector<int> nRefFineVect(syncRHS.size());
	//			levelProjPtr = this;
	//			for (int lev=m_level; lev<=finestLevel; lev++)
	//			{
	//				if (lev < finestLevel)
	//				{
	//					nRefFineVect[lev] = levelProjPtr->fineProjPtr()->nRefCrse();
	//				}
	//				else
	//				{
	//					nRefFineVect[lev] = 1;
	//				}
	//				levelProjPtr = levelProjPtr->fineProjPtr();
	//			}
	//
	//			Interval sumComps(0,0);
	//			Real sumRHS;
	//			sumRHS =  computeSum(syncRHS, nRefFineVect, m_dx, sumComps, m_level);
	//			pout() << "Sum(RHS) for sync solve = "
	//					<< setiosflags(ios::scientific) << sumRHS << endl;
	//
	//			// now solve
	//			a_solver.solve(syncCorr, syncRHS, vectorSize-1, m_level,
	//					true); // initialize syncCorr to zero
	//
	//			// apply physical boundary conditions before taking gradients
	//			levelProjPtr = this;
	//
	//			BCHolder bcHolder = m_physBCPtr->gradESyncFuncBC();
	//
	//			for (int lev=m_level; lev<=finestLevel; lev++)
	//			{
	//
	//				// loop over grids
	//				LevelData<FArrayBox>& level_eSync = *syncCorr[lev];
	//				const ProblemDomain& levelDomain = levelProjPtr->dProblem();
	//				Real levelDx = levelProjPtr->dx();
	//				const DisjointBoxLayout& levelGrids = levelProjPtr->getBoxes();
	//				DataIterator ditLev = level_eSync.dataIterator();
	//
	//				for (ditLev.reset(); ditLev.ok(); ++ditLev)
	//				{
	//					bcHolder.operator()(level_eSync[ditLev],
	//							levelGrids[ditLev],
	//							levelDomain,
	//							levelDx,
	//							false); // not homogeneous
	//				}
	//				// exchange here should be unnecessary
	//				//level_eSync.exchange(corrComps);
	//				levelProjPtr = levelProjPtr->fineProjPtr();
	//			}
	//
	//			// now apply sync correction
	//			Real scale = 1.0;
	//			LevelData<FArrayBox>* crseCorrPtr = nullptr;
	//			if (m_level > 0) crseCorrPtr = syncCorr[m_level-1];
	//
	//			applySyncCorrection(a_velocity, scale, crseCorrPtr);
	//
	//			// cleaning up after we're done --
	//			// now need to rescale eSync to look like other pressures
	//			// start with level-1 because we rescaled it for the BC
	//			// also need to delete syncRHS
	//
	//			Real corrScale = 1.0/a_dtSync;
	//			for (int lev=m_level-1; lev<=finestLevel; lev++)
	//			{
	//				if (lev >=0)
	//				{
	//					LevelData<FArrayBox>& levelCorr = *syncCorr[lev];
	//					DataIterator dit = levelCorr.dataIterator();
	//					for (dit.reset(); dit.ok(); ++dit)
	//					{
	//						levelCorr[dit] *= corrScale;
	//					}
	//				} // end if level > 0
	//			} // end loop over levels for rescaling eSync
	//
	//			// loop over levels to delete syncRHS
	//			for (int lev= m_level; lev<=finestLevel; lev++)
	//			{
	//				delete syncRHS[lev];
	//				syncRHS[lev] = nullptr;
	//			}
	//		} // end if do sync in the first place

	//	MayDay::Error("doSyncProjection - not implemented");
}

// ---------------------------------------------------------------
void CCProjectorComp::applySyncCorrection(Vector<LevelData<FArrayBox>* >& a_velocities,
		const Real a_scale,
		LevelData<FArrayBox>* crseCorrPtr)
{
	MayDay::Error("applySyncCorrection - not implemented");
	//	if (m_applySyncCorrection)
	//	{
	//		const DisjointBoxLayout& levelGrids = getBoxes();
	//
	//		LevelData<FArrayBox> levelGradCorr(levelGrids, SpaceDim);
	//		LevelData<FArrayBox>* fineCorrPtr = nullptr;
	//		int nRefFine = -1;
	//
	//		if (!isFinestLevel())
	//		{
	//			fineCorrPtr = &(m_fineProjPtr->eSync());
	//			nRefFine = m_fineProjPtr->nRefCrse();
	//		}
	//
	//		Gradient::compGradientCC(levelGradCorr, m_eSync,
	//				crseCorrPtr, fineCorrPtr,
	//				m_dx, nRefFine, m_domain,
	//				m_gradIVS,
	//				m_cfInterp);
	//
	//		// now loop over grids and correct
	//		LevelData<FArrayBox>& levelVel = *a_velocities[m_level];
	//
	//		DataIterator dit = levelVel.dataIterator();
	//		for (dit.reset(); dit.ok(); ++dit)
	//		{
	//			levelGradCorr[dit] *= a_scale;
	//			levelVel[dit] -= levelGradCorr[dit];
	//		}
	//
	//		// now recursively call finer level
	//		if (!isFinestLevel())
	//		{
	//			m_fineProjPtr->applySyncCorrection(a_velocities, a_scale,
	//					&m_eSync);
	//
	//			// average down finer level to this level
	//			const DisjointBoxLayout& fineGrids = m_fineProjPtr->getBoxes();
	//			int nRefFine = m_fineProjPtr->nRefCrse();
	//			CoarseAverage avgDownThingy(fineGrids, SpaceDim, nRefFine);
	//
	//			avgDownThingy.averageToCoarse(*a_velocities[m_level],
	//					*a_velocities[m_level+1]);
	//		}
	//	}
}

// ---------------------------------------------------------------
void CCProjectorComp::computeVDCorrection(Vector<LevelData<FArrayBox> * >& a_lambda,
		const Real a_newTime, const Real a_dtSync,
		AMRMultiGrid<LevelData<FArrayBox> >& a_solver,
		Vector<LevelData<FluxBox>* >& a_vel
)
{
	CH_assert (a_dtSync != 0);
	s_lambda_timestep = a_dtSync;
	Real lambdaMult = 1.0/a_dtSync;

	Vector<LevelData<FArrayBox>* > VDCorr(finest_level+1);

	// Construct RHS
	for (int lev = 0; lev<=finest_level; lev++)
	{
		// stuff eLambda into VDCorr array
		VDCorr[lev] = m_eLambda[lev];

		// use lambda array as RHS for solve
		LevelData<FArrayBox>& thisLambda = *(a_lambda[lev]);

		DataIterator dit = thisLambda.dataIterator();

		// RHS = (lambda-1)*eta/dtSync
		for (dit.reset(); dit.ok(); ++dit)
		{
			thisLambda[dit] -= 1.0;
			thisLambda[dit] *= lambdaMult;
		}
	}

	// Solve for each level on lev >= this level
	//	for (int lev=0; lev <= finest_level; lev++)
	//	{
	//		// debugging tool -- want sum(rhs) = 0
	//		Interval sumComps(0,0);
	//		Real sumRHS;
	//		sumRHS =  computeSum(a_lambda, m_refinement_ratios, m_amrDx[lev], sumComps, lev);
	//		pout() << "Time = " << a_newTime << " Sum(RHS) for VD solve = "
	//				<< setiosflags(ios::scientific) << sumRHS << endl;
	//
	//	}

	a_solver.solve(VDCorr, a_lambda, finest_level, 0,
			true); // initialize VDCorr to zero


	// apply all bc's here
	BCHolder bcHolder = m_physBCPtr->gradELambdaFuncBC();

	Interval corrComps(0,0);
	for (int lev=0; lev<=finest_level; lev++)
	{
		LevelData<FArrayBox>& levelCorr = *(VDCorr[lev]);
		const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
		DataIterator ditLev = levelCorr.dataIterator();
		for (ditLev.reset(); ditLev.ok(); ++ditLev)
		{
			bcHolder.operator()(levelCorr[ditLev],
					levelGrids[ditLev],
					m_amrDomains[lev],
					m_amrDx[lev],
					false); // not homogeneous
		}
		levelCorr.exchange(corrComps);
	}

	// compute correction field
	computeGrad_eLambda();

	// Add correction to solution
	for (int m_lev = 0; m_lev <= finest_level; m_lev++)
	{
		DataIterator dit = a_vel[m_lev]->dataIterator();
		for (dit.reset(); dit.ok(); ++dit)
		{
//			FArrayBox& eLambda = (*m_eLambda[m_lev])[dit];
			FluxBox& corr = (*m_grad_eLambda[m_lev])[dit];
//			FluxBox& vel = (*a_vel[m_lev])[dit];

//			FluxBox& gradP = (*m_gradPressureEdge[m_lev])[dit];
//			FArrayBox& P = (*m_pressure[m_lev])[dit];

			Box addbox = corr.box();
//			int num_ghost = m_grad_eLambda[m_lev]->ghostVect()[0];
			//			addbox.grow(-(num_ghost-1));

			for (int idir=0; idir<SpaceDim; idir++)
			{
				//				corr[idir].mult(0.5);
			}

			(*a_vel[m_lev])[dit].plus(corr, addbox, 0, 0);
//			int temp=0;
		}
	}

	// Fill ghost cells
	for (int lev = 1; lev <= finest_level; lev++)
	{
		m_amrFillPatchFace[lev]->fillInterp(*a_vel[lev], *a_vel[lev-1], *a_vel[lev-1], 0, 0, 0, 1);
	}

	//Average down to finer levels
	averageCoarseToFineFace(a_vel);

	for (int lev=0; lev<=finest_level; lev++)
	{
		a_vel[lev]->exchange();
	}

	// now return lambda to previous state
	for (int lev = 0; lev<=finest_level; lev++)
	{
		LevelData<FArrayBox>& levelLambda = *(a_lambda[lev]);
		DataIterator dit = levelLambda.dataIterator();

		for (dit.reset(); dit.ok(); ++dit)
		{
			levelLambda[dit] *= a_dtSync;
			levelLambda[dit] += 1.0;
		}
	}

}

// ---------------------------------------------------------------
void CCProjectorComp::computeGrad_eLambda()
{
	// do this recursively starting at the finest level
	//	MayDay::Error("computeGrad_eLambda - not implemented");
	for (int lev = finest_level; lev >=0; lev--)
	{

		// recall that composite MAC gradient is the same as the level
		// gradient, since finer level is not considered to be part
		// of this level (take care of covered regions with avgDown)

		//		if (lev > 0)
		//		{
		//			Gradient::levelGradientMAC(*m_grad_eLambda[lev], *m_eLambda[lev],
		//					m_eLambda[lev-1], m_amrDx[lev],
		//					*m_amrGradIVS[lev],
		//					*m_amrCFInterps[lev]);
		//		}
		//		else
		//		{
		//			QuadCFInterp interp;
		//			Gradient::levelGradientMAC(*m_grad_eLambda[lev], *m_eLambda[lev],
		//					nullptr, m_amrDx[lev],
		//					*m_amrGradIVS[lev],
		//					interp);
		//		}

		LevelData<FArrayBox>* a_phiCrse = nullptr;
		LevelData<FArrayBox>* a_phiFine = nullptr;
		int a_nRefFine;

		if (lev > 0)
		{
			a_phiCrse = m_eLambda[lev-1];
			//a_nRefCrse = m_refinement_ratios[lev-1];
		}
		if (lev < m_finest_level)
		{
			a_phiFine = m_eLambda[lev+1];
			a_nRefFine = m_refinement_ratios[lev];
		}

		//		Box domBox = m_amrDomains[lev].domainBox();

		Gradient::compGradientMAC(*m_grad_eLambda[lev], *m_eLambda[lev],
				&*a_phiCrse, &*a_phiFine, m_amrDx[lev], a_nRefFine,
				*m_amrGradIVS[lev], *m_amrCFInterps[lev]);

		Gradient::compGradientMAC(*m_gradPressureEdge[lev], *m_pressure[lev],
				&*a_phiCrse, &*a_phiFine, m_amrDx[lev], a_nRefFine,
				*m_amrGradIVS[lev], *m_amrCFInterps[lev]);
		//		Gradient::levelGradientMAC(*m_grad_eLambda[lev], *m_eLambda[lev],m_amrDx[lev]);

		if (!s_applyVDCorrection)
		{
			setValLevel(*m_grad_eLambda[lev], 0.0);
		}

		// now average down finer level grad_eLambda if necessary
		if (lev < finest_level)
		{
			const DisjointBoxLayout& fineGrids = m_amrGrids[lev+1];
			int nRefFine = m_refinement_ratios[lev];
			int nComp = 1;

			CoarseAverageEdge avgDownObject(fineGrids, nComp, nRefFine);

			const LevelData<FluxBox> & fineEdgeGrad = *m_grad_eLambda[lev+1];

			avgDownObject.averageToCoarse(*m_grad_eLambda[lev], fineEdgeGrad);
		}

	}
}

// ---------------------------------------------------------------
void CCProjectorComp::rescaleGrad_eLambda(Real a_dx_lbase)
{
	// do this recursively starting at the finest level
	MayDay::Error("rescaleGrad_eLambda - not implemented");
	//	if (!isFinestLevel())
	//	{
	//		m_fineProjPtr->rescaleGrad_eLambda(a_dx_lbase);
	//	}
	//
	//	// now rescale correction on this level.  we rescale correction
	//	// by (dt(lbase)/dt(this_level) = (dx(lbase)/m_dx for consistency
	//	Real scale = a_dx_lbase/m_dx;
	//	DataIterator dit = m_grad_eLambda.dataIterator();
	//	for (dit.begin(); dit.ok(); ++dit)
	//	{
	//		for (int dir=0; dir<SpaceDim; dir++)
	//		{
	//			m_grad_eLambda[dit][dir] *= scale;
	//		}
	//	}
}

// ---------------------------------------------------------------
void CCProjectorComp::initialLevelProject(LevelData<FArrayBox>& a_velocity,
		LevelData<FArrayBox>* a_crseVelPtr,
		const Real a_newTime, const Real a_dt)
{
	MayDay::Error("initialLevelProject - not implemented");
	//	LevelData<FArrayBox>* crseDataPtr=nullptr;
	//	LevelData<FArrayBox> levelProjRHS(getBoxes(),1);
	//
	//	// do coarse-fine BC
	//	if (m_level > 0)
	//	{
	//		if (doQuadInterp())
	//		{
	//			CH_assert(a_crseVelPtr != nullptr);
	//			CH_assert(m_crseProjPtr != nullptr);
	//
	//			// will need to add dt*gradPi to crse vel
	//			const DisjointBoxLayout& crseLevelGrids = m_crseProjPtr->getBoxes();
	//			crseDataPtr = new LevelData<FArrayBox>(crseLevelGrids, SpaceDim);
	//			LevelData<FArrayBox>& crseData = *crseDataPtr;
	//			m_crseProjPtr->gradPi(crseData);
	//
	//			// now multiply by dt and add crseVel
	//			DataIterator ditCrse = crseData.dataIterator();
	//			LevelData<FArrayBox>& crseVel = *a_crseVelPtr;
	//
	//			// inherent assumption that crseVel & crseGradPi share same layout
	//			for (ditCrse.reset(); ditCrse.ok(); ++ditCrse)
	//			{
	//				FArrayBox& thisCrseVel = crseVel[ditCrse];
	//				FArrayBox& thisCrseGradPi = crseData[ditCrse];
	//
	//				// note that we multiply by this level's dt
	//				thisCrseGradPi *= a_dt;
	//				// now add crseVel
	//				thisCrseGradPi.plus(thisCrseVel);
	//			}
	//
	//		}
	//		else
	//		{
	//			// alternative is not implemented at this point
	//			MayDay::Error("non-quadratic interpolation BC's not implemented");
	//		}
	//	} // end if coarser level exists
	//
	//	// just to be safe...
	//	Interval velComps(0,SpaceDim-1);
	//	a_velocity.exchange(velComps);
	//
	//	// now compute RHS for projection
	//	Divergence::levelDivergenceCC(levelProjRHS, a_velocity, crseDataPtr,
	//			m_dx, doQuadInterp(), m_cfInterp);
	//
	//	// for proper scaling of pi, divide this by dt
	//	DataIterator dit = levelProjRHS.dataIterator();
	//	Real dtScale = 1.0/a_dt;
	//
	//	for (dit.reset(); dit.ok(); ++dit)
	//	{
	//		levelProjRHS[dit] *= dtScale;
	//	}
	//
	//	// set up coarse BC's for solve, then solve
	//	const LevelData<FArrayBox>* crsePiPtr = nullptr;
	//	if (m_crseProjPtr != nullptr) crsePiPtr = &(m_crseProjPtr->Pi());
	//
	//	// now solve for Pi
	//	solveMGlevel(m_Pi, crsePiPtr, levelProjRHS);
	//
	//	if (m_crseProjPtr != nullptr)
	//	{
	//		// reapply coarse-fine BC's here if necessary
	//		m_cfInterp.coarseFineInterp(m_Pi, *crsePiPtr);
	//	}
	//
	//	// apply appropriate physical BC's
	//	BCHolder bcHolder = m_physBCPtr->gradPiFuncBC();
	//	const DisjointBoxLayout& levelGrids = getBoxes();
	//
	//	// loop over grids...
	//	for (dit.reset(); dit.ok(); ++dit)
	//	{
	//		bcHolder.operator()(m_Pi[dit],
	//				levelGrids[dit],
	//				m_domain,
	//				m_dx,
	//				false); // not homogeneous
	//	}
	//
	//	Interval piComps(0,0);
	//	m_Pi.exchange(piComps);
	//
	//	dtScale = a_dt;
	//	correctCCVelocities(a_velocity,dtScale);
	//
	//	// clean up storage
	//	if (crseDataPtr!=nullptr)
	//	{
	//		delete crseDataPtr;
	//		crseDataPtr=nullptr;
	//	}
}

// ---------------------------------------------------------------
void CCProjectorComp::initialSyncProjection(Vector<LevelData<FArrayBox>* >& a_vel,
		const Real a_newTime, const Real a_dtSync,
		AMRMultiGrid<LevelData<FArrayBox> >& a_solver
)
{
	// actually, in the current algorithm, i _think_ this is the same
	// as the normal sync projection
	doSyncProjection(a_vel, a_newTime, a_dtSync, a_solver);
}

// ---------------------------------------------------------------
void CCProjectorComp::doInitialSyncOperations(Vector<LevelData<FArrayBox>* >& a_vel,
		Vector<LevelData<FArrayBox>* >& a_lambda,
		const Real a_newTime, const Real a_dtSync)
{
	MayDay::Error("doInitialSyncOperations not implemented");
	//
	//	RelaxSolver<LevelData<FArrayBox> >* bottomSolver = new RelaxSolver<LevelData<FArrayBox> >;
	//	bottomSolver->m_verbosity = s_verbosity;
	//
	//	// now can define multilevel solver
	//	AMRMultiGrid<LevelData<FArrayBox> >* ispSolverPtr = new
	//			AMRMultiGrid<LevelData<FArrayBox> >;
	//
	//	defineMultiGrid(*ispSolverPtr, bottomSolver, false);
	//
	//	// now call sync projection
	//	initialSyncProjection(a_vel, a_newTime, a_dtSync, *ispSolverPtr);
	//
	//	// Need to delete ispSolverPtr->m_bottomSolver from AMRMultiGrid.
	//	delete ispSolverPtr;
	//	AMRMultiGrid<LevelData<FArrayBox> >* cvdcSolverPtr = new
	//			AMRMultiGrid<LevelData<FArrayBox> >;
	//	defineMultiGrid(*cvdcSolverPtr, bottomSolver, true);
	//
	//	// now do freestream preservation solve
	//	computeVDCorrection(a_lambda, a_newTime, a_dtSync, *cvdcSolverPtr);
	//
	//	delete cvdcSolverPtr; //Kris R.
	//
	//	//Now delete the bottom solver since AMRMultiGrid won't do this
	//	delete bottomSolver;

}

void CCProjectorComp::doLambdaCorrection(Vector<LevelData<FluxBox>* >& a_vel,
		Vector<LevelData<FArrayBox> * >& a_lambda,
		Real a_newTime, Real a_dtSync)
{
	//
	RelaxSolver<LevelData<FArrayBox> >* bottomSolver = new RelaxSolver<LevelData<FArrayBox> >;
	bottomSolver->m_verbosity = s_verbosity;

	// now can define multilevel solver
	//	AMRMultiGrid<LevelData<FArrayBox> >* ispSolverPtr = new
	//			AMRMultiGrid<LevelData<FArrayBox> >;

	//defineMultiGrid(*ispSolverPtr, bottomSolver, a_vel, false);

	// now call sync projection
	//initialSyncProjection(a_vel, a_newTime, a_dtSync, *ispSolverPtr);

	// Need to delete ispSolverPtr->m_bottomSolver from AMRMultiGrid.
	//delete ispSolverPtr;

	AMRMultiGrid<LevelData<FArrayBox> >* cvdcSolverPtr = new
			AMRMultiGrid<LevelData<FArrayBox> >;
	defineMultiGrid(*cvdcSolverPtr, bottomSolver, true);

	// now do freestream preservation solve
	computeVDCorrection(a_lambda, a_newTime, a_dtSync, *cvdcSolverPtr, a_vel);

	delete cvdcSolverPtr; //Kris R.

	//Now delete the bottom solver since AMRMultiGrid won't do this
	delete bottomSolver;

}

// ---------------------------------------------------------------
void CCProjectorComp::initialVelocityProject(Vector<LevelData<FArrayBox>* >& a_vel,
		bool a_homogeneousCFBC)
{
	//
	RelaxSolver<LevelData<FArrayBox> >* bottomSolver = new RelaxSolver<LevelData<FArrayBox> >;
	bottomSolver->m_verbosity = s_verbosity;

	// first need to define solver
	AMRMultiGrid<LevelData<FArrayBox> > bigSolver;
	defineMultiGrid(bigSolver, bottomSolver, false);

	initialVelocityProject(a_vel, bigSolver, a_homogeneousCFBC);

	//New delete the bottom solver since AMRMultiGrid won't do this
	delete bottomSolver;


}

// ---------------------------------------------------------------
void CCProjectorComp::initialVelocityProject(Vector<LevelData<FArrayBox>* >& a_vel,
		AMRMultiGrid<LevelData<FArrayBox> >& a_bigSolver,
		bool a_homogeneousCFBC)
{
	MayDay::Error("initialVelocityProject - not implemented");
	//	// set up temp storage
	//	int vectorSize = a_vel.size();
	//
	//	Vector<LevelData<FArrayBox>* > projRHS(vectorSize,nullptr);
	//	Vector<LevelData<FArrayBox>* > projCorr(vectorSize,nullptr);
	//
	//	CCProjectorComp* levelProjPtr = this;
	//
	//	int finestLevel = vectorSize - 1;
	//
	//	// this is a bogus levelData pointer for the composite divergence
	//	LevelData<FArrayBox>* crseVelPtr = nullptr;
	//	LevelData<FArrayBox>* fineVelPtr = nullptr;
	//	Vector<int> nRefFineVect(vectorSize, -1);
	//
	//	// loop over levels to allocate storage and compute RHS
	//	for (int lev=m_level; lev <= finestLevel; lev++)
	//	{
	//		const DisjointBoxLayout& levelGrids = a_vel[lev]->getBoxes();
	//
	//		projRHS[lev] = new LevelData<FArrayBox>(levelGrids,1);
	//		projCorr[lev] = &(levelProjPtr->eSync());
	//
	//		Real dx = levelProjPtr->dx();
	//		const ProblemDomain& levelDomain = levelProjPtr->dProblem();
	//		int nRefCrse = levelProjPtr->nRefCrse();
	//		int nRefFine = -1;
	//
	//		if (lev > 0)
	//		{
	//			crseVelPtr = a_vel[lev-1];
	//		}
	//		if (lev < finestLevel)
	//		{
	//			fineVelPtr = a_vel[lev+1];
	//			nRefFine = levelProjPtr->fineProjPtr()->nRefCrse();
	//		}
	//		else
	//		{
	//			fineVelPtr = nullptr;
	//		}
	//		nRefFineVect[lev] = nRefFine;
	//
	//		// just in case...
	//		Interval velComps(0,SpaceDim-1);
	//		a_vel[lev]->exchange(velComps);
	//
	//		// now compute divergence
	//		Divergence::compDivergenceCC(*projRHS[lev],*a_vel[lev],crseVelPtr,
	//				fineVelPtr,dx,nRefCrse,nRefFine,
	//				levelDomain,doQuadInterp());
	//
	//		// increment levelProjPtr
	//		levelProjPtr= levelProjPtr->fineProjPtr();
	//
	//	} // end loop over levels for computation of RHS
	//
	//	// take care of C/F Boundary condition for lBase
	//
	//	LevelData<FArrayBox>* crseBCPtr = nullptr;
	//	if (m_level > 0)
	//	{
	//		CH_assert(m_crseProjPtr != nullptr);
	//		if (a_homogeneousCFBC)
	//		{
	//			// need to define coarse BC and set to 0
	//			const DisjointBoxLayout& crseGrids = m_crseProjPtr->getBoxes();
	//
	//			crseBCPtr = new LevelData<FArrayBox>(crseGrids,1);
	//
	//			setValLevel(*crseBCPtr, 0.0);
	//
	//			projCorr[m_level-1] = crseBCPtr;
	//		}
	//		else
	//		{
	//			// not sure if this should ever be called
	//			MayDay::Warning("inhomogeneous BC's on initialVelProject");
	//			projCorr[m_level-1] = &m_crseProjPtr->eSync();
	//		}
	//
	//	} // end if level > 0
	//
	//	Vector<DisjointBoxLayout> allGrids(vectorSize);
	//	for (int lev=m_level; lev<=finestLevel; lev++)
	//	{
	//		allGrids[lev] = projCorr[lev]->getBoxes();
	//	}
	//
	//	Interval sumComps(0,0);
	//	Real sumRHS;
	//	sumRHS =  computeSum(projRHS, nRefFineVect, m_dx, sumComps, m_level);
	//	pout() << "Sum(RHS) for IVP solve = "
	//			<< setiosflags(ios::scientific) << sumRHS << endl;
	//
	//	// now solve!
	//	a_bigSolver.solve(projCorr, projRHS, finestLevel, m_level,
	//			true); // initialize projCorr to zero
	//
	//	// apply physical boundary conditions before
	//	// taking gradients
	//	Interval corrComps(0,0);
	//	levelProjPtr = this;
	//	BCHolder bcHolder = m_physBCPtr->gradESyncFuncBC();
	//
	//	for (int lev=m_level; lev<=finestLevel; lev++)
	//	{
	//
	//		// loop over grids
	//		LevelData<FArrayBox>& levelCorr = *projCorr[lev];
	//		const ProblemDomain& levelDomain = levelProjPtr->dProblem();
	//		Real levelDx = levelProjPtr->dx();
	//		const DisjointBoxLayout& levelGrids = allGrids[lev];
	//		// levelProjPtr->getBoxes();
	//		DataIterator ditLev = levelCorr.dataIterator();
	//
	//		for (ditLev.reset(); ditLev.ok(); ++ditLev)
	//		{
	//			bcHolder.operator()(levelCorr[ditLev],
	//					levelGrids[ditLev],
	//					levelDomain,
	//					levelDx,
	//					false); // not homogeneous
	//		}
	//
	//		levelCorr.exchange(corrComps);
	//
	//		levelProjPtr = levelProjPtr->fineProjPtr();
	//	}
	//
	//	// now apply correction
	//	Real scale = 1.0;
	//
	//	applySyncCorrection(a_vel, scale, crseBCPtr);
	//
	//	// delete temp storage -- RHS and crse BC
	//	for (int lev=m_level; lev<=finestLevel; lev++)
	//	{
	//		delete projRHS[lev];
	//		projRHS[lev] = nullptr;
	//	}
	//
	//	if (crseBCPtr != nullptr)
	//	{
	//		delete crseBCPtr;
	//		crseBCPtr = nullptr;
	//	}
}

// ---------------------------------------------------------------
void CCProjectorComp::defineMultiGrid(AMRMultiGrid<LevelData<FArrayBox> >& a_solver,
		LinearSolver<LevelData<FArrayBox> >* a_bottomSolver,
		bool a_freestreamSolve)
{
	//	MayDay::Error("defineMultiGrid - not implemented");
	// assume that this level is lBase...


	// determine finest existing level
	//		CCProjectorComp* levelProj = this;

	//		int finestLevel = finest_level;

	// put level grids togther into a vector
	//		Vector<DisjointBoxLayout> allGrids(finest_level+1);
	//		Vector<ProblemDomain> allDomains(finest_level+1);
	//		Vector<Real> allDx(finest_level+1);
	//		Vector<int> allRefRatios(finest_level+1);

	//		levelProj = this; // Now levelProj is at level m_level.

	//Since we are always going to need the grids on this level do that now
	//		allGrids[m_level] = a_vel[m_level]->getBoxes();

	//		for (int lev = m_level+1; lev<=finestLevel; lev++)
	//		{
	//			levelProj = levelProj->fineProjPtr(); //levelProj now points to the same level as index lev
	//			allGrids[lev] = m_pressure[lev]->getBoxes();
	//			allRefRatios[lev-1] = levelProj->nRefCrse(); //ref ratio to next finer level belongs to that level for CCProjectorComp
	//
	//		}
	//		// Now levelProj is at level finestLevel.
	//
	//		//Kris R.
	//		levelProj = this; //now levelProj is at m_level

	//now deal with the case where there are coarser levels than this one
	//		for (int lev = m_level-1; lev >=0; lev--) {
	//
	//			allRefRatios[lev] = levelProj->nRefCrse(); //levelProj points to the level above lev
	//
	//			levelProj = levelProj->crseProjPtr(); //lecvelProj now points to the same level as lev
	//			allGrids[lev] = m_pressure[lev]->getBoxes();
	//		}

	//levelProj is now at level 0. This is crucial for making the LHS operator consistent
	//		CH_assert(levelProj->m_level == 0);

	//		ProblemDomain baseDomain = levelProj->dProblem(); //Kris R.
	//		Real baseDx = levelProj->dx(); //Kris R.

	AMRPoissonOpFactory localPoissonOpFactory;
	// physBCPtr = m_physBCPtr->FreestreamCorrBC();
	// physBCPtr = m_physBCPtr->SyncProjBC();
	localPoissonOpFactory.define(m_amrDomains[0],
			m_amrGrids,
			m_refinement_ratios,
			m_amrDx[0],
			((a_freestreamSolve) ?
					m_physBCPtr->FreestreamCorrFuncBC() :
					m_physBCPtr->SyncProjFuncBC() ));

	// You really need to delete this when you're done with a_solver.
	// But m_bottomSolver is a protected field of AMRMultiGrid.
	RelaxSolver<LevelData<FArrayBox> >* bottomSolverPtr = new
			RelaxSolver<LevelData<FArrayBox> >;
	bottomSolverPtr->m_verbosity = s_verbosity;

	// AMRMultiGrid<LevelData<FArrayBox> >& a_solver
	a_solver.define(m_amrDomains[0],
			localPoissonOpFactory,
			a_bottomSolver,
			finest_level+1);
	// We also want to use m_limitCoarsening:
	// if true, only do multigrid coarsening down to next coarser
	// AMR level (only coarsen by a_nRefCrse).
	// If false, coarsen as far as possible. Only relevant when lBase > 0.
	a_solver.m_verbosity = s_verbosity;
	a_solver.m_eps = 1e-10;
	if (s_solver_tolerance > 0)
	{
		a_solver.m_eps = s_solver_tolerance;
	}
	a_solver.m_pre = s_num_smooth_down; // smoothings before avging
	a_solver.m_post = s_num_smooth_up; // smoothings after avging
}


void CCProjectorComp::defineSolverMG()
{


	Vector<RefCountedPtr<LevelData<FArrayBox> >  > aCoef;
	aCoef.resize(finest_level+1);

	for (int lev=0; lev<=finest_level; lev++)
	{
		aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_amrGrids[lev], 1, num_ghost*IntVect::Unit));
		setValLevel(*aCoef[lev], 0.0);
	}

	Real alpha = 0, beta = -1;

	m_solverOpFact.define(m_amrDomains[0],
				m_amrGrids,
				m_refinement_ratios,
				m_amrDx[0],
				m_physBCPtr->LevelPressureFuncBC(),
				alpha, aCoef,
				beta, m_permeability);

	m_solverMG.define(m_amrDomains[0],
			m_solverOpFact,
			m_bottomSolver,
			finest_level + 1);

	m_solverMG.m_verbosity = 2;
	m_solverMG.m_eps = 1e-10;
	if (s_solver_tolerance > 0)
	{
		m_solverMG.m_eps = s_solver_tolerance;
	}
	m_solverMG.m_pre = s_num_smooth_down; // smoothings before avging
	m_solverMG.m_post = s_num_smooth_up; // smoothings after avging
	m_solverMG.m_iterMax = 300;
}

// -------------------------------------------------------------
void CCProjectorComp::defineSolverMGlevel(const DisjointBoxLayout& a_grids,
		Vector<RefCountedPtr<LevelData<FluxBox> > > a_permeability,
		Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef,
		const DisjointBoxLayout* a_crseGridsPtr,
		const int lev)
{
	AMRPoissonOpFactory localPoissonOpFactory;
	VCAMRPoissonOp2Factory vcPoissonOpFactory;

	BCHolder pressureBC = m_physBCPtr->LevelPressureFuncBC();

	Real alpha=0, beta=-1;

	ProblemDomain baseDomain(m_amrDomains[lev]); // on this level
	int numSolverLevels;
	if (a_crseGridsPtr != nullptr)
	{ // coarser level exists:  define solver on two levels
		numSolverLevels = 2;
		// this returns a null domain for me:
		// baseDomain = m_crseProjPtr->dProblem();
		baseDomain.coarsen(m_refinement_ratios[lev-1]);
		// baseDomain = ProblemDomain(boundingBox(*a_crseGridsPtr));
		// LayoutIterator litFirst = a_crseGridsPtr->layoutIterator();
		// const Box& bxFirst = a_crseGridsPtr->operator[](litFirst);
		// Box bxBound = bxFirst;
		// for (LayoutIterator lit = a_crseGridsPtr->layoutIterator(); lit.ok(); ++lit)
		// {
		// const Box& thisBox = a_crseGridsPtr->operator[](lit);
		// bxBound.minBox(thisBox);
		// }
		// baseDomain = ProblemDomain(bxBound);
		// }
		Vector<DisjointBoxLayout> allGrids(2);
		allGrids[0] = *a_crseGridsPtr;
		allGrids[1] = a_grids;
		Vector<int> refRatios(1, m_refinement_ratios[lev-1]);
		// this returns zero for me:
		// Real dxCrse = m_crseProjPtr->dx();
		Real dxCrse = m_refinement_ratios[lev-1] * m_amrDx[lev];
		localPoissonOpFactory.define(baseDomain,
				allGrids,
				refRatios,
				dxCrse,
				pressureBC);

		vcPoissonOpFactory.define(baseDomain,
				allGrids,
				refRatios,
				dxCrse,
				pressureBC,
				alpha, aCoef,
				beta, a_permeability);
	}
	else
	{ // no coarser level:  define solver on only one level
		numSolverLevels = 1;
		localPoissonOpFactory.define(m_amrDomains[lev],
				a_grids,
				m_amrDx[lev],
				pressureBC);

		Vector<DisjointBoxLayout> allGrids(1);
		allGrids[0] = a_grids;

		Vector<int> refRatios(1, 2);

		vcPoissonOpFactory.define(baseDomain,
				allGrids,
				refRatios,
				m_amrDx[0],
				pressureBC,
				alpha, aCoef,
				beta, a_permeability);
	}

	// You really need to delete this when you're done with a_solver.
	// But m_bottomSolver is a protected field of AMRMultiGrid.

	//Delete previous bottom solver -- Kris R.
	if(m_bottomSolverLevel != nullptr)
	{
		delete m_bottomSolverLevel;
		m_bottomSolverLevel = nullptr;
	}

	//		RelaxSolver<LevelData<FArrayBox> >* newBottomPtr = new RelaxSolver<LevelData<FArrayBox> >;
//	RelaxSolver<LevelData<FArrayBox> >* newBottomPtr = new RelaxSolver<LevelData<FArrayBox> >;
	BiCGStabSolver<LevelData<FArrayBox> >* newBottomPtr = new BiCGStabSolver<LevelData<FArrayBox> >;
	newBottomPtr->m_verbosity = 0;
	m_bottomSolverLevel = newBottomPtr;

	// AMRMultiGrid<LevelData<FArrayBox> >& a_solver
	m_solverMGlevel.define(baseDomain, // on either this level or coarser level
//			localPoissonOpFactory,
			vcPoissonOpFactory,
			m_bottomSolverLevel,
			numSolverLevels);
	m_solverMGlevel.m_verbosity = s_verbosity;
	m_solverMGlevel.m_eps = 1e-15; // was 1e-10
	if (s_solver_tolerance > 0)
	{
		m_solverMGlevel.m_eps = s_solver_tolerance;
	}
	m_solverMGlevel.m_pre = s_num_smooth_down; // smoothings before avging
	m_solverMGlevel.m_post = s_num_smooth_up; // smoothings after avging
	m_solverMGlevel.m_iterMax = 30;
}


// -------------------------------------------------------------
void CCProjectorComp::solveMGlevel(LevelData<FArrayBox>&   a_phi,
		const LevelData<FArrayBox>*   a_phiCoarsePtr,
		const LevelData<FArrayBox>&   a_rhs)
{

	// Initialize a_phi to zero.
	setValLevel(a_phi, 0.0);
	Vector< LevelData<FArrayBox>* > phiVect;
	Vector< LevelData<FArrayBox>* > rhsVect;
	LevelData<FArrayBox>& rhsRef =
			const_cast< LevelData<FArrayBox>& >(a_rhs);
	int maxLevel;
	if (a_phiCoarsePtr != nullptr)
	{
		maxLevel = 1;
		LevelData<FArrayBox>* phiCoarsePtrRef =
				const_cast< LevelData<FArrayBox>* >(a_phiCoarsePtr);
		phiVect.push_back(phiCoarsePtrRef);
		rhsVect.push_back(nullptr); // I don't think this will be used
	}
	else
	{
		maxLevel = 0;
	}
	phiVect.push_back(&a_phi);
	rhsVect.push_back(&rhsRef);

	// l_max = maxLevel, l_base = maxLevel
	m_solverMGlevel.solve(phiVect, rhsVect, maxLevel, maxLevel,
			false); // don't initialize to zero
	//  int exitStatus = m_solverMGlevel.m_exitStatus;
}

///Assumes BCs are set on U^*
void
CCProjectorComp::projectVelocity(Vector<LevelData<FArrayBox> *> a_U,
		Vector<RefCountedPtr<LevelData<FArrayBox> > > a_permeability,
		Vector<LevelData<FluxBox>* > a_Uedge,
		const Vector<LevelData<FArrayBox> *> a_Ustar,
		const int order)
{
	bool enforceUstar = false;
	bool enforceDivUstar = false;
	bool enforceP = false;
	bool enforceGradP = false;
	bool enforceUstarEdge = false;
	Real rayleighTemp = 100; Real perturbation = 0.1;

	EdgeVelBCHolder edgeVelBC(m_physBCPtr->edgeVelFuncBC(false));

	Vector<LevelData<FluxBox>*> a_UstarFace;
	a_UstarFace.resize(finest_level+1, nullptr);

	for (int lev=0; lev<=finest_level; lev++)
	{
		a_UstarFace[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, num_ghost*IntVect::Unit);
		CellToEdge((*a_Ustar[lev]), (*a_UstarFace[lev]));

		edgeVelBC.applyBCs(*a_UstarFace[lev], m_amrGrids[lev],
							m_amrDomains[lev], m_amrDx[lev],
							false); // inhomogeneous
	}

	if (enforceUstar)
	{
		for (int lev=0; lev<=finest_level; lev++)
		{
			for (DataIterator dit = a_Ustar[lev]->dataIterator(); dit.ok(); ++dit)
			{
				Box b = (*a_Ustar[lev])[dit].box();
				Box interiorBox = b;
				//					interiorBox.grow(-(num_ghost+1));

				b.grow(-num_ghost);
				int idir=1;

				Box lo = adjCellLo(b, idir);
				Box hi = adjCellHi(b, idir);

				// Fill ghost cells in y direction
				for (BoxIterator bit=BoxIterator(lo); bit.ok(); ++bit)
				{
					IntVect iv = bit();

					Real v_i = (*a_Ustar[lev])[dit](iv, idir);
					Real v_ii = (*a_Ustar[lev])[dit](iv + BASISV(idir), idir);
					Real v_iii = (*a_Ustar[lev])[dit](iv + 2*BASISV(idir), idir);

					(*a_Ustar[lev])[dit](iv - BASISV(idir), 0) = 3*(v_i - v_ii) + v_iii;
				}

				for (BoxIterator bit=BoxIterator(hi); bit.ok(); ++bit)
				{
					IntVect iv = bit();

					Real v_i = (*a_Ustar[lev])[dit](iv, idir);
					Real v_ii = (*a_Ustar[lev])[dit](iv - BASISV(idir), idir);
					Real v_iii = (*a_Ustar[lev])[dit](iv -  2*BASISV(idir), idir);

					(*a_Ustar[lev])[dit](iv + BASISV(idir), 0) = 3*(v_i - v_ii) + v_iii;
				}


				b.grow(num_ghost);

				FArrayBox err(b, 1);
				err.copy((*a_Ustar[lev])[dit], 1, 0, 1);

				for (BoxIterator bit(b); bit.ok(); ++bit)
				{
					IntVect iv = bit();

					//Only fills ghost cells for now
					//						if (!interiorBox.contains(iv))
					//						{

					RealVect loc = iv;
					loc *= m_amrDx[lev];
					loc += 0.5*m_amrDx[lev]*RealVect::Unit;
					Real x = loc[0];
					Real y = loc[1];


					(*a_Ustar[lev])[dit](iv, 0) = 0;
					(*a_Ustar[lev])[dit](iv, 1) = rayleighTemp * (y + perturbation*cos(M_PI*(x))*sin(M_PI*(y)));

					err(iv, 0) -= rayleighTemp * (y + perturbation*cos(M_PI*(x))*sin(M_PI*(y)));
					//						}

				}

//				int temp=0;
			}
		}

	}

	//First calculate div(U^*)
	LevelData<FArrayBox>* a_uCrsePtr;
//	LevelData<FArrayBox>* a_uFinePtr;
	LevelData<FluxBox>* a_uFaceFinePtr = nullptr;
	Real* dxFine = nullptr;

	for (int lev=0; lev<=finest_level; lev++)
	{
		a_uCrsePtr = nullptr; // a_uFinePtr = nullptr;
		int refFine = 2; //refCrse = 2,

		if (lev > 0)
		{
			a_uCrsePtr = a_Ustar[lev-1];
//			refCrse = m_refinement_ratios[lev-1];
		}

		if (lev < finest_level)
		{
//			a_uFinePtr = a_Ustar[lev+1];
			a_uFaceFinePtr = a_UstarFace[lev+1];
			refFine = m_refinement_ratios[lev];
			dxFine = &(m_amrDx[lev+1]);
		}
		else
		{
			 a_uFaceFinePtr = nullptr;
		}



		QuadCFInterp velCFInterp;
		if (lev > 0)
		{
			// define two-component CF-interp object for velocities
			velCFInterp.define(a_Ustar[lev]->getBoxes(),
					&(a_uCrsePtr->getBoxes()),
					m_amrDx[lev], m_refinement_ratios[lev-1], SpaceDim, m_amrDomains[lev]);
		}

		Divergence::compDivergenceMAC(*m_divUstar[lev], *a_UstarFace[lev], a_uFaceFinePtr,
				m_amrDx[lev], dxFine, refFine, m_amrDomains[lev]);

	}

	if (enforceDivUstar)
	{
		for (int lev=0; lev<=finest_level; lev++)
		{
			for (DataIterator dit = m_divUstar[lev]->dataIterator(); dit.ok(); ++dit)
			{
				Box b = (*m_divUstar[lev])[dit].box();
				FArrayBox err(b, 1);
				for (BoxIterator bit(b); bit.ok(); ++bit)
				{
					IntVect iv = bit();
					RealVect loc = iv;
					loc *= m_amrDx[lev];
					loc += 0.5*m_amrDx[lev]*RealVect::Unit;
					Real x = loc[0];
					Real y = loc[1];
					Real div = rayleighTemp *(-1+perturbation*M_PI*cos(M_PI*(x))*cos(M_PI*(y)));

					err(iv, 0) = (*m_divUstar[lev])[dit](iv, 0) - div;

					(*m_divUstar[lev])[dit](iv, 0) =  div;
				}

//				int temp=0;
			}
		}
	}

	//Then solve for Pressure

	// Initialise to zero
	//todo - remove this!
//	for (int lev = 0; lev <= finest_level; lev++)
//	{
//		setValLevel(*m_pressure[lev], 0.0);
//	}

	Real sum = computeSum(m_divUstar, m_refinement_ratios, m_amrDx[0]);
	pout() << "     CCProjectorComp - sum of div(U^*) = " << sum << endl;

//	for (DataIterator dit = m_pressure[0]->dataIterator(); dit.ok(); ++dit)
//				{
//					FArrayBox& divU = (*m_divUstar[0])[dit];
//					FluxBox& perm = (*m_permeability[0])[dit];
//					int temp=0;
//
//				}

	m_solverMG.solve(m_pressure, m_divUstar, finest_level, 0, false);

	if (enforceP)
	{
		for (int lev=0; lev<=finest_level; lev++)
		{
			for (DataIterator dit = m_pressure[lev]->dataIterator(); dit.ok(); ++dit)
			{
				Box b = (*m_pressure[lev])[dit].box();
				for (BoxIterator bit(b); bit.ok(); ++bit)
				{
					IntVect iv = bit();
					RealVect loc = iv;
					loc *= m_amrDx[lev];
					loc += 0.5*m_amrDx[lev]*RealVect::Unit;
					Real x = loc[0];
					Real y = loc[1];
					(*m_pressure[lev])[dit](iv,0) = rayleighTemp *(y-0.5*y*y - (perturbation/(2*M_PI))*cos(M_PI*(x))*cos(M_PI*(y)));
				}
			}
		}
	}

	//Enforce physical BCs on Pressure
	//	BCHolder bcHolder = m_physBCPtr->gradMacPressureFuncBC();
	//	for (int lev=0; lev<=finest_level; lev++)
	//	{
	//	  const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
	//
	//	  DataIterator dit = m_pressure[lev]->dataIterator();
	//	  for (dit.reset(); dit.ok(); ++dit)
	//	    {
	//	      bcHolder.operator()((*m_pressure[lev])[dit],
	//	                          levelGrids[dit],
	//	                          m_amrDomains[lev],
	//	                          m_amrDx[lev],
	//	                          false); // not homogeneous
	//	    }
	//	}


	//Then take grad(P)

	for (int lev=0; lev <= finest_level; lev++)
	{
		//		 Gradient::levelGradientCCNewCentred(*m_gradPressure[lev],
		//						 *m_pressure[lev],
		//				    				m_amrDx[lev]);

		LevelData<FArrayBox>* a_phiCrse = nullptr;
		LevelData<FArrayBox>* a_phiFine = nullptr;
		int a_nRefFine; //a_nRefCrse

		if (lev > 0)
		{
			a_phiCrse = m_pressure[lev-1];
//			a_nRefCrse = m_refinement_ratios[lev-1];
		}
		if (lev < m_finest_level)
		{
			a_phiFine = m_pressure[lev+1];
			a_nRefFine = m_refinement_ratios[lev];
		}

		//		Box domBox = m_amrDomains[lev].domainBox();

		if (a_phiFine)
		{
                  Gradient::compGradientMAC(*m_gradPressureEdge[lev], *m_pressure[lev],
                                  &*a_phiCrse, &*a_phiFine, m_amrDx[lev], a_nRefFine,
                                  *m_amrGradIVS[lev], *m_amrCFInterps[lev]);
		}
		else if (a_phiCrse)
		{
		  Gradient::levelGradientMAC(*m_gradPressureEdge[lev], *m_pressure[lev],
		                             &*a_phiCrse, m_amrDx[lev],
		                             *m_amrGradIVS[lev], *m_amrCFInterps[lev]);
		}
		else
		{
		  Gradient::levelGradientMAC(*m_gradPressureEdge[lev], *m_pressure[lev],
		                             m_amrDx[lev]);
		}




		//		Gradient::levelGradientMACNew(*m_gradPressureEdge[lev], *m_pressure[lev], m_amrDx[lev]);

		//m_amrDomains[lev]);


//		for (DataIterator dit = m_gradPressureEdge[lev]->dataIterator(); dit.ok(); ++dit)
//		{
//
//			FluxBox& gradP = (*m_gradPressureEdge[lev])[dit];
//
//
//			int temp=0;
//
//		}

	} // end loop over levels


	//todo - remove testing - enforcing analytic grad(P)
	if (enforceGradP)
	{
		for (int lev=0; lev<=finest_level; lev++)
		{
			for (DataIterator dit = m_gradPressure[lev]->dataIterator(); dit.ok(); ++dit)
			{

				//				Box b = (*m_gradPressure[lev])[dit].box();

				for (int idir=0; idir<SpaceDim; idir++)
				{
					Box b = (*m_gradPressureEdge[lev])[dit][idir].box();


					for (BoxIterator bit(b); bit.ok(); ++bit)
					{
						IntVect iv = bit();
						RealVect loc = iv;
						loc *= m_amrDx[lev];
						loc += 0.5*m_amrDx[lev]*RealVect::Unit;
						Real xEdge = loc[0];
						Real yEdge = loc[1];



						Real gradPedge = 0;
						if (idir == 0)
						{
							xEdge = xEdge - m_amrDx[lev]/2;
							gradPedge = rayleighTemp *(0.5*perturbation*sin(M_PI*(xEdge))*cos(M_PI*(yEdge)));
						}
						else if (idir == 1)
						{
							yEdge = yEdge - m_amrDx[lev]/2;
							gradPedge = rayleighTemp *(1-yEdge+0.5*perturbation*cos(M_PI*(xEdge))*sin(M_PI*(yEdge)));
						}
						(*m_gradPressureEdge[lev])[dit][idir](iv,0) = gradPedge;

					}
				}
//				int temp=0;
			}
		}
	} // end if enforce grad(P)



	for (int lev=0; lev<=finest_level; lev++)
	{
		//Try enforcing Ustaredge
		if (enforceUstarEdge)
		{
			for (DataIterator dit = a_UstarFace[lev]->dataIterator(); dit.ok(); ++dit)
			{

				//				Box b = (*m_gradPressure[lev])[dit].box();

				for (int idir=0; idir<SpaceDim; idir++)
				{
					Box b = (*a_UstarFace[lev])[dit][idir].box();

					for (BoxIterator bit(b); bit.ok(); ++bit)
					{
						IntVect iv = bit();
						RealVect loc = iv;
						loc *= m_amrDx[lev];
						loc += 0.5*m_amrDx[lev]*RealVect::Unit;
						Real xEdge = loc[0];
						Real yEdge = loc[1];

						Real uEdgeDir = 0;
						if (idir == 0)
						{
							uEdgeDir = 0;
						}
						else if (idir == 1)
						{
							yEdge = yEdge - m_amrDx[lev]/2;
							Real theta = 1-yEdge + perturbation*cos(M_PI*xEdge)*sin(M_PI*yEdge);
							uEdgeDir = rayleighTemp *theta;
						}

						(*a_UstarFace[lev])[dit][idir](iv, 0) = uEdgeDir;


					}
				}
			}
		}

	}


	//Finally construct actual velocity U = - permeability*grad(P) + U^*
	for (int lev=0; lev<=finest_level; lev++)
	{
		for (DataIterator dit = m_gradPressure[lev]->dataIterator(); dit.ok(); ++dit)
		{
			//			(*a_U[lev])[dit].setVal(0);
			//			(*a_U[lev])[dit] -= (*m_gradPressure[lev])[dit];
			//			(*a_U[lev])[dit] += (*a_Ustar[lev])[dit];
//			FluxBox& gradP = (*m_gradPressureEdge[lev])[dit];
//			FArrayBox& P = (*m_pressure[lev])[dit];

			(*a_Uedge[lev])[dit].setVal(0);  // U=0
			(*a_Uedge[lev])[dit] -= (*m_gradPressureEdge[lev])[dit];  // U=-grad(P)
			(*a_Uedge[lev])[dit].mult((*m_permeability[lev])[dit], (*m_permeability[lev])[dit].box(), 0,0); // U = - permeability*grad(P)

			(*a_Uedge[lev])[dit] += (*a_UstarFace[lev])[dit]; // U = -permeability*grad(P) + Ra*theta
//			int temp=0;
		}
	}

	for (int lev=0; lev<=finest_level; lev++)
	{
		//Apply BCs
		edgeVelBC.applyBCs(*a_Uedge[lev], m_amrGrids[lev],
				m_amrDomains[lev], m_amrDx[lev],
				false); // inhomogeneous


		//Fill ghost cells
		if (lev > 0)
		{
			//NodeQCFI qfInterp(m_amrGrids[lev], m_amrDx[lev], m_amrDomains[lev], loCFIVS, hiCFIVS, m_refinement_ratios[lev-1], )
			int interpRadius = 2; // Must be integer multiple of refinement ratio
			int nComps = 1;
			PiecewiseLinearFillPatchFace filler(m_amrGrids[lev], m_amrGrids[lev-1], nComps, m_amrDomains[lev-1], m_refinement_ratios[lev-1], interpRadius);
			filler.fillInterp(*a_Uedge[lev], *a_Uedge[lev-1],*a_Uedge[lev-1],
					0, // time interp coeff
					0, 0, 1);

		}

		a_Uedge[lev]->exchange();

	}

	// This seems to make things worse?
	averageCoarseToFineFace(a_Uedge);

	for (int lev=0; lev<=finest_level; lev++)
	{
		for (int idir=0; idir<SpaceDim; idir++)
		{
			//			LevelData<FArrayBox>* a_UCCdir(a_U[lev]->disjointBoxLayout(), 1);
			EdgeToCell((*a_Uedge[lev]), (*a_U[lev]));
		}
	}


	//Clean up memory
	for (int lev=0; lev<=m_finest_level; lev++)
	{

		if(a_UstarFace[lev] != nullptr)
		{
			delete a_UstarFace[lev];
			a_UstarFace[lev]= nullptr;
		}
	}
}

void CCProjectorComp::averageCoarseToFineFace(Vector<LevelData<FluxBox>* > a_phi)
{
	for (int lev = finest_level-1; lev >= 0; lev--)
	{

		const DisjointBoxLayout& fineGrids = m_amrGrids[lev+1];
		int nRefFine = m_refinement_ratios[lev];
		int nComp = 1;

		CoarseAverageEdge avgDownObject(fineGrids, nComp, nRefFine);

		avgDownObject.averageToCoarse(*a_phi[lev], *a_phi[lev+1]);
	}

}

void CCProjectorComp::getDivUStar(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_div)
{
	for (int lev=0; lev<=finest_level; lev++)
	{
		m_divUstar[lev]->copyTo(*a_div[lev]);
	}
}
