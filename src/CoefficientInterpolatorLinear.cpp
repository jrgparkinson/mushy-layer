#include "CoefficientInterpolatorLinear.H"
#include "timeInterp.H"
#include "Interval.H"


 void CoefficientInterpolatorLinear::define(RefCountedPtr<LevelData<FArrayBox> > a_coefOld,
		 RefCountedPtr<LevelData<FArrayBox> > a_coefNew,
	   Real a_timeOld, Real a_timeNew)
  {
	  m_coefOld = a_coefOld;
	  m_coefNew = a_coefNew;
	  m_timeOld = a_timeOld;
	  m_timeNew = a_timeNew;

	  CoefficientInterpolator(m_coefOld->nComp());

  }

void CoefficientInterpolatorLinear::interpolate(LevelData<FArrayBox>& a_result,
        Real a_time)
{
	Interval interv(0,numComps()-1);
	timeInterp(a_result, a_time, *m_coefOld, m_timeOld,
					*m_coefNew, m_timeNew, interv);
}

void CoefficientInterpolatorLinear::interpolate(LevelData<FArrayBox>& a_result,
	                           const LevelData<FArrayBox>& a_solution,
	                           Real a_time)
{
	interpolate(a_result, a_time);
}

bool CoefficientInterpolatorLinear::dependsUponSolution() const
{
	return false;
}

void CoefficientInterpolatorLinear::interpolatePrime(LevelData<FArrayBox>& a_prime, const LevelData<FArrayBox>& a_solution, Real a_time)
{
MayDay::Error("Not implemented");
}

void CoefficientInterpolatorLinear::solve(LevelData<FArrayBox>&a_phi, const LevelData<FArrayBox>& a_f,
		Real a_time, const LevelData<FArrayBox>& a_phi0, Real a_tolerance)
{
	MayDay::Error("Not implemented");
}
