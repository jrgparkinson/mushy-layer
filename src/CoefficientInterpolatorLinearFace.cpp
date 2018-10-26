#include "CoefficientInterpolatorLinearFace.H"
#include "timeInterp.H"
#include "Interval.H"
#include "FArrayBox.H"
#include "FluxBox.H"


 void CoefficientInterpolatorLinearFace::define(RefCountedPtr<LevelData<FluxBox> > a_coefOld,
		 RefCountedPtr<LevelData<FluxBox> > a_coefNew,
	   Real a_timeOld, Real a_timeNew)
  {
	  m_coefOld = a_coefOld;
	  m_coefNew = a_coefNew;
	  m_timeOld = a_timeOld;
	  m_timeNew = a_timeNew;

	  CoefficientInterpolator(m_coefOld->nComp());

  }

void CoefficientInterpolatorLinearFace::interpolate(LevelData<FluxBox>& a_result,
        Real a_time)
{
	Interval interv(0,numComps()-1);

	// Can't do this with fluxboxes
	//	timeInterp(a_result, a_time, *m_coefOld, m_timeOld,
	//					*m_coefNew, m_timeNew, interv);


	Real alpha = (a_time - m_timeOld)/(m_timeNew-m_timeOld);

	if (alpha < 0 or alpha > 1)
	{
		MayDay::Error("CoefficientInterpolatorLinearFace - alpha not between 0 and 1");
	}

	for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit)
	{

//		FluxBox& oldData = (*m_coefOld)[dit];
		FluxBox& newData = (*m_coefNew)[dit];

		for (int dir=0; dir<SpaceDim; dir++)
		{

			FArrayBox newDataTemp(newData[dir].box(), numComps());
			newDataTemp.copy(newData[dir]);
			newDataTemp.mult(1-alpha);

			a_result[dit][dir].copy((*m_coefOld)[dit][dir]);
			a_result[dit][dir].mult(alpha);
			a_result[dit][dir].plus(newData[dir]);


		}
	}
}

void CoefficientInterpolatorLinearFace::interpolate(LevelData<FluxBox>& a_result,
	                           const LevelData<FluxBox>& a_solution,
	                           Real a_time)
{
	interpolate(a_result, a_time);
}

bool CoefficientInterpolatorLinearFace::dependsUponSolution() const
{
	return false;
}

void CoefficientInterpolatorLinearFace::interpolatePrime(LevelData<FluxBox>& a_prime, const LevelData<FluxBox>& a_solution, Real a_time)
{
MayDay::Error("Not implemented");
}

void CoefficientInterpolatorLinearFace::solve(LevelData<FluxBox>&a_phi, const LevelData<FluxBox>& a_f,
		Real a_time, const LevelData<FluxBox>& a_phi0, Real a_tolerance)
{
	MayDay::Error("Not implemented");
}
