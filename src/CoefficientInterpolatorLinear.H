#ifndef _COEFFICIENTINTERPOLATORLINEAR_H__
#define _COEFFICIENTINTERPOLATORLINEAR_H__

#include "CoefficientInterpolator.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

/// Linear coefficient interpolator
/**
 * Class to get time depdendent coefficients by interpolating linearly in time.
 */
class CoefficientInterpolatorLinear : public CoefficientInterpolator<LevelData<FArrayBox>, LevelData<FArrayBox> >
{
public:

//	virtual ~CoefficientInterpolator();

	/// default constructor
	CoefficientInterpolatorLinear(int a_numComps=1) : CoefficientInterpolator(a_numComps)
{
	  m_timeNew = -1;
	  m_timeOld = -1;
}
//	CoefficientInterpolatorLinear(LevelData<FArrayBox>* a_coefOld,
//			  LevelData<FArrayBox>* a_coefNew,
//			   Real a_timeOld, Real a_timeNew)
//{
//		define(a_coefOld, a_coefNew, a_timeOld, a_timeNew);
//}

	/// Interpolate in time
	 virtual void interpolate(LevelData<FArrayBox>& a_result,
	                           Real a_time);

	  //! Interpolate the coefficient at the given time, placing the result in the given
	  //! LevelData object. This method assumes that the coefficient depends upon the
	  //! solution to the partial differential equation in question, so the solution
	  //! is passed into it as an argument.
	  //! \param a_result The LevelData object that will store the result.
	  //! \param a_solution The solution to the partial differential equation.
	  //! \param a_time The time at which the coefficient is to be evaluated.
	  virtual void interpolate(LevelData<FArrayBox>& a_result,
	                           const LevelData<FArrayBox>& a_solution,
	                           Real a_time);

	  //! Returns true if the coefficient depends on the solution to the
	  //! partial differential equation (rendering it nonlinear), false
	  //! otherwise. By default, the coefficient is assumed not to depend
	  //! upon the solution.
	  virtual bool dependsUponSolution() const;

	  //! Computes the derivative of the coefficient with respect to the
	  //! solution at the desired time. By default, this sets \a a_deriv to 0.
	  //! \param a_prime The coefficient derivative data will be stored here.
	  //! \param a_solution The solution to the partial differential equation.
	  //! \param a_time The time at which to compute the coefficient data.
	  virtual void interpolatePrime(LevelData<FArrayBox>& a_prime,
	                                const LevelData<FArrayBox>& a_solution,
	                                Real a_time);

	  //! This virtual void method performs the iterative nonlinear solve
	  //! \f$A(\phi) \phi - f(\vec{x}) = 0\f$ for \f$\phi\f$.
	  //! \param a_phi The solution to the equation, \f$\phi\f$, will be stored here.
	  //! \param a_f The term \f$f(\vec{x})\f$ in the equation.
	  //! \param a_time The time at which the equation is solved.
	  //! \param a_phi0 The initial estimate for \f$\phi\f$.
	  //! \param a_tolerance The threshold for the error in the equation that
	  //!                    dictates when iteration should cease.
	  virtual void solve(LevelData<FArrayBox>& a_phi,
	                     const LevelData<FArrayBox>& a_f,
	                     Real a_time,
	                     const LevelData<FArrayBox>& a_phi0,
	                     Real a_tolerance);

	  /// Define coefficients
	  void define(RefCountedPtr<LevelData<FArrayBox> > a_coefOld,
			  RefCountedPtr<LevelData<FArrayBox> > a_coefNew,
	   Real a_timeOld, Real a_timeNew);

	  /// Coefficient at old time
	  RefCountedPtr<LevelData<FArrayBox> > m_coefOld;

	  /// Coefficient at new time
	  RefCountedPtr<LevelData<FArrayBox> > m_coefNew;

	  /// Old time
	  Real m_timeOld,

	  /// New time
	  m_timeNew;

private:
protected:

};

#include "NamespaceFooter.H"

#endif
