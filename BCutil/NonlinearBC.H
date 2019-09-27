/*
 * NonlinearBC.H
 *
 *  Created on: 14 Sep 2019
 *      Author: parkinsonjl
 */

#ifndef BCUTIL_NONLINEARBC_H_
#define BCUTIL_NONLINEARBC_H_


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
 *  \f$ r =  \left[F_0 - b( T_b - T_{ref}) \right] - \frac{T_g - T_i)}{\Delta x} \f$
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
//    m_cooling = a_cooling;
  }

  virtual void computeResidual(Real& a_residual, const Real a_ghost_enthalpy, const Real a_ghost_bulk_concentration, const Real a_interior_temperature)
  {
    Real a_ghost_T = m_params.computeTemperature(a_ghost_enthalpy, a_ghost_bulk_concentration);

//    Real Tflux = (a_ghost_T - a_interior_temperature)/m_dx;
//    Real Tboundary = (a_ghost_T + a_interior_temperature)/2;

    Real predicted_ghost = (a_interior_temperature*m_a/m_dx + m_F + m_b*a_interior_temperature/2 - m_T_ref*m_b) / (m_a/m_dx - m_b/2);

//    a_residual = (m_F - m_b*(Tboundary - m_T_ref)) - m_a*Tflux;

    // weighting factor
//    a_residual = a_residual / (m_a/m_dx - m_b/2);

//    if (m_cooling)
//    {
//      predicted_ghost = min(predicted_ghost, m_T_ref);
//    }

    a_residual = ( predicted_ghost - a_ghost_T);

    // Some basic error checking
    // should never have a residual this big?
    if (abs(a_residual) > 10.0)
    {
      MayDay::Warning("NonlinearTemperatureBCRobin::computeResidual - BC residual too large");
      a_residual = 0.0;
    }
  }


protected:
  Real m_a, m_b, m_F, m_T_ref;
//  bool m_cooling; // Is this a cooling boundary?
};

/// Base class for solving for nonlinear ghost cell values
/**
 * Given a NonlinearTemperatureBC object which computes the residual, find the root
 */
class NonlinearBCSolver
{

public:

  /// Destructor
  virtual ~NonlinearBCSolver()
  {

  }

  /// Default constructor
  NonlinearBCSolver()
  {
    m_bc = NULL;
    m_max_resid = 1e-5;
    m_max_iter = 3;
    m_relax_coeff = 0.5;
  }

  NonlinearBCSolver(NonlinearTemperatureBC* a_bc,
                    const int a_max_iter,
                    const Real a_max_residual,
                    const Real a_relax_coeff = 0.5)
  {
    m_bc = a_bc;
    m_max_resid = a_max_residual;
    m_max_iter = a_max_iter;
    m_relax_coeff = a_relax_coeff;
  }

  /// Method for solving the system - must be implemented by inherited classes
  virtual void solve(Real& a_ghost_enthalpy,
                     const Real a_interior_temperature,
                     const Real a_ghost_bulk_conc)
  {

  }

  NonlinearTemperatureBC* m_bc;
  int m_max_iter;
  Real m_max_resid;
  Real m_relax_coeff;

};

/// Picard solver
class NonlinearBCSolverPicard: public NonlinearBCSolver
{
public:

  /// Destructor
    virtual ~NonlinearBCSolverPicard()
    {

    }

    /// Full constructor
    NonlinearBCSolverPicard(NonlinearTemperatureBC* a_bc,
                            const int a_max_iter,
                            const Real a_max_residual,
                            const Real a_relax_coeff) : NonlinearBCSolver(a_bc, a_max_iter, a_max_residual, a_relax_coeff)
    {
    }

    /// Method for solving the system - must be implemented by inherited classes
    virtual void solve(Real& ghost_enthalpy,
                       const Real a_interior_temperature,
                       const Real a_ghost_bulk_conc)
    {
      CH_TIME("NonlinearBCSolverPicard::solve");

      Real resid = 0.0;

//      Real relax_coeff = 0.5;

      m_bc->computeResidual(resid,
                            ghost_enthalpy,  // enthalpy in ghost cell
                            a_ghost_bulk_conc,  // bulk conc in ghost cell
                            a_interior_temperature // temperature at interior grid cell
      );

      // Compute new estimate for the enthalpy in the ghost cell, by incrementing by the residual
      ghost_enthalpy = ghost_enthalpy + m_relax_coeff*resid;

      // 10^-2 is some arbitrary criteria
      int num_iter = 0;
      while (abs(resid) > m_max_resid && num_iter < m_max_iter)
      {

        // Recompute residual
        m_bc->computeResidual(resid,
                              ghost_enthalpy,  // enthalpy in ghost cell
                              a_ghost_bulk_conc,  // bulk conc in ghost cell
                              a_interior_temperature // temperature at interior grid cell
        );



        // Compute new estimate for the enthalpy in the ghost cell, by incrementing with the residual
        ghost_enthalpy = ghost_enthalpy + m_relax_coeff*resid;

        num_iter++;

      }


    }

};


/// Newton solver
class NonlinearBCSolverNewton: public NonlinearBCSolver
{
public:

  /// Destructor
    virtual ~NonlinearBCSolverNewton()
    {

    }

    /// Full constructor
    NonlinearBCSolverNewton(NonlinearTemperatureBC* a_bc,
                            const int a_max_iter,
                            const Real a_max_residual) : NonlinearBCSolver(a_bc, a_max_iter, a_max_residual)
    {
    }

    /// Method for solving the system - must be implemented by inherited classes
    virtual void solve(Real& ghost_enthalpy,
                       const Real a_interior_temperature,
                       const Real a_ghost_bulk_conc)
    {
      CH_TIME("NonlinearBCSolverNewton::solve");

      // This is somewhat arbitrary
      Real dH = 0.01;

      Real res_old = m_max_resid*10, res_new = m_max_resid*10;
      int num_iter = 0;

      while (abs(res_old) > m_max_resid && num_iter < m_max_iter)
      {

        // Need residual at H and H+dH
        m_bc->computeResidual(res_old,
                              ghost_enthalpy,  // enthalpy in ghost cell
                              a_ghost_bulk_conc,  // bulk conc in ghost cell
                              a_interior_temperature // temperature at interior grid cell
        );

        m_bc->computeResidual(res_new,
                              ghost_enthalpy+dH,  // enthalpy in ghost cell
                              a_ghost_bulk_conc,  // bulk conc in ghost cell
                              a_interior_temperature // temperature at interior grid cell
        );

        // Near the eutectic, can have the same temperature for different enthalpies and therefore
        // the same residual. Want to avoid dividing by 0, so just do picard iterations in this case.
        if (abs(res_new - res_old) < m_max_resid)
        {
          ghost_enthalpy = ghost_enthalpy + res_old;
          // Try increasing dH to stop this happening next time
          dH = 2*dH;
        }
        else
        {
          ghost_enthalpy = ghost_enthalpy - res_old * (dH)/(res_new - res_old);
        }

        num_iter++;
      }

    }



};



#endif /* BCUTIL_NONLINEARBC_H_ */
