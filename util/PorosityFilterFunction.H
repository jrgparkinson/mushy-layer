/*
 * PorosityFilterFunction.H
 *
 *  Created on: 1 Mar 2020
 *      Author: parkinsonjl
 */

#ifndef UTIL_POROSITYFILTERFUNCTION_H_
#define UTIL_POROSITYFILTERFUNCTION_H_


/// Abstract porosity filter function
class PorosityFilterFunction
{
public:
  PorosityFilterFunction() {};
  virtual bool operator() (Real porosity) = 0;
};

/// Filter where lowerLimit <= porosity <= upperLimit
class PorosityFilterRange : public PorosityFilterFunction
{
public:
  Real m_lowerLimit;
  Real m_upperLimit;
  bool m_inclusive;
  PorosityFilterRange(Real lowerLimit, Real upperLimit, bool inclusive)
  : m_lowerLimit(lowerLimit), m_upperLimit(upperLimit), m_inclusive(inclusive) {};
  virtual bool operator() (Real porosity)
  {
    if (m_inclusive)
    {
      return (porosity >= m_lowerLimit && porosity <= m_upperLimit);
    }
    else
    {
      return (porosity > m_lowerLimit && porosity < m_upperLimit);
    }
  }
};

class PorosityFilterGreaterThan : public PorosityFilterFunction
{
public:
  Real m_lowerLimit;
  bool m_inclusive;
  PorosityFilterGreaterThan(Real lowerLimit, bool inclusive)
  : m_lowerLimit(lowerLimit), m_inclusive(inclusive) {};
  virtual bool operator() (Real porosity)
  {
    if (m_inclusive)
    {
      return (porosity >= m_lowerLimit);
    }
    else
    {
      return (porosity > m_lowerLimit);
    }
  }
};


#endif /* UTIL_POROSITYFILTERFUNCTION_H_ */
