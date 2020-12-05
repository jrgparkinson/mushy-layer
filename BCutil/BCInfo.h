/*
 * BCInfo.h
 *
 *  Created on: 21 Apr 2019
 *      Author: parkinsonjl
 */

#ifndef BCUTIL_BCINFO_H_
#define BCUTIL_BCINFO_H_

#include "LoHiSide.H"
#include "Logging.H"
#include "ParmParse.H"
#include "UsingNamespace.H"

/// General class to contain some boundary condition information
/**
 * Hold information for some particular variable for each direction/side/component
 * Provides methods for easy access to each component
 *
 * Loads values from the inputs file of the expected format:
 *
 * bc.[name](Lo/Hi) = x-dir y-dir (z-dir, if applicable)
 *
 * e.g.
 *
 * bc.fixedFluxLo = 0.0 -1.0
 */
class BCInfo

{
public:

  // Default constructor
  BCInfo()
{
    m_name = "";
    m_required = false;
}

  /// Destructor
  virtual ~ BCInfo()
  {

  }

  /// Standard constructor
  BCInfo(const string a_name,
         const bool a_required = false,
         const Real a_default_value = 0.0)
  {
    m_name = a_name;
    m_required = a_required;


    // Setup vectors
    // Resize for two sides
    m_val.resize(2);

    Vector<string> side_labels = Vector<string>(2, string(""));
    side_labels[0] = string("Lo");
    side_labels[1] = string("Hi");

    ParmParse pp("bc");
    std::vector<float>  temp = std::vector<float>();

    // Loop over sides and load values
    for (int iside = 0; iside <= 1; iside++)
    {
      m_val[iside].resize(SpaceDim);

      string bc_name = m_name + side_labels[iside];

      if (pp.contains(bc_name.c_str()))
      {
        pp.getarr(bc_name.c_str(), temp, 0, SpaceDim);

        for (int idir=0; idir<SpaceDim; idir++)
        {
          m_val[iside][idir] = temp[idir];
        }
      }
      else
      {
        if (m_required)
        {
          LOG("Can't find BC: " << a_name);
          MayDay::Error("Couldn't find BC");
        }
        else
        {
          // Initialise to default value
          for (int idir=0; idir<SpaceDim; idir++)
          {
            m_val[iside][idir] = a_default_value;
          }
        }
      }


    } // finish looping over sides

  } // finish constructor


  Real getBC(int a_dir, Side::LoHiSide a_side)
  {
    CH_TIME("BCInfo::getBC");

    int side = int(a_side);
    return getBC(a_dir, side);
  }

  Real getBC(int a_dir, int a_side)
  {
    CH_TIME("BCInfo::getBC");

    CH_assert(a_side == 0 || a_side == 1);
    return m_val[a_side][a_dir];
  }

protected:

  string m_name;
  Vector<Vector<Real>> m_val;
  bool m_required;
};

#endif /* BCUTIL_BCINFO_H_ */
