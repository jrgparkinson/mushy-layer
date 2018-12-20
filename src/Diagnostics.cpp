/*
 * Diagnostics.cpp
 *
 *  Created on: 30 May 2017
 *      Author: parkinsonjl
 */

#include "Diagnostics.h"
#include <iostream>
#include "NamespaceHeader.H"


Diagnostics::Diagnostics ()
{
  m_diagnostics.resize(m_numDiagnostics, NULL);
  m_diagnosticNames.resize(m_numDiagnostics);

  for (int i = 0; i < m_numDiagnostics; i++)
  {
    m_diagnostics[i] = new Vector<Real>;
    m_diagnosticNames[i] = "";
  }

  m_diagnosticNames[m_dt] = "dt";
  m_diagnosticNames[m_Nu] = "Nusselt";
  m_diagnosticNames[m_NuLeft] = "NusseltLeft";
  m_diagnosticNames[m_NuRight] = "NusseltRight";
  m_diagnosticNames[m_maxVhalf] = "MaxVHalf";
  m_diagnosticNames[m_maxUhalf] = "MaxUHalf";



  m_diagnosticNames[m_soluteFluxTop] = "Fs_top";
  m_diagnosticNames[m_soluteFluxBottom] = "Fs_bottom";
  m_diagnosticNames[m_dUdt] = "dUdt";
  m_diagnosticNames[m_dSdt] = "dSdt";
  m_diagnosticNames[m_dTdt] = "dTdt";
  m_diagnosticNames[m_time] = "time";
  m_diagnosticNames[m_HorizAvSalinity0] = "horizontallyAveragedSalinity_0";
  m_diagnosticNames[m_HorizAvSalinity20] = "horizontallyAveragedSalinity_20";
  m_diagnosticNames[m_HorizAvSalinity40] = "horizontallyAveragedSalinity_40";
  m_diagnosticNames[m_HorizAvSalinity60] = "horizontallyAveragedSalinity_60";
  m_diagnosticNames[m_avSalinity] = "averageLiquidSalinity";
  m_diagnosticNames[m_mushDepth] = "mushyLayerDepth";

  m_diagnosticNames[m_soluteFluxSponge] = "Fs_sponge";
  m_diagnosticNames[m_heatFluxBottom] = "Fh_bottom";
  m_diagnosticNames[m_heatFluxTop] = "Fh_top";
  m_diagnosticNames[m_chimneySpacing] = "chimneySpacing";
  m_diagnosticNames[m_chimneyWidth] = "chimneyWidth";

  m_diagnosticNames[m_heatFluxAbsMismatch] = "Fh_abs_mismatch";
  m_diagnosticNames[m_saltFluxAbsMismatch] = "Fs_abs_mismatch";
  m_diagnosticNames[m_heatFluxRelMismatch] = "Fh_rel_mismatch";
  m_diagnosticNames[m_saltFluxRelMismatch] = "Fs_rel_mismatch";

  m_diagnosticNames[m_averageVerticalSaltFlux] = "Fs_vertical_av";

  m_diagnosticNames[m_L2FsVertDiffusion] = "L2FsVertDiffusion";
  m_diagnosticNames[m_L2FsVertFluid] = "L2FsVertFluid";
  m_diagnosticNames[m_L2FsVertFrame] = "L2FsVertFrame";

  m_diagnosticNames[m_L1FsVertDiffusion] = "L1FsVertDiffusion";
  m_diagnosticNames[m_L1FsVertFluid] = "L1FsVertFluid";
  m_diagnosticNames[m_L1FsVertFrame] = "L1FsVertFrame";

  m_diagnosticNames[m_L0FsVertDiffusion] = "L0FsVertDiffusion";
  m_diagnosticNames[m_L0FsVertFluid] = "L0FsVertFluid";
  m_diagnosticNames[m_L0FsVertFrame] = "L0FsVertFrame";

  m_diagnosticNames[m_maxLambda] = "LambdaMax";
  m_diagnosticNames[m_sumLambda] = "LambdaSum";
  m_diagnosticNames[m_maxVel] = "maxVel";

  // Sanity check to ensure we have given all the diagnostics names
  for (int i = 0; i < m_numDiagnostics; i++)
  {
    if (m_diagnosticNames[i].empty())
    {
      MayDay::Error("Diagnostics::Diagnostics() - a diagnostic doesn't have a name specified.");
    }

    // Default: print all diagnostics
    m_diagsToPrint.push_back(i);
  }

  movingAverageTimescale = 0;
  m_verbosity = 0;
  m_convergenceCriteria = 1e-4;


  m_defined = false;

}

void Diagnostics::define(Real a_movingAverageTimescale, int a_verbosity, Real a_convCrit)
{
  movingAverageTimescale = a_movingAverageTimescale;
  m_verbosity = a_verbosity;
  m_convergenceCriteria = a_convCrit;

  m_diagnosticsFile.open("diagnostics.out", std::ios_base::app);


  if (m_verbosity > 2)
  {
    pout() << "Diagnostics::define with timescale = " << a_movingAverageTimescale << std::endl;
  }

  m_defined = true;
}

Diagnostics::~Diagnostics ()
{
  // trying to solve memory leak
  for (int i = 0; i < m_numDiagnostics; i++)
    {
      delete m_diagnostics[i];
      m_diagnostics[i] = NULL;
    }
}

void Diagnostics::addDiagnostic(int a_diagnostic, Real a_time, Real value)
{
  int index = getIndex(a_time);

  // Add this time if we haven't already
  if (index == -1)
  {
    m_times.push_back(a_time);
    index = m_times.size() - 1;
    (*m_diagnostics[m_time]).push_back(a_time);
  }

  // Check if we've already added a value for this diagnostic at this time
  if (index < m_diagnostics[a_diagnostic]->size())
  {
    (*m_diagnostics[a_diagnostic])[index] = value;

  }
  else
  {
    // Just in case we've got out of sync, add NaN entries to get vectors up to size
    while (index > m_diagnostics[a_diagnostic]->size())
    {
      Real NaN;
      m_diagnostics[a_diagnostic]->push_back(NaN);
    }

    m_diagnostics[a_diagnostic]->push_back(value);

  }

}

Real Diagnostics::getDiagnostic(int a_diagnostic, Real a_time, int timestepOffset)
{
  Real val = 1.0e200;

  int index = getIndex(a_time);

  index = index + timestepOffset;

  if (index > -1)
  {
    if (m_diagnostics[a_diagnostic]->size() <= index)
    {
      val = std::nan("1");
    }
    else
    {
      val = (*m_diagnostics[a_diagnostic])[index];
    }
  }

  return val;
}

Real Diagnostics::getMovingAverage(int a_diagnostic, Real a_endTime, Real a_timeSpan)
{
  Real movingAverage = 0;

  int endIndex = getIndex(a_endTime);
  int Ntimesteps = 0;

  for (int k = endIndex; k >= 0; k--)
  {
    if ((a_endTime - m_times[k]) < a_timeSpan)
    {
      movingAverage += (*m_diagnostics[a_diagnostic])[k];
      Ntimesteps++;
    }
  }

  movingAverage = movingAverage/Ntimesteps;

  return movingAverage;
}

Real Diagnostics::getRateOfChange(int a_diagnostic, Real a_endTime, Real a_dt)
{
  Real rateOfChange = 1.0e200;

  int index = getIndex(a_endTime);

  Vector<Real> thisDiag = *m_diagnostics[a_diagnostic];

  if (thisDiag.size() > index && index > 0)
  {
    rateOfChange = (thisDiag[index] - thisDiag[index-1])/a_dt;
  }

  return rateOfChange;
}

Real Diagnostics::getSecondRateOfChange(int a_diagnostic, Real a_endTime, Real a_dt)
{
  Real rateOfChange = 1.0e200;

  int index = getIndex(a_endTime);

  Vector<Real> thisDiag = *m_diagnostics[a_diagnostic];

  if (thisDiag.size() > index && index > 1)
  {
    rateOfChange = (thisDiag[index] + thisDiag[index-2] - thisDiag[index-1])/(a_dt*a_dt);
  }

  return rateOfChange;
}

bool Diagnostics::movingAverageHasConverged(int a_diagnostic, Real a_time, Real a_dt)
{
  //  Real d2dt2 = a_dt*(getSecondRateOfChange(a_diagnostic, a_time, a_dt));
  //  Real d_dt = (getRateOfChange(a_diagnostic, a_time, a_dt));
  //  Real delta = a_dt*d_dt;
  //  Real prevDiag  = getDiagnostic(a_diagnostic, a_time, -1);

  //  if (a_diagnostic == m_Nu
  //      && (prevDiag < 1.01 && prevDiag > 0.99) // Don't think we've converged if the nusselt number is still 1
  //      && std::abs(delta) < 1e-5
  //      && std::abs(d2dt2) < 1e-3
  //  )
  //  {
  //    return  true;
  //  }

  // Consider the moving average
  if( (a_time - m_times[0]) > 2*movingAverageTimescale)
  {
    Real movingAverageDiff = std::abs(getMovingAverage(a_diagnostic, a_time, movingAverageTimescale) - getMovingAverage(a_diagnostic, a_time, 2*movingAverageTimescale));

    if (m_verbosity > 0)
    {
      pout() << "Difference in moving average (" << m_diagnosticNames[a_diagnostic] << ") = " << movingAverageDiff << std::endl;
    }

    if ( movingAverageDiff < m_convergenceCriteria )
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  else
  {
    // We haven't got enough data to decide if we've converged
    return false;
  }
}

void Diagnostics::printHeader()
{
  printHeader(m_diagnosticsFile);
}

void Diagnostics::printHeader(std::ofstream& a_file)
{
  for (int i = 0; i < m_diagsToPrint.size(); i++)
  {
    int diag_i = m_diagsToPrint[i];

    a_file << m_diagnosticNames[diag_i];

    // Add a comma to separate entries unless it's the final entry
    if (i < m_diagsToPrint.size()-1)
    {
      a_file << ",";
    }
  }

  a_file << endl;

}

void Diagnostics::printDiagnostics(Real a_time)
{
  CH_TIME("Diagnostics::printDiagnostics");

  printDiagnostics(a_time, m_diagnosticsFile);

  // Open latest diags file and delete contents, then rewrite
  m_diagnosticsFileLatest.open("diagnosticsLatest.out", std::ofstream::out | std::ofstream::trunc);
  printHeader(m_diagnosticsFileLatest);
  printDiagnostics(a_time, m_diagnosticsFileLatest);
  m_diagnosticsFileLatest.close();

}
void Diagnostics::printDiagnostics(Real a_time, std::ofstream& a_file)
{

  for (int i = 0; i < m_diagsToPrint.size(); i++)
  {
    int diag_i = m_diagsToPrint[i];

    Real diag = getDiagnostic(diag_i, a_time);

    a_file << setprecision(10) << diag;

    // Add a comma to separate entries unless it's the final entry
    if (i < m_diagsToPrint.size()-1)
    {
      a_file << ",";
    }
  }

  a_file << endl;
}


int Diagnostics::getIndex(Real a_time)
{
  Real index = -1;

  for (int i = 0; i < m_times.size(); i++)
  {
    if (m_times[i] == a_time)
    {
      return i;
    }
  }

  return index;
}

void Diagnostics::setPrintDiags(Vector<int> a_diagsToPrint)
{
  m_diagsToPrint = a_diagsToPrint;
}

#include "NamespaceFooter.H"
