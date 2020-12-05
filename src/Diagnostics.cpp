/*
 * Diagnostics.cpp
 *
 *  Created on: 30 May 2017
 *      Author: parkinsonjl
 */

#include "Diagnostics.h"
#include "Logging.H"
#include "NamespaceHeader.H"
#include <iostream>

Diagnostics::Diagnostics()
{
  m_diagnostics.resize(numDiagnostics, nullptr);
  m_diagnosticNames.resize(numDiagnostics);

  for (int i = 0; i < numDiagnostics; i++)
  {
    m_diagnostics[i] = new Vector<Real>;
    m_diagnosticNames[i] = "";
  }

  m_diagnosticNames[DiagnosticNames::diag_dt] = "dt";
  m_diagnosticNames[DiagnosticNames::diag_timestep] = "timestep";
  m_diagnosticNames[DiagnosticNames::diag_level] = "level";
  m_diagnosticNames[DiagnosticNames::diag_Nu] = "Nusselt";
  m_diagnosticNames[DiagnosticNames::diag_NuLeft] = "NusseltLeft";
  m_diagnosticNames[DiagnosticNames::diag_NuMiddle] = "NusseltMiddle";
  m_diagnosticNames[DiagnosticNames::diag_NuRight] = "NusseltRight";
  m_diagnosticNames[DiagnosticNames::diag_maxVhalf] = "MaxVHalf";
  m_diagnosticNames[DiagnosticNames::diag_maxUhalf] = "MaxUHalf";

  m_diagnosticNames[DiagnosticNames::diag_soluteFluxTop] = "Fs_top";
  m_diagnosticNames[DiagnosticNames::diag_soluteFluxBottom] = "Fs_bottom";
  m_diagnosticNames[DiagnosticNames::diag_dUdt] = "dUdt";
  m_diagnosticNames[DiagnosticNames::diag_dSdt] = "dSdt";
  m_diagnosticNames[DiagnosticNames::diag_dTdt] = "dTdt";
  m_diagnosticNames[DiagnosticNames::diag_time] = "time";
  m_diagnosticNames[DiagnosticNames::diag_HorizAvSalinity0] =
      "horizontallyAveragedSalinity_0";
  m_diagnosticNames[DiagnosticNames::diag_HorizAvSalinity20] =
      "horizontallyAveragedSalinity_20";
  m_diagnosticNames[DiagnosticNames::diag_HorizAvSalinity40] =
      "horizontallyAveragedSalinity_40";
  m_diagnosticNames[DiagnosticNames::diag_HorizAvSalinity60] =
      "horizontallyAveragedSalinity_60";
  m_diagnosticNames[DiagnosticNames::diag_avSalinity] = "averageLiquidSalinity";
  m_diagnosticNames[DiagnosticNames::diag_mushDepth] = "mushyLayerDepth";

  m_diagnosticNames[DiagnosticNames::diag_soluteFluxSponge] = "Fs_sponge";
  m_diagnosticNames[DiagnosticNames::diag_heatFluxBottom] = "Fh_bottom";
  m_diagnosticNames[DiagnosticNames::diag_heatFluxTop] = "Fh_top";
  m_diagnosticNames[DiagnosticNames::diag_chimneySpacing] = "chimneySpacing";
  m_diagnosticNames[DiagnosticNames::diag_chimneyWidth] = "chimneyWidth";

  m_diagnosticNames[DiagnosticNames::diag_heatFluxAbsMismatch] =
      "Fh_abs_mismatch";
  m_diagnosticNames[DiagnosticNames::diag_saltFluxAbsMismatch] =
      "Fs_abs_mismatch";
  m_diagnosticNames[DiagnosticNames::diag_heatFluxRelMismatch] =
      "Fh_rel_mismatch";
  m_diagnosticNames[DiagnosticNames::diag_saltFluxRelMismatch] =
      "Fs_rel_mismatch";

  m_diagnosticNames[DiagnosticNames::diag_averageVerticalSaltFlux] =
      "Fs_vertical_av";
  m_diagnosticNames[DiagnosticNames::diag_Fs10] = "Fs_vert_10pc";
  m_diagnosticNames[DiagnosticNames::diag_Fs20] = "Fs_vert_20pc";
  m_diagnosticNames[DiagnosticNames::diag_Fs30] = "Fs_vert_30pc";
  m_diagnosticNames[DiagnosticNames::diag_Fs40] = "Fs_vert_40pc";
  m_diagnosticNames[DiagnosticNames::diag_Fs50] = "Fs_vert_50pc";

  m_diagnosticNames[DiagnosticNames::diag_L2FsVertDiffusion] =
      "L2FsVertDiffusion";
  m_diagnosticNames[DiagnosticNames::diag_L2FsVertFluid] = "L2FsVertFluid";
  m_diagnosticNames[DiagnosticNames::diag_L2FsVertFrame] = "L2FsVertFrame";

  m_diagnosticNames[DiagnosticNames::diag_L1FsVertDiffusion] =
      "L1FsVertDiffusion";
  m_diagnosticNames[DiagnosticNames::diag_L1FsVertFluid] = "L1FsVertFluid";
  m_diagnosticNames[DiagnosticNames::diag_L1FsVertFrame] = "L1FsVertFrame";

  m_diagnosticNames[DiagnosticNames::diag_L0FsVertDiffusion] =
      "L0FsVertDiffusion";
  m_diagnosticNames[DiagnosticNames::diag_L0FsVertFluid] = "L0FsVertFluid";
  m_diagnosticNames[DiagnosticNames::diag_L0FsVertFrame] = "L0FsVertFrame";

  m_diagnosticNames[DiagnosticNames::diag_maxLambda] = "LambdaMax";
  m_diagnosticNames[DiagnosticNames::diag_sumLambda] = "LambdaSum";
  m_diagnosticNames[DiagnosticNames::diag_postRegridLambda] =
      "LambdaPostRegrid";
  m_diagnosticNames[DiagnosticNames::diag_maxVel] = "maxVel";
  m_diagnosticNames[DiagnosticNames::diag_maxFreestreamCorrection] =
      "maxFreestreamCorrection";

  m_diagnosticNames[DiagnosticNames::diag_mushyAverageBulkConc] =
      "mushAvBulkConc";
  m_diagnosticNames[DiagnosticNames::diag_mushyAveragePorosity] =
      "mushAvPorosity";
  m_diagnosticNames[DiagnosticNames::diag_mushyVol] = "mushVol";

  // Sanity check to ensure we have given all the diagnostics names
  for (int i = 0; i < numDiagnostics; i++)
  {
    if (m_diagnosticNames[i].empty())
    {
      MayDay::Error("Diagnostics::Diagnostics() - a diagnostic doesn't have a "
                    "name specified.");
    }

    // Default: print all diagnostics
    m_diagsToPrint.push_back(DiagnosticNames(i));
  }

  movingAverageTimescale = 0;
  m_verbosity = 0;
  m_convergenceCriteria = 1e-4;
  m_level = 0;

  m_defined = false;
}

void Diagnostics::define(Real a_movingAverageTimescale, int a_verbosity,
                         Real a_convCrit, int a_level, bool a_printAllLevels)
{
  movingAverageTimescale = a_movingAverageTimescale;
  m_verbosity = a_verbosity;
  m_convergenceCriteria = a_convCrit;
  m_level = a_level;
  m_printAllLevels = a_printAllLevels;

  if (m_printAllLevels)
  {
    m_diagnosticsFileName = "diagnostics." + std::to_string(a_level) + ".csv";
  }
  else
  {
    m_diagnosticsFileName = "diagnostics.csv";
  }

  if (m_verbosity > 2)
  {
    LOG("Diagnostics::define with timescale = " << a_movingAverageTimescale);
  }

  m_defined = true;
}

Diagnostics::~Diagnostics()
{
  // trying to solve memory leak
  for (int i = 0; i < numDiagnostics; i++)
  {
    delete m_diagnostics[i];
    m_diagnostics[i] = nullptr;
  }
}

void Diagnostics::addDiagnostic(DiagnosticNames a_diagnostic, Real a_time,
                                Real value)
{
  int index = getIndex(a_time);

  // Add this time if we haven't already
  if (index == -1)
  {
    m_times.push_back(a_time);
    index = m_times.size() - 1;
    (*m_diagnostics[diag_time]).push_back(a_time);
    (*m_diagnostics[diag_level]).push_back(m_level);
  }

  // Check if we've already added a value for this diagnostic at this time
  if (index < m_diagnostics[a_diagnostic]->size())
  {
    (*m_diagnostics[a_diagnostic])[index] = value;
  }
  else
  {
    // Just in case we've got out of sync, add NaN entries to get vectors up to
    // size
    while (index > m_diagnostics[a_diagnostic]->size())
    {
      Real NaN;
      m_diagnostics[a_diagnostic]->push_back(NaN);
    }

    m_diagnostics[a_diagnostic]->push_back(value);
  }
}

Real Diagnostics::getDiagnostic(DiagnosticNames a_diagnostic, Real a_time,
                                int timestepOffset)
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

Real Diagnostics::getMovingAverage(DiagnosticNames a_diagnostic, Real a_endTime,
                                   Real a_timeSpan)
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

  movingAverage = movingAverage / Ntimesteps;

  return movingAverage;
}

Real Diagnostics::getRateOfChange(DiagnosticNames a_diagnostic, Real a_endTime,
                                  Real a_dt)
{
  Real rateOfChange = 1.0e200;

  int index = getIndex(a_endTime);

  Vector<Real> thisDiag = *m_diagnostics[a_diagnostic];

  if (thisDiag.size() > index && index > 0)
  {
    rateOfChange = (thisDiag[index] - thisDiag[index - 1]) / a_dt;
  }

  return rateOfChange;
}

Real Diagnostics::getSecondRateOfChange(DiagnosticNames a_diagnostic,
                                        Real a_endTime, Real a_dt)
{
  Real rateOfChange = 1.0e200;

  int index = getIndex(a_endTime);

  Vector<Real> thisDiag = *m_diagnostics[a_diagnostic];

  if (thisDiag.size() > index && index > 1)
  {
    rateOfChange =
        (thisDiag[index] + thisDiag[index - 2] - thisDiag[index - 1]) /
        (a_dt * a_dt);
  }

  return rateOfChange;
}

bool Diagnostics::movingAverageHasConverged(DiagnosticNames a_diagnostic,
                                            Real a_time, Real a_dt)
{

  // Consider the moving average
  if ((a_time - m_times[0]) > 2 * movingAverageTimescale)
  {
    Real movingAverageDiff = std::abs(
        getMovingAverage(a_diagnostic, a_time, movingAverageTimescale) -
        getMovingAverage(a_diagnostic, a_time, 2 * movingAverageTimescale));

    LOG("Difference in moving average (" << m_diagnosticNames[a_diagnostic]
                                         << ") = " << movingAverageDiff);

    if (movingAverageDiff < m_convergenceCriteria)
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

bool Diagnostics::diagnosticIsIncluded(const DiagnosticNames a_diag)
{
  std::vector<DiagnosticNames> stdVect = m_diagsToPrint.stdVector();
  if (std::find(stdVect.begin(), stdVect.end(), a_diag) != stdVect.end())
  {
    return true;
  }
  else
  {
    return false;
  }
}

void Diagnostics::printHeader()
{
  if (m_printAllLevels || m_level == 0)
  {
    m_diagnosticsFile.open(m_diagnosticsFileName,
                           std::ios_base::app); // open file in append mode
    printHeader(m_diagnosticsFile);
    m_diagnosticsFile.close();
  }
}

void Diagnostics::printHeader(std::ofstream &a_file)
{
  for (int i = 0; i < m_diagsToPrint.size(); i++)
  {
    int diag_i = m_diagsToPrint[i];

    a_file << m_diagnosticNames[diag_i];

    // Add a comma to separate entries unless it's the final entry
    if (i < m_diagsToPrint.size() - 1)
    {
      a_file << ",";
    }
  }

  a_file << endl;
}

void Diagnostics::printDiagnostics(Real a_time)
{
  CH_TIME("Diagnostics::printDiagnostics");

  if (m_printAllLevels || m_level == 0)
  {

    LOG("Print diagnostics (level " << m_level << ", time " << a_time << ")");

    m_diagnosticsFile.open(m_diagnosticsFileName,
                           std::ios_base::app); // open file in append mode
    printDiagnostics(a_time, m_diagnosticsFile);
    m_diagnosticsFile.close();

    // Open latest diags file and delete contents, then rewrite
    std::string latestFilename = "diagnosticsLatest.csv";
    if (m_printAllLevels)
    {
      latestFilename = "diagnosticsLatest." + std::to_string(m_level) + ".csv";
    }
    m_diagnosticsFileLatest.open(latestFilename,
                                 std::ofstream::out | std::ofstream::trunc);
    printHeader(m_diagnosticsFileLatest);
    printDiagnostics(a_time, m_diagnosticsFileLatest);
    m_diagnosticsFileLatest.close();
  }
}
void Diagnostics::printDiagnostics(Real a_time, std::ofstream &a_file)
{
  for (int i = 0; i < m_diagsToPrint.size(); i++)
  {
    DiagnosticNames diag_i = m_diagsToPrint[i];

    Real diag = getDiagnostic(diag_i, a_time);

    a_file << setprecision(10) << diag;

    // Add a comma to separate entries unless it's the final entry
    if (i < m_diagsToPrint.size() - 1)
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

void Diagnostics::setPrintDiags(Vector<DiagnosticNames> a_diagsToPrint)
{
  m_diagsToPrint = a_diagsToPrint;
}

#include "NamespaceFooter.H"
