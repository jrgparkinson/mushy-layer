/*
 * Diagnostics.h
 *
 *  Created on: 30 May 2017
 *      Author: parkinsonjl
 */

#ifndef SRC_DIAGNOSTICS_H_
#define SRC_DIAGNOSTICS_H_

//#include <cmath>
#include "FArrayBox.H"
#include "Vector.H"
#include "Misc.H"
#include "parstream.H"
#include <iomanip>
#include <fstream>
#include "mushyLayerOpt.h"

#include "NamespaceHeader.H"


/// Different diagnostics we consider
enum DiagnosticNames{
  diag_time,
  diag_dt,

  diag_averageVerticalSaltFlux,
  diag_L2FsVertDiffusion,
  diag_L2FsVertFluid,
  diag_L2FsVertFrame,
  diag_L1FsVertDiffusion,
  diag_L1FsVertFluid,
  diag_L1FsVertFrame,
  diag_L0FsVertDiffusion,
  diag_L0FsVertFluid,
  diag_L0FsVertFrame,
  diag_soluteFluxTop,
  diag_soluteFluxBottom,
  diag_soluteFluxSponge,
  diag_heatFluxBottom,
  diag_heatFluxTop,
  diag_dUdt,
  diag_dSdt,
  diag_dTdt,
  diag_chimneySpacing,
  diag_chimneyWidth,
  diag_HorizAvSalinity0,
  diag_HorizAvSalinity20,
  diag_HorizAvSalinity40,
  diag_HorizAvSalinity60,
  diag_avSalinity,
  diag_mushDepth,
  diag_Nu,
  diag_NuRight,
  diag_NuLeft,
  diag_NuMiddle,
  diag_maxVhalf,
  diag_maxUhalf,
  diag_heatFluxAbsMismatch,
  diag_saltFluxAbsMismatch,
  diag_heatFluxRelMismatch,
  diag_saltFluxRelMismatch,

  diag_maxLambda,
  diag_sumLambda,
  diag_maxVel,

  // Make sure this comes last in this list!
  numDiagnostics
};

/// Class to contain diagnostics
/**
 * This class manages various diagnostics that we want to track during simulations
 */
class Diagnostics
{
public:
  /// Default constructore
  Diagnostics ();

  /// Define object
  void define (Real a_movingAverageTimescale, int a_verbosity, Real a_convCrit);

  /// Destructor
  virtual  ~Diagnostics ();

  /// Add a diagnostic
  void addDiagnostic(DiagnosticNames a_diagnostic, Real a_time, Real value);

  /// Get the value of a diagnostic, \f$ \alpha \f$
  Real getDiagnostic(DiagnosticNames a_diagnostic, Real a_time, int timestepOffset = 0);

  /// Get the moving average of a diagnostic
  Real getMovingAverage(DiagnosticNames a_diagnostic, Real a_endTime, Real a_timeSpan);

  /// Get \f$ \frac{d \alpha}{d t} \f$
  Real getRateOfChange(DiagnosticNames a_diagnostic, Real a_endTime, Real a_dt);

  /// Calculate \f$ \frac{d^2 \alpha}{dt^2} \f$
  Real getSecondRateOfChange(DiagnosticNames a_diagnostic, Real a_endTime, Real a_dt);

  /// Determine if the moving average has reached steady state
  bool movingAverageHasConverged(DiagnosticNames a_diagnostic, Real m_time, Real a_dt);

  /// Print header of all diagnostic names
  void printHeader();

  /// Print header to specified file
  void printHeader(std::ofstream& a_file);

  /// Print diagnostics at specified time
  void printDiagnostics(Real a_time);

  /// Print diagnostics at given time to a certain file
  void printDiagnostics(Real a_time, std::ofstream& a_file);

  /// Returns whether or not the specified diagnostic is one that's in one list of diagnostics to print
  bool diagnosticIsIncluded(const DiagnosticNames a_diag);



  /// Specify which diagnostics we should print out
  void setPrintDiags(Vector<DiagnosticNames> a_diagsToPrint);

private:

  /// Times at which diagnostics have been calculated
  Vector<Real> m_times;

  /// Diagnostics at each time
  Vector<Vector<Real>* > m_diagnostics;

  /// Names of diagnostics
  Vector<string> m_diagnosticNames;

  /// Diagnostics to print out
  Vector<DiagnosticNames> m_diagsToPrint;

  /// Timescale for computing moving averages
  Real movingAverageTimescale;

  /// Criteria for determining convergence
  Real m_convergenceCriteria;

  /// Verbosity
  int m_verbosity;

  /// Is object defined properly
  bool m_defined;

  /// File to which all diagnostics are written
  std::ofstream m_diagnosticsFile,

  /// File to which the latest the diagnostics are written
  m_diagnosticsFileLatest;

  /// Get the index for a certain timestep
  int getIndex(Real a_time);



};

#include "NamespaceFooter.H"

#endif /* SRC_DIAGNOSTICS_H_ */
