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

#include "NamespaceHeader.H"

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
  void addDiagnostic(int a_diagnostic, Real a_time, Real value);

  /// Get the value of a diagnostic, \f$ \alpha \f$
  Real getDiagnostic(int a_diagnostic, Real a_time, int timestepOffset = 0);

  /// Get the moving average of a diagnostic
  Real getMovingAverage(int a_diagnostic, Real a_endTime, Real a_timeSpan);

  /// Get \f$ \frac{d \alpha}{d t} \f$
  Real getRateOfChange(int a_diagnostic, Real a_endTime, Real a_dt);

  /// Calculate \f$ \frac{d^2 \alpha}{dt^2} \f$
  Real getSecondRateOfChange(int a_diagnostic, Real a_endTime, Real a_dt);

  /// Determine if the moving average has reached steady state
  bool movingAverageHasConverged(int a_diagnostic, Real m_time, Real a_dt);

  /// Print header of all diagnostic names
  void printHeader();

  /// Print header to specified file
  void printHeader(std::ofstream& a_file);

  /// Print diagnostics at specified time
  void printDiagnostics(Real a_time);

  /// Print diagnostics at given time to a certain file
  void printDiagnostics(Real a_time, std::ofstream& a_file);

  /// Different diagnostics we consider
  enum diagnosticNames{
    m_time,
    m_dt,

    m_averageVerticalSaltFlux,
    m_L2FsVertDiffusion,
    m_L2FsVertFluid,
    m_L2FsVertFrame,
    m_L1FsVertDiffusion,
    m_L1FsVertFluid,
    m_L1FsVertFrame,
    m_L0FsVertDiffusion,
    m_L0FsVertFluid,
    m_L0FsVertFrame,
    m_soluteFluxTop,
    m_soluteFluxBottom,
    m_soluteFluxSponge,
    m_heatFluxBottom,
    m_heatFluxTop,
    m_dUdt,
    m_dSdt,
    m_dTdt,
    m_chimneySpacing,
    m_chimneyWidth,
    m_HorizAvSalinity0,
    m_HorizAvSalinity20,
    m_HorizAvSalinity40,
    m_HorizAvSalinity60,
    m_avSalinity,
    m_mushDepth,
    m_Nu,
    m_heatFluxAbsMismatch,
    m_saltFluxAbsMismatch,
    m_heatFluxRelMismatch,
    m_saltFluxRelMismatch,

    m_maxLambda,
    m_sumLambda,
    m_maxVel,


    // Make sure this comes last in this list!
    m_numDiagnostics
  };

private:

  /// Times at which diagnostics have been calculated
  Vector<Real> m_times;

  /// Diagnostics at each time
  Vector<Vector<Real>* > m_diagnostics;

  /// Names of diagnostics
  Vector<string> m_diagnosticNames;

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

  int getIndex(Real a_time);

};

#include "NamespaceFooter.H"

#endif /* SRC_DIAGNOSTICS_H_ */
