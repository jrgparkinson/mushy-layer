#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTIBC_H_
#define _ADVECTIBC_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "PhysIBC.H"

#include "UsingNamespace.H"

/// IBC for simple advection
/**
   Parameters:

 */
class AdvectIBC:public PhysIBC
{
public:
  /// Null constructor
  /**
   */
  AdvectIBC();

  AdvectIBC(int a_numComp);

  /// Destructor
  /**
   */
  ~AdvectIBC()
  {
  }

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false
     (i.e., its define() must be called before it is used) and
     m_isFortranCommonSet set to value of m_isFortranCommonset in the
     current (factory) object.
   */
  PhysIBC *new_physIBC();

  /// set default boundary values
  void setDefaultValues(int a_numComp = 1);

  /// Set boundary fluxes
  /**
   */
  void primBC(FArrayBox&            a_WGdnv,
              const FArrayBox&      a_Wextrap,
              const FArrayBox&      a_W,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const Real&           a_time);

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
   */
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  /// Set up initial conditions
  /**
   */
  void initialize(LevelData<FArrayBox>& a_U);

  /// set velocity
  void advectionVel(const RealVect& a_advVel);

  /// advection velocity. No idea what this is here for.
  const RealVect& advectionVel() const;

  /// set problem type
  void probType(const int a_probtype);

  /// set boundary value (default is 0)
  void setBoundaryValue(Real a_bcVal, int a_bcType, int a_dir,
                        Side::LoHiSide a_hiLo, int a_comp=0);

  /// Set boundary values (N dimensional version)
 void  setBoundaryValues(RealVect bcVals, IntVect bcType,
                              Side::LoHiSide a_hiLo, int a_comp = 0);

  /// Set plume details
  void setPlume(Vector<Real> a_plumeVals, Vector<Real> plumeBounds);

  /// Set the type of condition on each boundary
  void setBCType(int a_bcType, int a_dir,
                              Side::LoHiSide a_hiLo);

  /// access boundary value
  Real getBoundaryValue(int a_dir, Side::LoHiSide a_hiLo, int a_comp=0) const;

  /// Set pointer to advection velocity
  void setAdvVel(LevelData<FluxBox>* a_advVel);

  /// Set pointer to advection velocity for the box we're integrating over
  void setAdvVel(FluxBox* a_advVel);

  /// set slope value (default is 0)
  void setSlopeValue(Real a_slopeVal, int a_dir,
                     Side::LoHiSide a_hiLo);

  /// access slope value
  Real getSlopeValue(int a_dir, Side::LoHiSide a_hiLo) const;

  /// Apply artifical viscosity BCs
  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time);

  /// accessor
  int probType() const;

  /// Possible BC types
  enum bcType{
     m_dirichlet,
     m_extrap,
     m_inflowOutflow,
     m_plumeInflow,
     m_neumann,

     m_numBCTypes
   };

protected:


  /// Don't know what this is for
  RealVect m_velocity;

  /// boundary value at the plume
  Real m_bcValPlume;
  Real m_inflowVal;
  Vector<Real> m_plumeBounds;

  /// advection velocity field for the whole level
  LevelData<FluxBox>* m_advVel;

  /// advection velocity for the specific box we're integrating over
  FluxBox* m_advVelBox;

  /// Boundary value
  Vector<Real> m_bcVal[SpaceDim][2];

  /// Boundary slope value
  Real m_slopeVal[SpaceDim][2];

  Vector<Real> m_plumeVals;

  /// Boundary condition type
  Vector<int> m_bcType[SpaceDim][2];

  /// Whether or not BCs have been set
  bool m_isBCvalSet;

  /// Whether or not boundary slopes have been set
  bool m_isSlopeValSet;

  /// Whether or not BC types have been set
  bool m_isBCtypeSet;

  /// number of components
  int m_numComps;



};

#endif
