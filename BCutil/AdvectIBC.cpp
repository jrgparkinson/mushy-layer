#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "SolidF_F.H"
#include "AdvectIBC.H"
#include "Logging.H"

// Null constructor
AdvectIBC::AdvectIBC()
{
  setDefaultValues();
}

AdvectIBC::AdvectIBC(int a_numComp)
{
  setDefaultValues(a_numComp);
}

void
AdvectIBC::setDefaultValues(int a_numComp)
{
  m_numComps = a_numComp;

  // default values
  for (int dir=0; dir<SpaceDim; dir++)
  {
    m_bcVal[dir][0].resize(m_numComps);
    m_bcVal[dir][1].resize(m_numComps);
    m_bcType[dir][0].resize(m_numComps);
    m_bcType[dir][1].resize(m_numComps);

    for (int comp=0; comp<m_numComps; comp++)
    {
      m_bcVal[dir][0][comp] = 0.0;
      m_bcVal[dir][1][comp] = 0.0;
      m_slopeVal[dir][0] = 0.0;
      m_slopeVal[dir][1] = 0.0;

      m_bcType[dir][0][comp] = m_extrap;
      m_bcType[dir][1][comp] = m_extrap;
    }
  }
  m_isBCvalSet = false;
  m_isSlopeValSet = false;
  m_isBCtypeSet = false;
  m_advVel=nullptr;
  m_advVelBox = nullptr;
}

// Factory method - this object is its own factory:
// Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
// its define() must be called before it is used) and m_isFortranCommonSet
// set to value of m_isFortranCommonset in the current (factory) object.
PhysIBC* AdvectIBC::new_physIBC()
{
  AdvectIBC* retval = new AdvectIBC();
  retval->m_velocity = m_velocity;
  retval->m_isBCvalSet = m_isBCvalSet;
  retval->m_isBCtypeSet = m_isBCtypeSet;
  retval->m_isSlopeValSet = m_isSlopeValSet;
  retval->m_advVel = m_advVel;
  retval->m_advVelBox = m_advVelBox;
  retval->m_bcValPlume = m_bcValPlume;
  retval->m_plumeBounds = m_plumeBounds;
  retval->m_plumeVals = m_plumeVals;
  retval->m_numComps = m_numComps;

  for (int dir=0; dir<SpaceDim; dir++)
  {
    (retval->m_bcVal[dir][0]).resize(m_numComps);
    (retval->m_bcVal[dir][1]).resize(m_numComps);
    (retval->m_bcType[dir][0]).resize(m_numComps);
    (retval->m_bcType[dir][1]).resize(m_numComps);

    for (int comp=0; comp<m_numComps; comp++)
    {
      if (m_isBCvalSet)
      {
        retval->m_bcVal[dir][0][comp] = m_bcVal[dir][0][comp] ;
        retval->m_bcVal[dir][1][comp]  = m_bcVal[dir][1][comp] ;
      }

      if (m_isBCtypeSet)
      {
        retval->m_bcType[dir][0][comp] = m_bcType[dir][0][comp] ;
        retval->m_bcType[dir][1][comp]  = m_bcType[dir][1][comp] ;
      }

      if (m_isSlopeValSet)
      {
        retval->m_slopeVal[dir][0] = m_slopeVal[dir][0];
        retval->m_slopeVal[dir][1] = m_slopeVal[dir][1];
      }
    }
  }

  return static_cast<PhysIBC*>(retval);
}

// Set boundary fluxes
void AdvectIBC::primBC(FArrayBox&            a_WGdnv,
                       const FArrayBox&      a_Wextrap,
                       const FArrayBox&      a_W,
                       const int&            a_dir,
                       const Side::LoHiSide& a_side,
                       const Real&           a_time)
{
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
  {
    // This needs to be fixed.
    CH_assert(a_WGdnv.nComp() == m_numComps);

    for (int comp = 0; comp < m_numComps; comp++)
    {
      int lohisign;
      Box tmp = a_WGdnv.box() & m_domain;
      Real bcVal;
      int bcType;
      Real plumeVal = m_plumeVals[comp];

      // Determine which side and thus shifting directions
      if (a_side == Side::Lo)
      {
        lohisign = -1;
        bcVal = m_bcVal[a_dir][0][comp];
        bcType = m_bcType[a_dir][0][comp];
      }
      else
      {
        lohisign = 1;
        bcVal = m_bcVal[a_dir][1][comp];
        bcType = m_bcType[a_dir][1][comp];
      }

      tmp.shiftHalf(a_dir,lohisign);

      // Is there a domain boundary next to this grid
      if (!m_domain.contains(tmp))
      {
        tmp &= m_domain;

        Box boundaryBox;

        // Find the strip of cells next to the domain boundary
        if (a_side == Side::Lo)
        {
          boundaryBox = bdryLo(tmp,a_dir);
        }
        else
        {
          boundaryBox = bdryHi(tmp,a_dir);
        }

        // Let's get the FAB for the velocity at this point
        FluxBox* advVel = nullptr;
        if (m_advVelBox)
        {
          advVel = m_advVelBox;
        }
        else
        {
          if (m_advVel == nullptr)
          {
            LOG("AdvectIBC error - haven't specified the velocity field");
          }

          DataIterator dit = m_advVel->dataIterator();
          for(dit.reset(); dit.ok(); ++dit)
          {
            Box advVelBox = (*m_advVel)[dit].box();
            if (advVelBox.bigEnd() >= boundaryBox.bigEnd() &&
                advVelBox.smallEnd() <= boundaryBox.smallEnd())
            {
              break;
            }
          }

          // What if we don't find a data index?
          if (!dit.ok())
          {
            LOG("AdvectIBC error - velocity field doesn't include boundary");
          }

          advVel = &(*m_advVel)[dit];
        } // end search for advection velocity

        const FArrayBox& advVelDir = (*advVel)[a_dir];

        switch(bcType)
        {

          case m_dirichlet:

            // Set variable value at the boundary
            /* e.g. if we're calculating the advective flux u.s, set the value of s*/
            FORT_SOLIDBCF(CHF_FRA(a_WGdnv),
                          CHF_CONST_FRA(a_Wextrap),
                          CHF_CONST_REAL(bcVal),
                          CHF_CONST_INT(lohisign),
                          CHF_CONST_REAL(m_dx),
                          CHF_CONST_INT(a_dir),
                          CHF_BOX(boundaryBox),
                          CHF_CONST_INT(comp));

            break;

          case m_neumann:
            FORT_NEUMANN(CHF_FRA(a_WGdnv),
                         CHF_CONST_REAL(bcVal),
                         CHF_CONST_INT(lohisign),
                         CHF_CONST_REAL(m_dx),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(boundaryBox),
                         CHF_CONST_INT(comp));

            break;

          case m_inflowOutflow:
            // For inflow - flux = boundary value
            // for outflow - flux = extrapolation
            FORT_ADVINFLOWOUTFLOW(CHF_FRA(a_WGdnv),
                                  CHF_CONST_FRA(advVelDir),
                                  CHF_CONST_FRA(a_Wextrap),
                                  CHF_CONST_REAL(bcVal),
                                  CHF_CONST_INT(lohisign),
                                  CHF_CONST_REAL(m_dx),
                                  CHF_CONST_INT(a_dir),
                                  CHF_BOX(boundaryBox),
                                  CHF_CONST_INT(comp));
            break;

          case m_extrap:

            // Extrapolate from interior
            FORT_SOLIDEXTRAPBCF(CHF_FRA(a_WGdnv),
                                CHF_CONST_FRA(a_Wextrap),
                                CHF_CONST_REAL(bcVal),
                                CHF_CONST_INT(lohisign),
                                CHF_CONST_REAL(m_dx),
                                CHF_CONST_INT(a_dir),
                                CHF_BOX(boundaryBox),
                                CHF_CONST_INT(comp));
            break;

          case m_plumeInflow:
            Real plumeStart = m_plumeBounds[0];
            Real plumeEnd = m_plumeBounds[1];
            FORT_ADVINFLOWPLUME(CHF_FRA(a_WGdnv),
                                CHF_CONST_FRA(a_Wextrap),
                                CHF_CONST_REAL(bcVal),
                                CHF_CONST_REAL(plumeVal),
                                CHF_CONST_REAL(plumeStart),
                                CHF_CONST_REAL(plumeEnd),
                                CHF_CONST_INT(lohisign),
                                CHF_CONST_REAL(m_dx),
                                CHF_CONST_INT(a_dir),
                                CHF_BOX(boundaryBox),
                                CHF_CONST_INT(comp));
            break;
        } // end loop over different bc types

      } // end if domain contains box

    } // end loop over components

  } // end if not periodic
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void AdvectIBC::setBdrySlopes(FArrayBox&       a_dW,
                              const FArrayBox& a_W,
                              const int&       a_dir,
                              const Real&      a_time)
{
  // We don't know what the slopes should be, so do nothing here
}

// Set up initial conditions
void AdvectIBC::initialize(LevelData<FArrayBox>& a_U)
{
  /// shouldn't be in this function
  MayDay::Error("Shouldn't be in AdvectIBC::initialize");
}

/// set velocity
void
AdvectIBC::advectionVel(const RealVect& a_advVel)
{
  m_velocity = a_advVel;
}

/// get advection velocity
const RealVect&
AdvectIBC::advectionVel() const
{
  return m_velocity;
}
void
AdvectIBC::setBoundaryValues(RealVect a_bcVals, IntVect a_bcType,
                            Side::LoHiSide a_hiLo, int a_comp)
{
  for (int dir = 0; dir<SpaceDim; dir++)
  {
  if (a_hiLo == Side::Lo)
    {

      m_bcVal[dir][0][a_comp] = a_bcVals[dir];
      m_bcType[dir][0][a_comp] = a_bcType[dir];
    }
    else
    {
      m_bcVal[dir][1][a_comp] = a_bcVals[dir];
      m_bcType[dir][1][a_comp] = a_bcType[dir];
    }
  }

    m_isBCvalSet = true;
    m_isBCtypeSet = true;
}

void
AdvectIBC::setBoundaryValue(Real a_bcVal, int a_bcType, int a_dir,
                            Side::LoHiSide a_hiLo, int a_comp)
{
  if (a_hiLo == Side::Lo)
  {
    m_bcVal[a_dir][0][a_comp] = a_bcVal;
    m_bcType[a_dir][0][a_comp] = a_bcType;
  }
  else
  {
    m_bcVal[a_dir][1][a_comp] = a_bcVal;
    m_bcType[a_dir][1][a_comp] = a_bcType;
  }

  m_isBCvalSet = true;
  m_isBCtypeSet = true;
}

void AdvectIBC::setPlume(Vector<Real> a_plumeVals, Vector<Real> plumeBounds)
{
  m_plumeBounds = plumeBounds;
  m_plumeVals = a_plumeVals;
}

void
AdvectIBC::setBCType(int a_bcType, int a_dir,
                     Side::LoHiSide a_hiLo)
{
  if (a_hiLo == Side::Lo)
  {
    m_bcType[a_dir][0] = a_bcType;
  }
  else
  {
    m_bcType[a_dir][1] = a_bcType;
  }

  m_isBCtypeSet = true;
}

void AdvectIBC::setAdvVel(LevelData<FluxBox>* a_advVel)
{
  m_advVel = a_advVel;
}
void AdvectIBC::setAdvVel(FluxBox* a_advVel)
{
  m_advVelBox = a_advVel;
}

Real
AdvectIBC::getBoundaryValue(int a_dir, Side::LoHiSide a_hiLo, int a_comp) const
{
  Real bcval;
  if (a_hiLo == Side::Lo)
  {
    bcval = m_bcVal[a_dir][0][a_comp];
  }
  else
  {
    bcval = m_bcVal[a_dir][1][a_comp];
  }

  return bcval;
}

void
AdvectIBC::setSlopeValue(Real a_slopeVal, int a_dir,
                         Side::LoHiSide a_hiLo)
{
  // do nothing
}

Real
AdvectIBC::getSlopeValue(int a_dir, Side::LoHiSide a_hiLo) const
{
  return -1.0;
}

void
AdvectIBC::artViscBC(FArrayBox&       a_F,
                     const FArrayBox& a_U,
                     const FArrayBox& a_divVel,
                     const int&       a_dir,
                     const Real&      a_time)
{
  MayDay::Error("AdvectIBC::artViscBC - not implemented");
}
