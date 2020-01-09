#include "PhysBCUtil.H"

#include "MushyLayerUtils.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "VelIBC.H"
#include "AdvectIBC.H"
#include "BCFunctions.H"
#include "computeSum.H"
#include "NonlinearBC.H"
#include "viscousBCF_F.H"

//#include "NamespaceHeader.H"

static void getDomainFacePosition(RealVect&             a_retval,
                                  const IntVect&        a_validIV,
                                  const Real&           a_dx,
                                  const int&            a_dir,
                                  const Side::LoHiSide& a_side)
{
  Real* dataPtr = a_retval.dataPtr();

  D_TERM( dataPtr[0] = a_dx*(a_validIV[0] + 0.5);,\
  dataPtr[1] = a_dx*(a_validIV[1] + 0.5);,\
  dataPtr[2] = a_dx*(a_validIV[2] + 0.5);)

  int isign = sign(a_side);
  dataPtr[a_dir] += 0.5*Real(isign)*a_dx;
}


// ---------------------------------------------------------------
extern "C"
{

  static void ExtraBC(FArrayBox&     a_state,
                      const Box&     a_valid,
                      int            a_dir,
                      Side::LoHiSide a_side,
                      int            a_order,
                      int            a_comp = 0)
  {
    // Fortran version
    CH_assert(a_comp < a_state.nComp());
    ExtrapBC(a_state, a_valid, a_dir, a_side, a_order, a_comp);
  }


}

/// Constant value boundary condition function
/**
 * Returns a single scalar value, independent of side or position along side,
 * which can be different for different components.
 */
class ConstValueFunction: public BCValueFunction
{
public:

  /// The constant value to apply
  Real m_value;

  /// The component to apply this BC to
  int m_nComp;

  /// Default constructor
  ConstValueFunction()
  :
    m_value(0), m_nComp(0)
  {
  }

  /// Constructor
  ConstValueFunction(Real a_value,int a_nComp)
  : m_value(a_value),
    m_nComp(a_nComp)
  {
  }

  /// Apply BC
  virtual void operator()(Real*           a_pos,
                          int*            a_dir,
                          Side::LoHiSide* a_side,
                          Real*           a_value)
  {
    if (m_nComp > 0)
    {
      for (int comp = 0; comp < m_nComp; comp++)
      {
        a_value[comp] = m_value;
      }
    }
    else
    {
      MayDay::Error("undefined ConstValueFunction object");
    }
  }
};

/// Return one value for inflow, and another if not
/*
 * Assumes plume along top boundary
 */
class InflowValueFunction: public BCValueFunction
{
public:

  /// The constant value to apply
  Real m_value;
  /// The component to apply this BC to
  int m_nComp;

  /// Value during inflow
  Real m_inflowValue;

  /// Position where plume starts
  Real m_plumeStart;

  /// Position where plume ends
  Real m_plumeEnd;


  /// Default constructor
  InflowValueFunction()
  :
    m_value(0), m_nComp(0), m_inflowValue(0.0),m_plumeStart(0.0), m_plumeEnd(0.0)
  {
  }

  /// Constructor
  InflowValueFunction(Real a_value, int a_nComp,
                      Real a_inflowValue, Real a_plumeStart, Real a_plumeEnd)
  : m_value(a_value),
    m_nComp(a_nComp),
    m_inflowValue(a_inflowValue),
    m_plumeStart(a_plumeStart),
    m_plumeEnd(a_plumeEnd)
  {
  }

  /// Apply BC
  virtual void operator()(Real*           a_pos,
                          int*            a_dir,
                          Side::LoHiSide* a_side,
                          Real*           a_value)
  {
    if (m_nComp > 0)
    {
      for (int comp = 0; comp < m_nComp; comp++)
      {
        if (*a_pos > m_plumeStart && *a_pos < m_plumeEnd)
        {
          a_value[comp] = m_inflowValue;
        }
        else
        {
          a_value[comp] = m_value;
        }
      }
    }
    else
    {
      MayDay::Error("undefined ConstValueFunction object");
    }
  }
};

/// Computes values for the pressure at boundaries to enforce inflow
class PressureInflowValueFunction: public BCValueFunction
{
public:

  /// The constant value to apply
  Real m_value;
  /// The component to apply this BC to
  int m_nComp;

  /// Boundary value for inflow
  Real m_inflowValue;

  /// Position along domain edge where plume starts
  Real m_plumeStart;

  /// Position along domain edge where plume ends
  Real m_plumeEnd;

  /// Params object
  MushyLayerParams* m_params;

  /// Default constructor
  PressureInflowValueFunction()
  :
    m_value(0), m_nComp(0), m_inflowValue(0.0), m_plumeStart(0.0), m_plumeEnd(0.0),
    m_params(nullptr)
  {
  }

  /// Constructor
  PressureInflowValueFunction(Real a_value, int a_nComp,
                              Real a_inflowValue, Real a_plumeStart, Real a_plumeEnd,
                              MushyLayerParams* a_params)
  : m_value(a_value),
    m_nComp(a_nComp),
    m_inflowValue(a_inflowValue),
    m_plumeStart(a_plumeStart),
    m_plumeEnd(a_plumeEnd),
    m_params(a_params)

  {
  }

  /// Apply BC
  virtual void operator()(Real*           a_pos,
                          int*            a_dir,
                          Side::LoHiSide* a_side,
                          Real*           a_value)
  {
    if (m_nComp > 0)
    {
      for (int comp = 0; comp < m_nComp; comp++)
      {
        if (*a_pos > m_plumeStart && *a_pos < m_plumeEnd)
        {
          //          a_value[comp] = m_inflowValue;
          Real x = (*a_pos - m_plumeStart)/(m_plumeEnd-m_plumeStart);

          Real w = m_params->inflowVelocity * ( 0.5*0.5 - (x-0.5)*(x-0.5)  ) / (0.5*0.5);
          w = m_params->inflowVelocity;
          a_value[comp] = -  w + m_params->rayleighComposition * m_params->ThetaPlumeInflow;
        }
        else
        {
          //          a_value[comp] = m_value;
          //          a_value[comp] = m_params->rayleighComposition * m_params->ThetaTop;
          a_value[comp] = 0;
        }
      }
    }
    else
    {
      MayDay::Error("undefined ConstValueFunction object");
    }
  }
};


/// Abstract boundary condition for cell centres
class AbstractScalarBCFunction: public BCFunction
{
public:

  /// Is object defined?
  bool m_isDefined;
  /// Physical parameters for the problem
  MushyLayerParams m_params;
  /// Always enforce homogeneous BCs?
  bool m_homogeneous;
  /// Advection velocity - for inflow/outflow BCs
  LevelData<FluxBox>* m_advVel;
  /// Grid spacing for m_advVel
  Real m_dx;

  /// Interval to fill
  Interval m_interval;




  /// Default constructor
  AbstractScalarBCFunction()
  :
    m_isDefined(false), m_homogeneous(false), m_advVel(nullptr), m_dx(-1), m_interval(Interval(0,0))
  {
  }

  /// Full constructor
  AbstractScalarBCFunction(bool a_isDefined,
                           MushyLayerParams a_params,
                           bool a_homogeneous,
                           LevelData<FluxBox>* a_advVel,
                           Real a_dx,
                           Interval a_interval = Interval(0,0))
  :
    m_isDefined(a_isDefined),
    m_params(a_params),
    m_homogeneous(a_homogeneous),
    m_advVel(a_advVel),
    m_dx(a_dx),
    m_interval(a_interval)
  {
  }

  /// Apply BC
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous) = 0;

};


/// Abstract boundary condition for cell faces
class AbstractFaceBCFunction: public BCFunction
{
public:
  /// Always enforce homogeneous form of BCs?
  bool m_isHomogeneous;
  /// Component of velocity to apply BCs to
  int m_comp;
  /// Physical parameters for the problem
  MushyLayerParams m_params;

  /// Default constructor
  AbstractFaceBCFunction()
  :
    m_isHomogeneous(false), m_comp(-1)
  {
  }
  /// Full constructor
  AbstractFaceBCFunction(bool a_isHomogeneous,
                         int  a_comp,
                         MushyLayerParams a_params):
                           m_isHomogeneous(a_isHomogeneous),
                           m_comp(a_comp),
                           m_params(a_params)
  {
  }
  /// Apply BC
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous) = 0;

  /// Apply dirichlet BCs (edge centred)
  /// apply the same BC value at all points
  void DiriEdgeBC(FArrayBox&           a_state,
                  const Box&           a_valid,
                  Real                 a_dx,
                  bool                 a_homogeneous,
                  BCValueHolder        a_value,
                  int                  a_idir,
                  Side::LoHiSide       a_side)
  {
    // a_state is FACE-centered on face idir;
    // a_valid is CELL-centered.
    Box face(a_valid);
    face.surroundingNodes(a_idir);
    int coord = face.sideEnd(a_side)[a_idir];
    face.setRange(a_idir, coord);

    if (a_homogeneous)
    {
      a_state.setVal(0.);
    }
    else
    {
      RealVect facePos;
      int ncomp = a_state.nComp();
      Real* value = new Real[ncomp];
      RealVect junk;
      a_value(junk.dataPtr(), &a_idir, &a_side, value);
      for (int comp = 0; comp < ncomp; comp++)
      {
        a_state.setVal(value[comp], face, comp);
      }

      delete[] value;
      value = nullptr;
    }
  }

  /// Apply dirichlet BCs (edge centred), allow different value at different locations
  void DiriEdgeVariableBC(FArrayBox&           a_state,
                          const Box&           a_valid,
                          Real                 a_dx,
                          bool                 a_homogeneous,
                          BCValueHolder        a_value,
                          int                  a_idir,
                          Side::LoHiSide       a_side)
  {

    // a_state is FACE-centered on face idir;
    // a_valid is CELL-centered.
    Interval a_interval = a_state.interval();
    Box toRegion(a_valid);
    toRegion.surroundingNodes(a_idir);
    int coord = toRegion.sideEnd(a_side)[a_idir];
    toRegion.setRange(a_idir, coord);

    int isign = sign(a_side);
    RealVect facePos;

    toRegion &= a_state.box();

    Box fromRegion = toRegion;
    fromRegion.shift(a_idir, -isign);

    a_state.copy(a_state, fromRegion, 0, toRegion, 0, a_state.nComp());

    if (!a_homogeneous)
    {
      Real* value = new Real[a_state.nComp()];

      for (BoxIterator bit(toRegion); bit.ok(); ++bit)
      {
        const IntVect& ivTo = bit();
        IntVect ivClose = ivTo - isign*BASISV(a_idir);
        IntVect ivFar = ivTo - 2*isign*BASISV(a_idir);


        getDomainFacePosition(facePos, ivClose, a_dx, a_idir, a_side);

        a_value(facePos.dataPtr(), &a_idir, &a_side, value);
        for (int icomp = a_interval.begin(); icomp <= a_interval.end(); icomp++)
        {

          //Enforce the exact value on the boundary face
          a_state(ivTo, icomp) = value[icomp];

        }
      }

      delete[] value;
    }


  }

  /// Apply neumann BCs (edge centred)
  virtual void NeumEdgeBC(FArrayBox&      a_state,
                          const Box&      a_valid,
                          Real            a_dx,
                          bool            a_homogeneous,
                          const BCValueHolder&   a_valueA,
                          int             a_dir,
                          Side::LoHiSide  a_side,
                          Interval        a_interval)
  {
    // a_state is FACE-centered on face idir;
    // a_valid is CELL-centered.
    Box toRegion(a_valid);
    toRegion.surroundingNodes(a_dir);
    int coord = toRegion.sideEnd(a_side)[a_dir];
    toRegion.setRange(a_dir, coord);

    BCValueHolder& a_value = (BCValueHolder&)a_valueA;
    int isign = sign(a_side);
    RealVect facePos;

    //    Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
    //      toRegion.shift(a_dir, isign);
    toRegion &= a_state.box();

    Box fromRegion = toRegion;
    fromRegion.shift(a_dir, -isign);

    a_state.copy(a_state, fromRegion, 0, toRegion, 0, a_state.nComp());

    if (!a_homogeneous)
    {
      Real* value = new Real[a_state.nComp()];

      for (BoxIterator bit(toRegion); bit.ok(); ++bit)
      {
        const IntVect& ivTo = bit();
        IntVect ivClose = ivTo - isign*BASISV(a_dir);
        IntVect ivFar = ivTo - 2*isign*BASISV(a_dir);


        getDomainFacePosition(facePos, ivClose, a_dx, a_dir, a_side);

        a_value(facePos.dataPtr(), &a_dir, &a_side, value);
        for (int icomp = a_interval.begin(); icomp <= a_interval.end(); icomp++)
        {
          //            Real phi1 = a_state(ivClose, icomp);
          //            Real phi2 = a_state(ivFar, icomp);

          if (value[icomp] == 0)
          {
            // can do higher order stuff here because I've done the maths
            //              a_state(ivTo, icomp) = (1/5)*(phi2 - 4*phi1);
            a_state(ivTo, icomp) += Real(isign)*a_dx*value[icomp];
          }
          else
          {
            a_state(ivTo, icomp) += Real(isign)*a_dx*value[icomp];
          }
        }
      }

      delete[] value;
    }
  }
};



/// Boundary condition function for fields which are advected and diffused
/**
 * All these fields have the same type of conditions at each boundary, but just different values.
 * This boundary condition can also handle a plume inflow, if required.
 */
class AdvectDiffuseScalarBC: public AbstractScalarBCFunction
{
public:

  /// Each element of the vector corresponds to a different interval
  Vector<Real> m_plumeVal;

  /// BC type on the lo side. First index refers to direction, second to component
  Vector<Vector<int> > m_customLoBC,
  /// BC type on the hi side. First index refers to direction, second to component
  m_customHiBC;

  /// BC value on the lo side. First index refers to direction, second to component
  Vector<Vector<Real> > m_customLoBCVal,

  /// BC value on the hi side. First index refers to direction, second to component
  m_customHiBCVal;



  /// Constructor
  AdvectDiffuseScalarBC(bool a_isDefined,
                        MushyLayerParams a_params,
                        bool a_homogeneous,
                        LevelData<FluxBox>* a_advVel,
                        Real a_dx,
                        Vector<Real>  a_plumeVals,
                        Interval a_interval,
                        Vector<Vector<int> >& a_customLoBC,
                        Vector<Vector<int> >&  a_customHiBC,
                        Vector<Vector<Real> >& a_customLoBCVal,
                        Vector<Vector<Real> >&  a_customHiBCVal) : AbstractScalarBCFunction(a_isDefined,
                                                                                            a_params, a_homogeneous,
                                                                                            a_advVel, a_dx, a_interval),
                                                                                            m_plumeVal(a_plumeVals),
                                                                                            m_customLoBC(a_customLoBC),
                                                                                            m_customHiBC(a_customHiBC),
                                                                                            m_customLoBCVal(a_customLoBCVal),
                                                                                            m_customHiBCVal(a_customHiBCVal)
  { }

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (!m_isDefined)
    {
      MayDay::Error("undefined AdvectDiffuseScalarBC object");
    }

    a_homogeneous = a_homogeneous || m_homogeneous;


    //int order = m_params.m_BCAccuracy;
    // always use 1st order for stability
    int order = 1;

    const Box& domainBox = a_domain.domainBox();

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (!a_domain.isPeriodic(idir))
      {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();
          //          int isign = sign(side);

          if (a_valid.sideEnd(side)[idir] ==
              domainBox.sideEnd(side)[idir])
          {

            //Iterate over components
            for (int comp = m_interval.begin(); comp <= m_interval.end(); comp++)
            {
              Real boundaryValue;


              int bcType;

              if (side == Side::Lo)
              {

                if (m_customLoBC.size() > 0)
                {
                  bcType = m_customLoBC[idir][comp];
                }
                else
                {
                  bcType = m_params.bcTypeScalarLo[idir] ;
                }

                boundaryValue = m_customLoBCVal[idir][comp];

              }
              else
              {
                if (m_customHiBC.size() > 0)

                {
                  bcType = m_customHiBC[idir][comp];
                }
                else
                {
                  bcType = m_params.bcTypeScalarHi[idir];
                }

                boundaryValue = m_customHiBCVal[idir][comp];
              }

              if (bcType == PhysBCUtil::FixedTemperature || bcType == PhysBCUtil::TemperatureFlux || bcType == PhysBCUtil::TemperatureFluxRadiation)
              {
//                CH_assert(a_state.nComp() == 2);

                if (a_state.nComp() == 1)
                {
                  int a_comp = 0;
                  // if for some reason we only have one component, just apply extrapolation bcs to a_state

                  ExtraBC(a_state, a_valid,
                                      idir, side, order, a_comp);

                }
                else
                {

                  // This is a special case

                  int Hcomp = 0;
                  int Ccomp = 1;


                  CH_TIME("BCFunctions::TemperatureBC");
                  int isign = sign(side);

                  NonlinearTemperatureBC* residual;

                  // default values
                  Real a = 0.0, b=0.0, T_ref=0.0, no_flux_position_limit=0.0, flux=0.0;

                  no_flux_position_limit = m_params.m_bc_noFluxLimit.getBC(idir,  side);

                  switch (bcType)
                  {
                    case PhysBCUtil::FixedTemperature:
                      T_ref = boundaryValue;
                      b = 1.0;  // the (+) sign is important here in terms of driving the bc lower if our estimated boundary temperature is too high (and vice versa)

                      break;

                    case PhysBCUtil::TemperatureFlux:

                      flux = boundaryValue;
                      a = 1.0;

                      break;

                    case PhysBCUtil::TemperatureFluxRadiation:
                      a = m_params.m_bc_a.getBC(idir, side);;
                      b = m_params.m_bc_b.getBC(idir, side);
                      flux = boundaryValue;
                      T_ref = m_params.m_bc_bTref.getBC(idir, side);

                      break;

                    default:

                      //MayDay::Error("Unknown bcType specified");

                      break;
                  }

                  residual = new NonlinearTemperatureBCRobin(m_params, m_dx,
                                                             a, // coefficient of dt/d(x, z)
                                                             b, // coefficient of thermal radiation term
                                                             flux, // flux
                                                             T_ref); //reference temperature for thermal radiation (dummy value)

                  NonlinearBCSolver* nonlinearSolver;

                  if (m_params.bc_nonlinear_solve_method == NonlinearBCSolveMethods::newton)
                  {
                    nonlinearSolver = new NonlinearBCSolverNewton(residual, m_params.max_bc_iter, m_params.max_bc_residual);
                  }
                  else if (m_params.bc_nonlinear_solve_method == NonlinearBCSolveMethods::picard)
                  {
                    nonlinearSolver = new NonlinearBCSolverPicard(residual, m_params.max_bc_iter, m_params.max_bc_residual, m_params.bc_relax_coeff);
                  }

                  Box grownBox(a_valid);


                  Box toRegion = adjCellBox(grownBox, idir, side, 1);
                  toRegion &= a_state.box();

                  for (BoxIterator bit = BoxIterator(toRegion); bit.ok(); ++bit)
                  {
                    CH_TIME("BCFunctions::TemperatureBC::loopOverBox");

                    IntVect ivto = bit();
                    IntVect iv_interior = ivto - isign*BASISV(idir);

                    // Assume first order BCs
                    // Assume a_state contains enthalpy and bulk concentration components
                    Real interior_enthalpy = a_state(iv_interior, Hcomp);
                    Real interior_bulk_concentration = a_state(iv_interior, Ccomp);
                    Real interior_porosity =  m_params.computePorosity(interior_enthalpy, interior_bulk_concentration);

                    // Compute temperature on the boundary
                    bool legacy_fixed_temp_bc = false;

                    if (bcType == PhysBCUtil::FixedTemperature && legacy_fixed_temp_bc)
                    {
                      CH_TIME("BCFunctions::FixedTemperatureBC");

                      Real boundary_temperature_value;
                      Real bcVal;

                      boundary_temperature_value = boundaryValue;
                      if (comp == Hcomp)
                      {
                        bcVal = interior_porosity*m_params.stefan + boundary_temperature_value;
                      }
                      else
                      {
                        // This is only true in a mushy layer!!
                        // need to think a bit more carefully about what the bulk concentration boundary condition
                        // should be given a fixed temperature/temperature flux
                        MayDay::Error("Can't apply temperature bcs to bulk concentration yet");
                        bcVal = interior_porosity*boundary_temperature_value + m_params.compositionRatio*(1-interior_porosity);
                      }

                      // enforce dirichlet condition
                      if (a_homogeneous)
                      {
                        bcVal = 0;
                      }

                      // Set boundary value as required
                      // This could be more accurate for flux BCs, maybe fix later
                      a_state(ivto, comp) = 2*bcVal - a_state(iv_interior, comp);


                    }
                    // if (bcType == PhysBCUtil::TemperatureFlux || bcType == PhysBCUtil::TemperatureFluxRadiation)
                    else
                    {
                      CH_TIME("BCFunctions::TemperatureFluxBC");

                      if (comp == Hcomp)
                      {

                        // define perpendicular direction = current dir + 1 unless
                        // current dir is the largest dimension, then go back to 0
                        int perp_dir = (idir == SpaceDim-1) ? 0 : idir + 1;

                        Real face_pos = m_dx*(iv_interior[perp_dir] + 0.5);
                        if (face_pos < no_flux_position_limit || a_homogeneous)
                        {
                          // no flux below z = some value (0.75) for now
  //                        bcVal = 0.0;

                          // No flux boundaries are easy:
                          a_state(ivto, comp) = a_state(iv_interior, comp);
                        }
                        else
                        {

                          CH_TIME("BCFunctions::FixedTempBC::nonlinearsolve");

                          // Compute BC valute to satisfy some boundary condition whose residual is computed by the residual object.

                          // Get the temperature at the interior cell (which won't change) so we can compute fluxes across the boundary
                          Real interior_temperature = m_params.computeTemperature(interior_enthalpy, interior_bulk_concentration);

                          // Initial guess of the enthalpy in the ghost cell is the 0th order extrapolation from inside the domain
                          Real ghost_enthalpy =  interior_enthalpy;

                          // Now solve for the enthalpy in the ghost cell
                          nonlinearSolver->solve(ghost_enthalpy, interior_temperature, interior_bulk_concentration);

                          // Since we have determined what the enthalpy should be in the ghost cell to enforce the flux condition
                          // on temperature, just enforce this value
                          a_state(ivto, comp) = ghost_enthalpy;

                        } // end if greater than/less than some z value
                      }
                      else
                      {

                        MayDay::Error("Can't apply temperature bcs to bulk concentration yet");

                      }


                    } // end if temperature flux condition



                  } // end loop over box

                  // Clean up memory
                  delete residual;


                }


              }
              else if (bcType == PhysBCUtil::MixedBCTemperatureSalinity)
              {

                  // This is a special case

                  int Tcomp = 0;
                  int Slcomp = 1;


                  CH_TIME("BCFunctions::MixedTemperatureSalinityBC");
                  int isign = sign(side);

                  //                                  NonlinearTemperatureBC* residual;

                  // default values
                  Real a = 0.0, b=0.0, T_ref=0.0, no_flux_position_limit=0.0, flux=0.0;

                  no_flux_position_limit = m_params.m_bc_noFluxLimit.getBC(idir,  side);

                  a = m_params.m_bc_a.getBC(idir, side);;
                  b = m_params.m_bc_b.getBC(idir, side);
                  flux = boundaryValue;
                  T_ref = m_params.m_bc_bTref.getBC(idir, side);

                  Box toRegion = adjCellBox(a_valid, idir, side, 1);
                  toRegion &= a_state.box();

                  Real ghostVal = 0.0;

                  for (BoxIterator bit = BoxIterator(toRegion); bit.ok(); ++bit)
                  {
                    // TODO - should write this in fortran for speed
                    CH_TIME("BCFunctions::MixedTemperatureSalinityBC::loopOverBox");

                    IntVect ivto = bit();
                    IntVect iv_interior = ivto - isign*BASISV(idir);

                    ghostVal = (a_state(iv_interior, Tcomp)*(a/m_dx + b/2) + flux - b*T_ref ) / (a/m_dx - b/2);
                    a_state(ivto, Tcomp) = ghostVal;
                    if (a_state.nComp() == 2)
                    {
                      a_state(ivto, Slcomp) = a_state(iv_interior, Slcomp);
                    }

                  }


              }

              else
              {
              // All other BCs (e.g. dirichlet, neumann) are applied here.

                applyCorrectScalarBC(a_state,
                                     m_advVel,
                                     a_valid,
                                     a_homogeneous,
                                     boundaryValue,
                                     m_plumeVal[comp],
                                     0,
                                     idir,
                                     side,
                                     a_dx,
                                     order,
                                     bcType,
                                     m_params.plumeBounds,
                                     m_dx,
                                     comp);

              } // end loop over different BCs

            } // loop over components

          } // end if this box is at the edge of the domain
        } // end loop over sides
      } // end if domain isn't periodic in this direction
    } // loop over directions

  } // end operator() function
};



/// Apply extrapolation boundary conditions to a flux component
/**
 * Extrapolation BCs for the specified component of a flux box.
 */
// ---------------------------------------------------------------
class BasicFluxExtrapBCFunction: public BCFunction
{
public:

  /// Component of velocity to apply BCs to
  int m_comp;

  // Default constructor
  BasicFluxExtrapBCFunction()
  :
    m_comp(-1)
  {
  }

  /// Full constructor
  BasicFluxExtrapBCFunction(
      int  a_comp)
  :
    m_comp(a_comp)
  {
  }
  /// Apply BC
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {

    int order = 2;

    // a_state is FACE-centered, now in m_comp direction;
    // a_valid is CELL-centered
    if (m_comp >= 0)
    {
      const Box& domainBox = a_domain.domainBox();

      // loop over directions
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (!a_domain.isPeriodic(idir))
        {
          SideIterator sit;
          for (sit.reset(); sit.ok(); ++sit)
          {
            Side::LoHiSide side = sit();
            if (a_valid.sideEnd(side)[idir] ==
                domainBox.sideEnd(side)[idir])
            {

              ExtraBC(a_state, a_valid,
                      idir, side, order, m_comp);
              //						  ExtrapBC(  a_state, a_valid,  idir,   side, order);

            } // if ends match
          } // end iteration over sides
        } // if not periodic in this direction
      } // end iteration over directions
    } // if m_comp >= 0
    else // m_comp < 0
    {
      MayDay::Error("undefined BasicFluxExtrapBCFunction object");
    }
  }

};



/// Boundary conditions for velocity (edge-centred)
/**
 * Sets specified boundary conditions for the velocity.
 * For use on an edge-centred velocity, i.e. a LevelData<FluxBox>
 */

// ---------------------------------------------------------------
class BasicECVelBCFunction: public AbstractFaceBCFunction
{
public:
  /// Is there non-zero viscosity?
  bool m_isViscous;

  /// Currently unused
  Interval m_interval;

  /// Velocity values to enforce on the domain boundaries
  LevelData<FluxBox>* m_velBCvals;

  /// Full constructor
  BasicECVelBCFunction(bool a_isHomogeneous,
                       bool a_isViscous,
                       int  a_comp,
                       const Interval& a_interval,
                       MushyLayerParams a_params,
                       LevelData<FluxBox>* a_velocityBCVals = nullptr) : AbstractFaceBCFunction(a_isHomogeneous,
                                                                                             a_comp,
                                                                                             a_params),
                                                                                             m_isViscous(a_isViscous),
                                                                                             m_interval(a_interval),
                                                                                             m_velBCvals(a_velocityBCVals)
  {
  }

  /// Apply BC
  void operator()(FArrayBox&           a_state,
                  const Box&           a_valid,
                  const ProblemDomain& a_domain,
                  Real                 a_dx,
                  bool                 a_homogeneous)
  {
    // Stick with 1st order for stability
    int order = 1;

    // a_state is FACE-centered, now in m_comp direction;
    // a_valid is CELL-centered
    if (m_comp >= 0)
    {
      const Box& domainBox = a_domain.domainBox();

      RefCountedPtr<ConstValueFunction>
      zeroFunc(new ConstValueFunction(0.0, a_state.nComp()));
      RefCountedPtr<ConstValueFunction>
      bcValueFunc(new ConstValueFunction(m_params.inflowVelocity, a_state.nComp()));

      RefCountedPtr<InflowValueFunction>
      inflowBCValueFunc(new InflowValueFunction(0.0, a_state.nComp(),
                                                m_params.inflowVelocity,
                                                m_params.plumeBounds[0], m_params.plumeBounds[1]));

      // loop over directions
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (!a_domain.isPeriodic(idir))
        {
          SideIterator sit;
          for (sit.reset(); sit.ok(); ++sit)
          {
            Side::LoHiSide side = sit();
            if (a_valid.sideEnd(side)[idir] ==
                domainBox.sideEnd(side)[idir])
            {

              int bctype = m_params.getVelBCType(idir, side);

              switch (bctype)
              {
                case PhysBCUtil::SolidWall :
                {
                  // normal velocity BCs:
                  // always no-flow
                  if (idir == m_comp)
                  {
                    DiriEdgeBC(a_state, a_valid, a_dx,
                               a_homogeneous,
                               BCValueHolder(zeroFunc),
                               idir, side);
                  }
                  else
                  {
                    // tangential BCs:
                    // no-slip if viscous, extrap if inviscid
                    if (m_isViscous)
                    {

                      // need to fake this a bit
                      // want to use Cell-centered
                      // BC function to set tangential BC's
                      // on face-centered data. Do this by
                      // shifting valid-region to
                      // face-centering.  DiriBC function
                      // only really cares that valid box
                      // and state box have the same centering
                      // (DFM -- 9/22/08)
                      Box validFace(a_valid);
                      validFace.surroundingNodes(m_comp);
                      DiriBC(a_state, validFace, a_dx,
                             a_homogeneous,
                             BCValueHolder(zeroFunc),
                             idir, side, order);
                    }
                    else // inviscid
                    {
                      order  = 2;
                      ExtraBC(a_state, a_valid,
                              idir, side, order);
                    }
                  } // end if tangential
                  break;
                }
                case PhysBCUtil::Inflow :
                {
                  // this will most likely get overwritten in a derived class
                  if (!m_isHomogeneous && idir == m_comp)
                  {
                    DiriEdgeBC(a_state, a_valid, a_dx,
                               a_homogeneous,
                               BCValueHolder(bcValueFunc),
                               idir, side);
                  }
                  else
                  {


                    // See DFM comment above about valid regions and box centering
                    Box validFace(a_valid);
                    validFace.surroundingNodes(m_comp);
                    DiriBC(a_state, validFace, a_dx,
                           a_homogeneous,
                           BCValueHolder(zeroFunc),
                           idir, side, order);

                  }
                  break;
                }
                case PhysBCUtil::VelInflowPlume :
                {

                  // normal velocity BCs:
                  // always no-flow
                  if (idir == m_comp)
                  {
                    //                                      DiriEdgeBC(a_state, a_valid, a_dx,
                    //                                                 a_homogeneous,
                    //                                                 BCValueHolder(zeroFunc),
                    //                                                 idir, side);

                    DiriEdgeVariableBC(a_state, a_valid, a_dx,
                                       a_homogeneous,
                                       BCValueHolder(inflowBCValueFunc),
                                       idir, side);

                    //                    int temp=0;
                  }
                  else
                  {
                    // tangential BCs:
                    // no-slip if viscous, extrap if inviscid
                    if (m_isViscous)
                    {

                      // need to fake this a bit
                      // want to use Cell-centered
                      // BC function to set tangential BC's
                      // on face-centered data. Do this by
                      // shifting valid-region to
                      // face-centering.  DiriBC function
                      // only really cares that valid box
                      // and state box have the same centering
                      // (DFM -- 9/22/08)
                      Box validFace(a_valid);
                      validFace.surroundingNodes(m_comp);
                      DiriBC(a_state, validFace, a_dx,
                             a_homogeneous,
                             BCValueHolder(zeroFunc),
                             idir, side, order);
                    }
                    else // inviscid
                    {
                      order  = 2;
                      ExtraBC(a_state, a_valid,
                              idir, side, order, m_comp);
                      //                                                                          ExtrapBC(  a_state, a_valid,  idir,   side, order);
                    }
                  } // end if tangential

                  break;
                }
                case PhysBCUtil::Outflow :
                {
                  // this is set to whatever it is set to by MAC projection
                  // NoOp


                  break;
                }
                case PhysBCUtil::OutflowPressureGrad :
                {

                  // Set a_state  = m_velBCvals on face
                  // a_state is FACE-centered on face idir;
                  // a_valid is CELL-centered.

                  // Need to find the correct dataiterator
                  DataIterator dit = m_velBCvals->dataIterator();

                  //                  FluxBox* velBCVals = nullptr;
                  FArrayBox* velBCValComp = nullptr;

                  Box domBox2(domainBox);
                  domBox2.surroundingNodes(idir);

                  Box domBox3(domainBox);
                  domBox3.surroundingNodes(m_comp);

                  Box stateBox = a_state.box();
                  stateBox &= domBox3;

                  for (dit.reset(); dit.ok(); ++dit)
                  {
                    Box velBox = (*m_velBCvals)[dit][m_comp].box();
                    velBox &= domBox3;

                    IntVect se = velBox.smallEnd();
                    IntVect be = velBox.bigEnd();
                    //                    pout() << velBox.smallEnd() << " - " << velBox.bigEnd() << endl;
                    //                    pout() << stateBox.smallEnd() << " - " << stateBox.bigEnd() << endl;

                    if (velBox == stateBox || stateBox.contains(velBox))
                    {
                      //                      velBCVals = &(*m_velBCvals)[dit];
                      velBCValComp = &(*m_velBCvals)[dit][m_comp];
                    }

                  }

                  if (velBCValComp != nullptr)
                  {

                    Box toRegion(a_valid);

                    if (idir == m_comp)
                    {
                      toRegion.surroundingNodes(idir);
                      int coord = toRegion.sideEnd(side)[idir];
                      toRegion.setRange(idir, coord);
                      toRegion &= a_state.box();
                    }
                    else
                    {
                      toRegion.surroundingNodes(m_comp);
                      toRegion = adjCellBox(toRegion, idir, side, 1);
                      toRegion &= a_state.box();
                    }

                    Box sbox = a_state.box();
                    Box bc_box = (*velBCValComp).box();



                    a_state.copy(*velBCValComp, toRegion, 0, toRegion, 0, a_state.nComp());


                  }

                  break;
                }

                case PhysBCUtil::VelInflowOutflow :
                case PhysBCUtil::OutflowNormal :
                case PhysBCUtil::PressureHead :
                {
                  // Normal component - zero gradient
                  if (idir == m_comp)
                  {
                    NeumEdgeBC(a_state, a_valid, a_dx,
                               a_homogeneous,
                               BCValueHolder(zeroFunc),
                               idir, side, Interval(0,0));

                  }
                  else
                  {
                    // tangential components = 0
                    Box validFace(a_valid);
                    validFace.surroundingNodes(m_comp);
                    DiriBC(a_state, validFace, a_dx,
                           a_homogeneous,
                           BCValueHolder(zeroFunc),
                           idir, side, order);
                  }
                  //

                  break;
                }
                case PhysBCUtil::noShear :
                case PhysBCUtil::Symmetry :
                {
                  // normal velocity BC's still no-flow
                  if (idir == m_comp)
                  {
                    DiriEdgeBC(a_state, a_valid, a_dx,
                               a_homogeneous,
                               BCValueHolder(zeroFunc),
                               idir, side);
                  }
                  else
                  {
                    // tangential BC is do-nothing BC???
                    // NoOp
                  } // end if tangential
                  break;
                }
                default :
                {
                  MayDay::Error("BasicECVelBCFunction - unknown BC type");
                }
              } // end switch

            } // if ends match
          } // end iteration over sides
        } // if not periodic in this direction
      } // end iteration over directions
    } // if m_comp >= 0
    else // m_comp < 0
    {
      MayDay::Error("undefined BasicECVelBCFunction object");
    }
  }


};


/// Boundary condition for source terms
/**
 *  sets boundary conditions on source terms. At the moment, this is just 0th order extrapolation from
 the interior, so it's real simple. Fills interior as well as domain ghost cells.
 */
class ExtrapolationBCFunction: public BCFunction
{
public:

  /// Is object defined?
  bool m_isDefined;

  /// Order of accuracy (1st, 2nd etc.)
  int m_order;

  // Default constructor
  ExtrapolationBCFunction()
  :
    m_isDefined(false), m_order(0)
  {
  }
  /// Full constructor
  ExtrapolationBCFunction(bool a_isDefined, int a_order=0)
  :
    m_isDefined(a_isDefined), m_order(a_order)
  {
  }
  /// Apply BC
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_isDefined)
    {
      int nComp = a_state.nComp();
      const Box& grownBox = a_state.box();

      for (int idir = 0; idir < SpaceDim; idir++)
      {
        // only apply BCs if not periodic!
        // changing this so we still apply BCs if periodic, in order to fill interior CF
        // boundaries in the periodic direction. Cells on the domain boundary are then
        // replaced later by exchange calls
        //        if (!a_domain.isPeriodic(idir))
        //        {
        // Huh, these will be undefined if conditions below don't hold.
        Box loBox, hiBox;

        if (grownBox.smallEnd(idir) < a_valid.smallEnd(idir))
        {
          loBox = grownBox;
          loBox.setBig(idir, a_valid.smallEnd(idir)-1);
        }
        if (grownBox.bigEnd(idir) > a_valid.bigEnd(idir))
        {
          hiBox = grownBox;
          hiBox.setSmall(idir, a_valid.bigEnd(idir)+1);
        }

        if (m_order == 0)
        {

          FORT_EXTRAPOLATIONBC(CHF_FRA(a_state),
                               CHF_BOX(a_valid),
                               CHF_BOX(loBox),
                               CHF_BOX(hiBox),
                               CHF_INT(idir),
                               CHF_INT(nComp));

        }
        else if (m_order == 1)
        {
          FORT_EXTRAPOLATIONBC1(CHF_FRA(a_state),
                                CHF_BOX(a_valid),
                                CHF_BOX(loBox),
                                CHF_BOX(hiBox),
                                CHF_INT(idir),
                                CHF_INT(nComp));
        }
      } // end loop over all directions
      //      }
    }
    else
    {
      MayDay::Error("undefined ViscousBCFunction object");
    }
  }
};

/// Boundary condition for velocity (cell-centred)
class BasicCCVelBCFunction: public BCFunction
{
public:

  /// Component of velocity to apply BCs to
  int m_comp;
  /// Always enforce homogeneous form of BCs?
  bool m_isHomogeneous;
  /// Is there non-zero viscosity?
  bool m_isViscous;
  /// Specifies where in a_state we will find the component we should apply BCs to
  Interval m_interval;
  /// Physical parameters for the problem
  MushyLayerParams m_params;

  /// Default constructor
  BasicCCVelBCFunction()
  :
    m_comp(-1),  m_isHomogeneous(false), m_isViscous(false)
  {
  }

  /// Full constructor
  BasicCCVelBCFunction(bool a_isHomogeneous,
                       bool a_isViscous,
                       int  a_comp,
                       const Interval& a_interval,
                       MushyLayerParams a_params)
  : m_comp(a_comp),
    m_isHomogeneous(a_isHomogeneous),
    m_isViscous(a_isViscous),
    m_interval(a_interval),
    m_params(a_params)
  {
  }
  /// Apply BC
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_comp >= 0)
    {
      const Box& domainBox = a_domain.domainBox();
      //			PhysBCUtil bcInfo(m_params); // sets BCs from ParmParse table

      FArrayBox aliasStateFab(m_interval, a_state);

      // aliasStateFab only has 1 component, in index 0
      int fab_comp = 0;

      RefCountedPtr<ConstValueFunction>
      zeroFunc(new ConstValueFunction(0.0, aliasStateFab.nComp()));
      RefCountedPtr<ConstValueFunction>
      bcValueFunc(new ConstValueFunction(m_params.inflowVelocity, aliasStateFab.nComp()));

      int order = m_params.m_BCAccuracy;
      //      order = 1;

      // loop over directions
      for (int idir = 0; idir < SpaceDim; idir++)
      {

        if (!a_domain.isPeriodic(idir))
        {
          SideIterator sit;
          for (sit.reset(); sit.ok(); ++sit)
          {
            Side::LoHiSide side = sit();
            if (a_valid.sideEnd(side)[idir] ==
                domainBox.sideEnd(side)[idir])
            {
              int bctype = m_params.getVelBCType(idir, side);

              switch (bctype)
              {
                case PhysBCUtil::SolidWall :
                {
                  // normal velocity BCs:
                  // always no-flow
                  if (idir == m_comp)
                  {
                    order = 1;
                    DiriBC(aliasStateFab, a_valid, a_dx,
                           a_homogeneous,
                           BCValueHolder(zeroFunc),
                           idir, side, order);
                  }
                  else
                  {
                    // tangential BCs:
                    // no-slip if viscous, extrap if inviscid
                    if (m_isViscous)
                    {
                      order = 1;
                      DiriBC(aliasStateFab, a_valid, a_dx,
                             a_homogeneous,
                             BCValueHolder(zeroFunc),
                             idir, side, order);
                    }
                    else // inviscid
                    {

                      order = 1;
                      ExtraBC(aliasStateFab, a_valid,
                              idir, side, order, fab_comp);
                      //									  ExtrapBC(aliasStateFab, a_valid, idir, side, order);
                    }
                  } // end if tangential
                  break;
                }
                case PhysBCUtil::Inflow :
                {
                  // this will most likely get overwritten in a derived class
                  if (!m_isHomogeneous && idir == m_comp)
                  {
                    order = 1;
                    DiriBC(aliasStateFab, a_valid, a_dx,
                           a_homogeneous,
                           BCValueHolder(bcValueFunc),
                           idir, side, order);
                  }
                  else
                  {
                    order = 1;
                    DiriBC(aliasStateFab, a_valid, a_dx,
                           a_homogeneous,
                           BCValueHolder(zeroFunc),
                           idir, side, order);
                  }
                  break;
                }
                case PhysBCUtil::VelInflowPlume :
                {

                  // normal velocity BCs:
                  // always no-flow
                  if (idir == m_comp)
                  {

                    //                    DiriBC(aliasStateFab, a_valid, a_dx,
                    //                           a_homogeneous,
                    //                           BCValueHolder(zeroFunc),
                    //                           idir, side, order);

                    order = 1;
                    OnlyInflowBC(aliasStateFab,
                                 a_valid,
                                 a_homogeneous,
                                 m_params.inflowVelocity, // vertical velocity for plume
                                 0, // vertical velocity away from plume (no flow)
                                 m_params.plumeBounds,
                                 a_dx,
                                 idir,
                                 side,
                                 order,
                                 m_comp);

                  }
                  else
                  {
                    // tangential BCs:
                    // no-slip if viscous, extrap if inviscid
                    if (m_isViscous)
                    {
                      order = 1;
                      DiriBC(aliasStateFab, a_valid, a_dx,
                             a_homogeneous,
                             BCValueHolder(zeroFunc),
                             idir, side, order);
                    }
                    else // inviscid
                    {
                      order = 1;
                      ExtraBC(aliasStateFab, a_valid,
                              idir, side, order, fab_comp);
                      //                                                                          ExtrapBC(aliasStateFab, a_valid, idir, side, order);
                    }
                  } // end if tangential




                  break;
                }
                case PhysBCUtil::Outflow :
                {

                  order = 1;
                  ExtrapBC(aliasStateFab, a_valid, idir, side, order, fab_comp);

                  break;
                }
                case PhysBCUtil::VelInflowOutflow :
                case PhysBCUtil::OutflowNormal :
                case PhysBCUtil::OutflowPressureGrad :
                case PhysBCUtil::PressureHead :
                {
                  // normal component - no gradient
                  if (idir == m_comp)
                  {
                    //                    ExtrapBC(aliasStateFab, a_valid, idir, side, order);
                    NeumBC(aliasStateFab, a_valid, a_dx,
                           a_homogeneous,
                           BCValueHolder(zeroFunc),
                           idir, side);
                  }
                  else
                  {
                    // tangential components = 0
                    order = 1;
                    DiriBC(aliasStateFab, a_valid, a_dx,
                           a_homogeneous,
                           BCValueHolder(zeroFunc),
                           idir, side, order);
                  }

                  break;
                }
                case PhysBCUtil::noShear :
                case PhysBCUtil::Symmetry :
                {
                  // normal velocity BC's still no-flow
                  if (idir == m_comp)
                  {
                    order = 1;
                    DiriBC(aliasStateFab, a_valid, a_dx,
                           a_homogeneous,
                           BCValueHolder(zeroFunc),
                           idir, side, order);
                  }
                  else
                  {
                    NeumBC(aliasStateFab, a_valid, a_dx,
                           a_homogeneous,
                           BCValueHolder(zeroFunc),
                           idir, side);
                  } // end if tangential
                  break;
                }
                default :
                {
                  MayDay::Error("BasicCCVelBCFunction - unknown BC type");
                }
              } // end switch


              // Fill outer ghost cells?
              Box grownValid(a_valid);
              grownValid.grow(1);

              Box shrunkFabBox = aliasStateFab.box();
              shrunkFabBox.grow(-1);

              // If shrunk fab box contains grownValid, that means we have a set of cells
              // outside grown valid which we need to fill
              while(shrunkFabBox.contains(grownValid))
              {
                ExtraBC(aliasStateFab, grownValid, idir, side,
                        0, fab_comp); //0th order
                grownValid.grow(1);
              }

            } // if ends match
          } // end iteration over sides
        } // if not periodic in this direction
      } // end iteration over directions
    } // if m_comp >= 0
    else // m_comp < 0
    {
      MayDay::Error("undefined BasicCCVelBCFunction object");
    }
  }
};




/// Boundary conditions for the streamfunction \f$ \psi \f$
class StreamFunctionBC: public AbstractScalarBCFunction
{
public:
  /// Apply BCs for streamfunction
  StreamFunctionBC(bool a_isDefined,
                   MushyLayerParams a_params,
                   bool a_homogeneous,
                   LevelData<FluxBox>* a_advVel,
                   Real a_dx) : AbstractScalarBCFunction(a_isDefined,
                                                         a_params,
                                                         a_homogeneous,
                                                         a_advVel,
                                                         a_dx)
{}


  /// Apply BC
  void operator()(FArrayBox&           a_state,
                  const Box&           a_valid,
                  const ProblemDomain& a_domain,
                  Real                 a_dx,
                  bool                 a_homogeneous)
  {
    if (m_isDefined)
    {
      const Box& domainBox = a_domain.domainBox();

      RefCountedPtr<ConstValueFunction>
      zeroFunc(new ConstValueFunction(0.0, a_state.nComp()));



      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (!a_domain.isPeriodic(idir))
        {
          SideIterator sit;
          for (sit.reset(); sit.ok(); ++sit)
          {
            Side::LoHiSide side = sit();
            if (a_valid.sideEnd(side)[idir] ==
                domainBox.sideEnd(side)[idir])
            {

              int bctype = m_params.getVelBCType(idir, side);
              switch(bctype)
              {
                case PhysBCUtil::SolidWall :
                case PhysBCUtil::Inflow :
                case PhysBCUtil::noShear :
                case PhysBCUtil::Symmetry : // think this is right for symmetry
                {
                  // Constant streamfunction on solid walls
                  int order = 1;
                  DiriBC(a_state, a_valid, a_dx,
                         a_homogeneous,
                         BCValueHolder(zeroFunc),
                         idir, side, order);
                  break;
                }
                case PhysBCUtil::Outflow :
                case PhysBCUtil::OutflowNormal :
                case PhysBCUtil::VelInflowOutflow:
                case PhysBCUtil::VelInflowPlume:
                case PhysBCUtil::OutflowPressureGrad :
                case PhysBCUtil::PressureHead :
                {
                  // not sure if this right or not
                  NeumBC(a_state, a_valid, a_dx,
                         a_homogeneous,
                         BCValueHolder(zeroFunc),
                         idir, side);
                  break;
                }



                default :
                {
                  MayDay::Error("StreamFunctionBC - unknown BC type");
                }
              } // end switch

            } // end condition of matching coordinates
          } // end loop over both sides
        } // end if not periodic in this direction
      } // end loop over all directions
    }
    else
    {
      MayDay::Error("undefined StreamFunctionBC object");
    }
  }
};

/// Boundary conditions for pressure \f$ P \f$ (subcycled version)
class BasicPressureBCFunctionSubcycled: public AbstractScalarBCFunction
{
public:
  /// Creator
  BasicPressureBCFunctionSubcycled(bool a_isDefined,
                                   MushyLayerParams a_params,
                                   bool a_homogeneous,
                                   LevelData<FluxBox>* a_advVel,
                                   Real a_dx) : AbstractScalarBCFunction(a_isDefined,
                                                                         a_params,
                                                                         a_homogeneous,
                                                                         a_advVel,
                                                                         a_dx)
{}


  /// Apply BC
  void operator()(FArrayBox&           a_state,
                  const Box&           a_valid,
                  const ProblemDomain& a_domain,
                  Real                 a_dx,
                  bool                 a_homogeneous)
  {
    if (m_isDefined)
    {
      const Box& domainBox = a_domain.domainBox();

      RefCountedPtr<ConstValueFunction> zeroFunc(new ConstValueFunction(0.0, a_state.nComp()));




      RefCountedPtr<InflowValueFunction> inflowFunc(new InflowValueFunction(0.0, a_state.nComp(),
                                                                            -m_params.inflowVelocity,
                                                                            m_params.plumeBounds[0],
                                                                            m_params.plumeBounds[1]));
      RefCountedPtr<PressureInflowValueFunction> inflowFuncNeum(new PressureInflowValueFunction(0.0, a_state.nComp(),
                                                                                                -m_params.inflowVelocity,
                                                                                                m_params.plumeBounds[0],
                                                                                                m_params.plumeBounds[1], &m_params));

      int order = m_params.m_pressureBCAccuracy;

      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (!a_domain.isPeriodic(idir))
        {
          SideIterator sit;
          for (sit.reset(); sit.ok(); ++sit)
          {
            Side::LoHiSide side = sit();
            if (a_valid.sideEnd(side)[idir] ==
                domainBox.sideEnd(side)[idir])
            {

              int bctype = m_params.getVelBCType(idir, side);
              switch(bctype)
              {
                case PhysBCUtil::SolidWall :
                case PhysBCUtil::Inflow :
                case PhysBCUtil::noShear :
                case PhysBCUtil::Symmetry :

                {
                  NeumBC(a_state, a_valid, a_dx,
                         a_homogeneous,
                         BCValueHolder(zeroFunc),
                         idir, side);
                  break;
                }
                case PhysBCUtil::Outflow :
                case PhysBCUtil::OutflowNormal :
                case PhysBCUtil::OutflowPressureGrad :
                {
                  // Fix as first order for stability ?
                      //                  int order = 1;
                  DiriBC(a_state, a_valid, a_dx,
                         a_homogeneous,
                         BCValueHolder(zeroFunc),
                         idir, side, order);


                  break;
                }
                case PhysBCUtil::VelInflowOutflow:
                {
                  //                  int order = 1;
                  PressureInflowOutflow(a_state, a_valid, a_dx,
                                        a_homogeneous,
                                        m_advVel, // advection velocity
                                        0, // dirichlet value (outflow)
                                        0, // neumann value (inflow)
                                        idir, side, order, m_dx);



                  break;
                }

                case PhysBCUtil::PressureHead :
                {
                  Real pressureVal = 0.0;

                  if (side == Side::Lo)
                  {
                    pressureVal = m_params.bcValPressureLo[idir];
                  }
                  else
                  {
                    pressureVal = m_params.bcValPressureHi[idir];
                  }

                  RefCountedPtr<ConstValueFunction> pressureHead(new ConstValueFunction(pressureVal, a_state.nComp()));

                  DiriBC(a_state, a_valid, a_dx,
                         a_homogeneous,
                         BCValueHolder(pressureHead),
                         idir, side, order);

                  break;
                }


                case PhysBCUtil::VelInflowPlume:
                {

                  NeumBC(a_state, a_valid, a_dx,
                         a_homogeneous,
                         //                         BCValueHolder(inflowFuncNeum),
                         BCValueHolder(zeroFunc),
                         idir, side);
                  break;
                }
                default :
                {
                  MayDay::Error("BasicPressureBCFunction - unknown BC type");
                }
              } // end switch

            } // end condition of matching coordinates
          } // end loop over both sides
        } // end if not periodic in this direction
      } // end loop over all directions
    }
    else
    {
      MayDay::Error("undefined BasicPressureBCFunction object");
    }
  }
};


// ---------------------------------------------------------------
/// Boundary condition for face centered porosity and permeability
class BasicPorosityPermeabilityFaceBCFunction: public AbstractFaceBCFunction
{
public:

  /// Value in plume
  Real m_plumeVal;

  /// BC type on the lo side.
  Vector<Vector <int> > m_bcTypeLo,

  /// BC type on the hi side.
  m_bcTypeHi;

  /// BC type on the lo side.
  Vector<Vector <Real> > m_bcValLo,

  /// BC type on the hi side.
  m_bcValHi;

  /// Full constructor
  BasicPorosityPermeabilityFaceBCFunction(bool a_isHomogeneous,
                                          int  a_comp,
                                          MushyLayerParams a_params,
                                          Real a_plumeVal,
                                          Vector<Vector <int> > a_bcTypeLo,
                                          Vector<Vector <int> > a_bcTypeHi,
                                          Vector<Vector <Real> > a_bcValLo,
                                          Vector<Vector <Real> > a_bcValHi
  ):
    AbstractFaceBCFunction(a_isHomogeneous,
                           a_comp,
                           a_params) , m_plumeVal(a_plumeVal),
                           m_bcTypeLo(a_bcTypeLo),
                           m_bcTypeHi(a_bcTypeHi),
                           m_bcValLo(a_bcValLo),
                           m_bcValHi(a_bcValHi)

  {  }

  /// Apply BC
  void operator()(FArrayBox&           a_state,
                  const Box&           a_valid,
                  const ProblemDomain& a_domain,
                  Real                 a_dx,
                  bool                 a_homogeneous)
  {


    // a_state is FACE-centered, now in m_comp direction;
    // a_valid is CELL-centered
    if (m_comp >= 0)
    {
      const Box& domainBox = a_domain.domainBox();


      // loop over directions
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (!a_domain.isPeriodic(idir))
        {
          SideIterator sit;
          for (sit.reset(); sit.ok(); ++sit)
          {
            Side::LoHiSide side = sit();

            if (a_valid.sideEnd(side)[idir] ==
                domainBox.sideEnd(side)[idir])
            {

              int bcType = (side == Side::Lo) ? m_bcTypeLo[idir][0] :  m_bcTypeHi[idir][0];
              Real bcVal = (side == Side::Lo) ? m_bcValLo[idir][0] :  m_bcValHi[idir][0];

              RefCountedPtr<ConstValueFunction> bc(new ConstValueFunction(bcVal, a_state.nComp()));


              if (bcType == PhysBCUtil::Dirichlet)
              {
                if (idir == m_comp)
                {
                  DiriEdgeBC(a_state, a_valid, a_dx,
                             a_homogeneous,
                             BCValueHolder(bc),
                             idir, side);
                }
                else
                {
                  // These BCs don't matter
                }

              }
              else if (bcType == PhysBCUtil::Neumann)
              {
                if (idir == m_comp)
                {
                  Box validFace(a_valid);
                  validFace.surroundingNodes(m_comp);

                  NeumEdgeBC(a_state, validFace, a_dx,
                             a_homogeneous,
                             BCValueHolder(bc),
                             idir, side, Interval(0,0));
                }
                else
                {

                  // Don't think these BCs matter?
                }
              }



            } // if ends match
          } // end iteration over sides
        } // if not periodic in this direction
      } // end iteration over directions
    } // if m_comp >= 0
    else // m_comp < 0
    {
      MayDay::Error("undefined BC object");
    }
  }
};



/// Boundary conditions for porosity, \f$ \chi \f$, and permeability \f$\Pi\f$
class BasicPorosityPermeabilityBCFunction: public AdvectDiffuseScalarBC
{


public:
  /// Constructor
  BasicPorosityPermeabilityBCFunction(bool a_isDefined,
                                      MushyLayerParams a_params,
                                      bool a_homogeneous,
                                      LevelData<FluxBox>* a_advVel,
                                      Real a_dx,
                                      Vector<Real> a_plumeVal,
                                      Interval intvl,
                                      Vector<Vector<int> >& a_customLoBC,
                                      Vector<Vector<int> >&  a_customHiBC,
                                      Vector<Vector<Real> >& a_customLoBCVal,
                                      Vector<Vector<Real> >&  a_customHiBCVal) : AdvectDiffuseScalarBC(a_isDefined,
                                                                                                       a_params,
                                                                                                       a_homogeneous,
                                                                                                       a_advVel,
                                                                                                       a_dx,
                                                                                                       a_plumeVal,
                                                                                                       intvl,
                                                                                                       a_customLoBC,
                                                                                                       a_customHiBC,
                                                                                                       a_customLoBCVal,
                                                                                                       a_customHiBCVal)
{ }


  /// Apply BC
  void operator()(FArrayBox&           a_state,
                  const Box&           a_valid,
                  const ProblemDomain& a_domain,
                  Real                 a_dx,
                  bool                 a_homogeneous)
  {

    // assuming we only have one component to fill here
    CH_assert(a_state.nComp() == 1);

    if (!m_isDefined)
    {
      MayDay::Error("undefined BasicPorosityPermeabilityBCFunction object");
    }

    if (m_homogeneous)
    {
      a_homogeneous = m_homogeneous;
    }
    //    RefCountedPtr<ConstValueFunction>
    //    topBC(new ConstValueFunction(permeabilityTop, a_state.nComp()));

    //    RefCountedPtr<ConstValueFunction>
    //    bottomBC(new ConstValueFunction(m_bottomVal[0], a_state.nComp()));
    //
    //    RefCountedPtr<ConstValueFunction>
    //    sideBC(new ConstValueFunction(0, a_state.nComp()));
    //
    //    RefCountedPtr<InflowValueFunction>
    //    topBC(new InflowValueFunction(m_topVal[0],
    //                                  a_state.nComp(),
    //                                  m_plumeVal[0],
    //                                  m_params.plumeBounds[0], m_params.plumeBounds[1]));

    int order = 1;

    const Box& domainBox = a_domain.domainBox();
    int comp = 0;

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (!a_domain.isPeriodic(idir))
      {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();
          //          int isign = sign(side);
          if (a_valid.sideEnd(side)[idir] ==
              domainBox.sideEnd(side)[idir])
          {

            int bcType = (side==Side::Lo) ? m_customLoBC[idir][comp] : m_customHiBC[idir][comp];
            Real bcVal = (side==Side::Lo) ? m_customLoBCVal[idir][comp] : m_customHiBCVal[idir][comp];

            RefCountedPtr<ConstValueFunction>
            bc(new ConstValueFunction(bcVal, a_state.nComp()));


            if (bcType == 0)
            {


              DiriBC(a_state,
                     a_valid,
                     a_dx,
                     a_homogeneous,
                     BCValueHolder(bc),
                     idir,
                     side,
                     order);

            }
            else if(bcType == 1)
            {

              NeumBC(a_state,
                     a_valid,
                     a_dx,
                     a_homogeneous,
                     BCValueHolder(bc),
                     idir,
                     side);
            }
            //            if (idir == 0)
              //            {
            //              //x direction
            //              NeumBC(a_state,
            //                     a_valid,
            //                     a_dx,
            //                     a_homogeneous,
            //                     BCValueHolder(sideBC),
            //                     idir,
            //                     side);
            //            }
            //            else if(idir == 1)
            //            {
            //
            //              //y direction
            //
            //              if (isign > 0)
            //              {
            //                DiriBC(a_state,
            //                       a_valid,
            //                       a_dx,
            //                       a_homogeneous,
            //                       BCValueHolder(topBC),
            //                       idir,
            //                       side,
            //                       order);
            //
            //              }
            //              else if (isign < 0)
            //              {
            //                DiriBC(a_state,
            //                       a_valid,
            //                       a_dx,
            //                       a_homogeneous,
            //                       BCValueHolder(bottomBC),
            //                       idir,
            //                       side,
            //                       order);
            //              }
          }

        }
      }
    }
  }


};


/// Boundary condition for reflux correction
/**
 * Sets correction = 0 on all boundaries
 */
class BasicRefluxCorrBCFunction: public BCFunction
{
public:

  /// Is object defined?
  bool m_isDefined;
  /// Physical parameters for the problem
  MushyLayerParams m_params;

  /// Default constructor
  BasicRefluxCorrBCFunction()
  :
    m_isDefined(false)
  {
  }
  /// Full constructor
  BasicRefluxCorrBCFunction(bool a_isDefined,
                            MushyLayerParams a_params)
  :
    m_isDefined(a_isDefined),
    m_params(a_params)
  {
  }


  /// Apply BC
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (!m_isDefined)
    {
      MayDay::Error("undefined BasicRefluxCorrBCFunction object");
    }

    RefCountedPtr<ConstValueFunction>
    zeroBC(new ConstValueFunction(0, a_state.nComp()));

    int order = 1;

    const Box& domainBox = a_domain.domainBox();

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (!a_domain.isPeriodic(idir))
      {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();
          //					int isign = sign(side);
          if (a_valid.sideEnd(side)[idir] ==
              domainBox.sideEnd(side)[idir])
          {
            DiriBC(a_state,
                   a_valid,
                   a_dx,
                   a_homogeneous,
                   BCValueHolder(zeroBC),
                   idir,
                   side,
                   order);

          }
        }
      }
    }

  }
};

/// Boundary condition for no normal flux
/**
 * Sets normal flux = 0 on all boundaries
 */
class NoFluxBCFunction: public BCFunction
{
public:

  /// Is object defined?
  bool m_isDefined;
  /// Physical parameters for the problem
  MushyLayerParams m_params;

  /// Default constructor
  NoFluxBCFunction()
  :
    m_isDefined(false)
  {
  }
  /// Full constructor
  NoFluxBCFunction(bool a_isDefined,
                   MushyLayerParams a_params)
  :
    m_isDefined(a_isDefined),
    m_params(a_params)
  {
  }


  /// Apply BC
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (!m_isDefined)
    {
      MayDay::Error("undefined BasicRefluxCorrBCFunction object");
    }

    RefCountedPtr<ConstValueFunction>
    zeroBC(new ConstValueFunction(0, a_state.nComp()));

    //    int order = 1;

    const Box& domainBox = a_domain.domainBox();

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (!a_domain.isPeriodic(idir))
      {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();
          //                                    int isign = sign(side);
          if (a_valid.sideEnd(side)[idir] ==
              domainBox.sideEnd(side)[idir])
          {
            NeumBC(a_state,
                   a_valid,
                   a_dx,
                   a_homogeneous,
                   BCValueHolder(zeroBC),
                   idir,
                   side);

          }
        }
      }
    }

  }
};


/// Boundary condition for the freestream correction
/**
 * sets normal gradient to be 0
 */
class FreestreamCorrBCFunction: public BCFunction
{
public:

  /// Is object defined?
  bool m_isDefined;
  /// Physical parameters for the problem
  MushyLayerParams m_params;

  /// Default constructor
  FreestreamCorrBCFunction()
  :
    m_isDefined(false)
  {
  }
  /// Full constructor
  FreestreamCorrBCFunction(bool a_isDefined,
                           MushyLayerParams a_params)
  :
    m_isDefined(a_isDefined),
    m_params(a_params)
  {
  }

  /// Apply BC
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_isDefined)
    {
      const Box& domainBox = a_domain.domainBox();
      //			PhysBCUtil bcInfo(m_params); // sets BCs from ParmParse table
      RefCountedPtr<ConstValueFunction>
      zeroFunc(new ConstValueFunction(0.0, a_state.nComp()));

      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (!a_domain.isPeriodic(idir))
        {
          SideIterator sit;
          for (sit.reset(); sit.ok(); ++sit)
          {
            Side::LoHiSide side = sit();
            if (a_valid.sideEnd(side)[idir] ==
                domainBox.sideEnd(side)[idir])
            {
              //							int bctype = bcInfo.getBC(idir, side);
              int bctype = m_params.getVelBCType(idir, side);
              switch(bctype)
              {
                //                case PhysBCUtil::Inflow :
                //                case PhysBCUtil::Outflow :
                //                case PhysBCUtil::VelInflowOutflow :
                //                case PhysBCUtil::OutflowNormal :
                case PhysBCUtil::SolidWall :
                case PhysBCUtil::noShear :
                  //                case PhysBCUtil::VelInflowPlume :
                case PhysBCUtil::Symmetry :
                {
                  NeumBC(a_state, a_valid, a_dx,
                         a_homogeneous,
                         BCValueHolder(zeroFunc),
                         idir, side);
                  break;
                }
                // For boundaries which allow normal flow, I don't think the above BCs are appropriate
                case PhysBCUtil::Inflow :
                case PhysBCUtil::Outflow :
                case PhysBCUtil::VelInflowOutflow :
                case PhysBCUtil::OutflowNormal :
                case PhysBCUtil::VelInflowPlume :
                case PhysBCUtil::OutflowPressureGrad :
                case PhysBCUtil::PressureHead :
                {
                  DiriBC(a_state, a_valid, a_dx,
                         a_homogeneous,
                         BCValueHolder(zeroFunc),
                         idir, side);
                  break;
                }
                default:
                {
                  MayDay::Error("FreestreamCorrBCFunction - unknown BC type");
                }
              } // end switch
            } // end condition of matching coordinates
          } // end loop over both sides
        } // end if not periodic in this direction
      } // end loop over all directions
    }
    else
    {
      MayDay::Error("undefined FreestreamCorrBCFunction object");
    }
  }
};



// ---------------------------------------------------------------
/// Boundary condition for the pressure gradient
/**
 * Extrapolation on all sides
 */
class BasicGradPressureBCFunction: public BCFunction
{
public:

  /// Is object defined?
  bool m_isDefined;
  /// Physical parameters for the problem
  MushyLayerParams m_params;

  /// Default constructor
  BasicGradPressureBCFunction()
  :
    m_isDefined(false)
  {
  }
  /// Full constructor
  BasicGradPressureBCFunction(bool a_isDefined,
                              MushyLayerParams a_params)
  :
    m_isDefined(a_isDefined),
    m_params(a_params)
  {
  }

  /// Apply BC
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_isDefined)
    {
      const Box& domainBox = a_domain.domainBox();

      RefCountedPtr<ConstValueFunction>
      zeroFunc(new ConstValueFunction(0.0, a_state.nComp()));

      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (!a_domain.isPeriodic(idir))
        {
          SideIterator sit;
          for (sit.reset(); sit.ok(); ++sit)
          {
            Side::LoHiSide side = sit();
            if (a_valid.sideEnd(side)[idir] ==
                domainBox.sideEnd(side)[idir])
            {

              int bctype = m_params.getVelBCType(idir, side);
              switch(bctype)
              {
                case PhysBCUtil::SolidWall :
                case PhysBCUtil::Inflow :
                case PhysBCUtil::VelInflowOutflow :
                case PhysBCUtil::noShear :
                case PhysBCUtil::VelInflowPlume :
                case PhysBCUtil::Symmetry :
                case PhysBCUtil::PressureHead :
                  //                case PhysBCUtil::Outflow :
                {
                  int order = 2;// equivalent of HOExtrapBC. This is crucial to get projection right.
                  for (int comp=0; comp<a_state.nComp(); comp++)
                  {
                    ExtraBC(a_state, a_valid,
                            idir, side, order, comp);
                  }
                  break;
                }
                case PhysBCUtil::Outflow :
                case PhysBCUtil::OutflowNormal :
                case PhysBCUtil::OutflowPressureGrad :
                {
                  //Constant pressure BC
                  // Doesn't work with higher order?
                  int order = 1;
                  DiriBC(a_state, a_valid, a_dx,
                         a_homogeneous,
                         BCValueHolder(zeroFunc),
                         idir, side, order);
                  break;
                }


                default:
                {
                  MayDay::Error("BasicGradPressureBCFunction - unknown BC type");
                }
              } // end switch
            } // end condition of matching coordinates
          } // end loop over both sides
        } // end if not periodic in this direction
      } // end loop over all directions
    }
    else
    {
      MayDay::Error("undefined BasicGradPressureBCFunction object");
    }
  }
};

// ---------------------------------------------------------------
/// Extrapolation boundary condition
/**
 * Extrapolation on all sides
 */
class DomainExtrapBCFunction: public BCFunction
{
public:

  /// Is object defined?
  bool m_isDefined;
  /// Physical parameters for the problem
  MushyLayerParams m_params;

  /// Default constructor
  DomainExtrapBCFunction()
  :
    m_isDefined(false)
  {
  }
  /// Full constructor
  DomainExtrapBCFunction(bool a_isDefined,
                         MushyLayerParams a_params)
  :
    m_isDefined(a_isDefined),
    m_params(a_params)
  {
  }

  /// Apply BC
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_isDefined)
    {
      const Box& domainBox = a_domain.domainBox();

      RefCountedPtr<ConstValueFunction>
      zeroFunc(new ConstValueFunction(0.0, a_state.nComp()));

      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (!a_domain.isPeriodic(idir))
        {
          SideIterator sit;
          for (sit.reset(); sit.ok(); ++sit)
          {
            Side::LoHiSide side = sit();
            if (a_valid.sideEnd(side)[idir] ==
                domainBox.sideEnd(side)[idir])
            {

              int order = 2;// equivalent of HOExtrapBC. This is crucial to get projection right.

              for (int comp=0; comp<a_state.nComp(); comp++)
              {
                ExtraBC(a_state, a_valid,
                        idir, side, order, comp);
              }



            } // end condition of matching coordinates
          } // end loop over both sides
        } // end if not periodic in this direction
      } // end loop over all directions
    }
    else
    {
      MayDay::Error("undefined DomainExtrapBCFunction object");
    }
  }
};


/// Apply extrapolation boundary condition on all sides
/**
 *
 */
class BasicExtrapBCFunction: public BCFunction
{
public:

  /// Is object defined?
  bool m_isDefined;

  /// Physical parameters for the problem
  MushyLayerParams m_params;

  /// Number of ghost cells to fill
  Interval m_ghost;

  /// Component of velocity to apply BCs to
  int m_comp;

  /// Default constructor
  BasicExtrapBCFunction()
  :
    m_isDefined(false), m_comp(0)
  {
  }
  /// Full constructor
  BasicExtrapBCFunction(bool a_isDefined,
                        MushyLayerParams a_params,
                        Interval a_ghost = Interval(0,0), // 0->0: fill first ghost cells
                        int a_comp = 1)
  :
    m_isDefined(a_isDefined),
    m_params(a_params),
    m_ghost(a_ghost),
    m_comp(a_comp)
  {
  }

  /// Apply BC
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_isDefined)
    {
      const Box& domainBox = a_domain.domainBox();
      // PhysBCUtil bcInfo(m_params); // sets BCs from ParmParse table

      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (!a_domain.isPeriodic(idir))
        {
          SideIterator sit;
          for (sit.reset(); sit.ok(); ++sit)
          {
            Side::LoHiSide side = sit();
            if (a_valid.sideEnd(side)[idir] ==
                domainBox.sideEnd(side)[idir])
            {

              Box grownBox(a_valid);


              // First expand box in case we're doing outer ghost cells
              // Expand once for Interval(0,x), twice for interval(1,x) etc.
              grownBox.grow(m_ghost.begin());


              //Grow the box by one cell in the direction perpendicular to the direction we are extrapolating in,
              //So that we fill enough ghost cells to be able to correctly calculate face values

              for (int i=0; i<SpaceDim; i++)
              {
                if (i != idir)
                {
                  grownBox.grow(i, 1);
                }
              }



              int order = 1;

              for (int ghost_i = m_ghost.begin(); ghost_i <=m_ghost.end(); ghost_i++)
              {
                //								ExtraBC(a_state, grownBox,
                //										idir, side, order);
                for (int comp=0; comp<a_state.nComp(); comp++)
                {
                  ExtrapBC(  a_state, grownBox,  idir,   side, order, comp);
                }

                // Grow box in the direction we're extrapolating in so we can fill more ghost cells
                grownBox.grow(idir, 1);
              }


            } // end condition of matching coordinates
          } // end loop over both sides
        } // end if not periodic in this direction
      } // end loop over all directions
    }
    else
    {
      MayDay::Error("undefined BasicExtrapBCFunction object");
    }
  }
};



/// This class only fills cells on the interior of the domain!
class BasicExtrapInteriorFunction: public BCFunction
{
public:

  /// Is object defined?
  bool m_isDefined;
  /// Physical parameters for the problem
  MushyLayerParams m_params;
  /// Number of ghost cells to fill
  int m_ghost;

  /// Default constructor
  BasicExtrapInteriorFunction()
  :
    m_isDefined(false),
    m_ghost(0)
  {
  }
  /// Full constructor
  BasicExtrapInteriorFunction(bool a_isDefined,
                              MushyLayerParams a_params,
                              int a_ghost)
  :
    m_isDefined(a_isDefined),
    m_params(a_params),
    m_ghost(a_ghost)
  {
  }

  /// Apply BC
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_isDefined)
    {
      // hard code component for now - later may change this
      int comp = 0;

      const Box& domainBox = a_domain.domainBox();
      //      PhysBCUtil bcInfo(m_params); // sets BCs from ParmParse table

      for (int idir = 0; idir < SpaceDim; idir++)
      {

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();
          if ( (a_valid.sideEnd(side)[idir] < domainBox.sideEnd(side)[idir] && side == Side::Hi)
              or (a_valid.sideEnd(side)[idir] > domainBox.sideEnd(side)[idir] && side == Side::Lo) )
          {
            int order = 2; // equivalent of HOExtrapBC

            //Grow the box by one cell in the direction perpendicular to the direction we are extrapolating in,
            //So that we fill enough ghost cells to be able to correctly calculate face values
            Box grownBox(a_valid);
            for (int i=0; i<SpaceDim; i++)
            {
              if (i != idir)
              {
                grownBox.grow(i, 1);
              }
            }

            for (int ghost_i = 1; ghost_i <=m_ghost; ghost_i++)
            {
              //							ExtraBC(a_state, grownBox,
              //									idir, side, order);

              ExtrapBC(  a_state, grownBox,  idir,   side, order, comp);

              // Grow box in the direction we're extrapolating in so we can fill more ghost cells
              grownBox.grow(idir, 1);
            }


          } // end condition of matching coordinates
        } // end loop over both sides

      } // end loop over all directions
    }
    else
    {
      MayDay::Error("undefined BasicExtrapInteriorFunction object");
    }
  }
};



Tuple<BCHolder, SpaceDim>
PhysBCUtil::velExtrapBC(Interval ghostInterval) const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  for (int idir=0; idir<SpaceDim; idir++)
  {
    //    Interval intvl(idir, idir);

    bcVec[idir] = extrapVelFuncBC(idir, ghostInterval);
  }
  return bcVec;
}

//Tuple<BCHolder, SpaceDim>
//PhysBCUtil::velExtrapInterior(int n_ghost) const
//{
//  Tuple<BCHolder, SpaceDim> bcVec;
//  bool isHomogeneous = false;
//  for (int idir=0; idir<SpaceDim; idir++)
//  {
//    Interval intvl(idir, idir);
//
//    bcVec[idir] = extrapVelFuncInterior(isHomogeneous, idir,
//                                        intvl, n_ghost);
//  }
//  return bcVec;
//}

BCHolder PhysBCUtil::viscousRefluxBC(int a_dir, bool a_isViscous) const
{
  //  bool isViscous = true;
  bool isHomogeneous = true;
  Interval intvl(0, 0); // 13 nov 2007:  not sure about this
  return basicCCVelFuncBC(isHomogeneous, a_isViscous, a_dir,
                          intvl);
}

// ---------------------------------------------------------------
//todo - should implement more than just solid walls BCs, for the cases where we have either slip or outflow bcs
PhysIBC*
PhysBCUtil::advectionVelIBC() const
{
  VelIBC* newIBCPtr = new VelIBC;
  // default to solid wall BC's everywhere
  // (i.e. set normal velocities to 0)
  // if something different (like inflow) is desired, do this in the
  // derived BC server class.
  for (int idir=0; idir<SpaceDim; idir++)
  {
    SideIterator sit;
    for (sit.reset(); sit.ok(); ++sit)
    {
      Side::LoHiSide side = sit();
      newIBCPtr->setNormalWallVel(0.0, idir, side);
    }
  }

  return newIBCPtr;
}


// ---------------------------------------------------------------
Tuple<BCHolder, SpaceDim>
PhysBCUtil::uStarFuncBC(bool a_isViscous) const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  bool isHomogeneous = false;
  for (int idir=0; idir<SpaceDim; idir++)
  {
    Interval intvl(idir, idir);
    bcVec[idir] = basicCCVelFuncBC(isHomogeneous, a_isViscous, idir,
                                   intvl);
  }
  return bcVec;
}

Tuple<BCHolder, SpaceDim>
PhysBCUtil::edgeVelFuncBC(bool a_isViscous,
                          LevelData<FluxBox>* velocityBCVals) const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  bool isHomogeneous = false;
  for (int idir=0; idir<SpaceDim; idir++)
  {
    // data will have one component in each space dimension
    Interval intvl(0, 0);
    bcVec[idir] = basicECVelFuncBC(isHomogeneous, a_isViscous, idir,
                                   intvl, velocityBCVals);
  }
  return bcVec;
}



Tuple<BCHolder, SpaceDim>
PhysBCUtil::porosityFaceBC() const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  for (int idir=0; idir<SpaceDim; idir++)
  {
    // data will have one component in each space dimension
    Interval intvl(0, 0);
    bcVec[idir] = BasicPorosityFaceFuncBC(idir);
  }
  return bcVec;
}

Tuple<BCHolder, SpaceDim>
PhysBCUtil::permeabilityFaceBC() const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  for (int idir=0; idir<SpaceDim; idir++)
  {
    // data will have one component in each space dimension
    Interval intvl(0, 0);
    bcVec[idir] = BasicPermeabilityFaceFuncBC(idir);
  }
  return bcVec;
}

// ---------------------------------------------------------------
Tuple<BCHolder, SpaceDim>
PhysBCUtil::tracingVelFuncBC() const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  bool isHomogeneous = false;
  bool isViscous = false;
  for (int idir=0; idir<SpaceDim; idir++)
  {
    Interval intvl(idir, idir);
    bcVec[idir] = basicCCVelFuncBC(isHomogeneous, isViscous, idir,
                                   intvl);
  }
  return bcVec;
}


// ---------------------------------------------------------------
Tuple<BCHolder, SpaceDim>
PhysBCUtil::advectionVelFuncBC(bool a_isViscous) const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  bool isHomogeneous = false;
  for (int idir=0; idir<SpaceDim; idir++)
  {
    // data will have one component in each space dimension
    Interval intvl(0, 0);
    bcVec[idir] = basicECVelFuncBC(isHomogeneous, a_isViscous, idir,
                                   intvl);
  }
  return bcVec;
}


Tuple<BCHolder, SpaceDim>
PhysBCUtil::fluxExtrapBC() const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  for (int idir=0; idir<SpaceDim; idir++)
  {
    // data will have one component in each space dimension
    RefCountedPtr<BasicFluxExtrapBCFunction>  bc (new BasicFluxExtrapBCFunction(idir));
    bcVec[idir] = BCHolder(bc);
  }

  return bcVec;
}




BCHolder PhysBCUtil::basicECVelFuncBC(bool a_isHomogeneous,
                                      bool a_isViscous,
                                      int  a_comp,
                                      const Interval& a_interval,
                                      LevelData<FluxBox>* a_velocityBCVals) const
{
  RefCountedPtr<BasicECVelBCFunction>
  basicECVelBCFunction(new BasicECVelBCFunction(a_isHomogeneous,
                                                a_isViscous,
                                                a_comp,
                                                a_interval, m_params, a_velocityBCVals));

  return BCHolder(basicECVelBCFunction);
}


BCHolder PhysBCUtil::streamFunctionBC() const
{
  RefCountedPtr<StreamFunctionBC>
  bc(new StreamFunctionBC(true, m_params,
                          false, // not homogeneous
                          m_advVel, m_dx));

  return BCHolder(bc);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::BasicGradPressureFuncBC() const
{
  RefCountedPtr<BasicGradPressureBCFunction>
  basicGradPressureBCFunction(new BasicGradPressureBCFunction(true, m_params));

  return BCHolder(basicGradPressureBCFunction);
}


/// For subcycled code
BCHolder PhysBCUtil::LevelPressureFuncBC() const
{
  return BasicPressureFuncBC(false);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::gradPiFuncBC() const
{
  return BasicGradPressureFuncBC();
}

BCHolder PhysBCUtil::extrapFuncBC() const
{
  RefCountedPtr<DomainExtrapBCFunction>
  bc(new DomainExtrapBCFunction(true, m_params));

  return BCHolder(bc);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::gradESyncFuncBC() const
{
  return BasicGradPressureFuncBC();
}

BCHolder PhysBCUtil::SyncProjFuncBC() const
{
  return BasicPressureFuncBC(true);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::gradMacPressureFuncBC() const
{
  return BasicGradPressureFuncBC();
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::gradELambdaFuncBC() const
{
  return BasicGradPressureFuncBC();
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::FreestreamCorrFuncBC() const
{
  RefCountedPtr<FreestreamCorrBCFunction>
  freestreamCorrBCFunction(new FreestreamCorrBCFunction(true, m_params));

  return BCHolder(freestreamCorrBCFunction);
}



/// For subcycled code
BCHolder PhysBCUtil::BasicPressureFuncBC(bool a_isHomogeneous) const
{

  RefCountedPtr<BasicPressureBCFunctionSubcycled>
  basicPressureBCFunction(new BasicPressureBCFunctionSubcycled(true,m_params, a_isHomogeneous,
                                                               m_advVel, m_dx));

  return BCHolder(basicPressureBCFunction);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::velFuncBC(int a_idir, bool a_viscous, Interval interval) const
{
  bool isHomogeneous = false;
  return basicCCVelFuncBC(isHomogeneous, a_viscous, a_idir,
                          interval);
}



// ---------------------------------------------------------------
/// Be very careful with this function - it fills CF ghost cells
BCHolder PhysBCUtil::extrapolationFuncBC(int a_order) const
{
  RefCountedPtr<ExtrapolationBCFunction>
  viscousBCFunction(new ExtrapolationBCFunction(true, a_order));

  return BCHolder(viscousBCFunction);
}


BCHolder PhysBCUtil::ThetaLFuncBC(bool a_homogeneous, LevelData<FluxBox>* a_advVel) const
{

  Vector<Vector<int> > customLoBC, customHiBC; // First index refers to direction, second to component
  Vector< Vector<Real> > customLoBCVal, customHiBCVal; // first index refers to direction, second to component
  Vector<Real> plumeVal(1, m_params.ThetaLPlumeInflow);

  for (int dir =0; dir<SpaceDim; dir++)
  {
    customLoBC.push_back(Vector<int>(1, m_params.bcTypeLiquidConcentrationLo[dir]));
    customHiBC.push_back(Vector<int>(1, m_params.bcTypeLiquidConcentrationHi[dir]));

    customLoBCVal.push_back(Vector<Real>(1, m_params.bcValLiquidConcentrationLo[dir]));
    customHiBCVal.push_back(Vector<Real>(1, m_params.bcValLiquidConcentrationHi[dir]));

  }

  RefCountedPtr<AdvectDiffuseScalarBC>
  basicThetaLBCFunction(new AdvectDiffuseScalarBC(true,
                                                  m_params, a_homogeneous,
                                                  m_advVel, m_dx,
                                                  plumeVal,
                                                  Interval(0,0),
                                                  customLoBC, customHiBC,
                                                  customLoBCVal, customHiBCVal));

  return BCHolder(basicThetaLBCFunction);
}


BCHolder PhysBCUtil::TracerBC(bool a_homogeneous, LevelData<FluxBox>* a_advVel, Real a_boundaryVal) const
{

  Vector<Vector<int> > customLoBC, customHiBC; // First index refers to direction, second to component
  Vector< Vector<Real> > customLoBCVal, customHiBCVal; // first index refers to direction, second to component
  Vector<Real> plumeVal(1, 1.0);

  for (int dir =0; dir<SpaceDim; dir++)
  {
//    customLoBC.push_back(Vector<int>(1, m_params.bcTypeLiquidConcentrationLo[dir]));
//    customHiBC.push_back(Vector<int>(1, m_params.bcTypeLiquidConcentrationHi[dir]));

    customLoBCVal.push_back(Vector<Real>(1, a_boundaryVal));
    customHiBCVal.push_back(Vector<Real>(1, a_boundaryVal));

  }

  RefCountedPtr<AdvectDiffuseScalarBC>
  bc(new AdvectDiffuseScalarBC(true,
                                                  m_params, a_homogeneous,
                                                  m_advVel, m_dx,
                                                  plumeVal,
                                                  Interval(0,0),
                                                  customLoBC, customHiBC,
                                                  customLoBCVal, customHiBCVal));

  return BCHolder(bc);
}




BCHolder PhysBCUtil::BasicthetaFuncBC(bool a_homogeneous, LevelData<FluxBox>* a_advVel) const
{

  Vector<Vector<int> >   customLoBC, customHiBC;       // First index refers to direction, second to component
  Vector< Vector<Real> > customLoBCVal, customHiBCVal; // first index refers to direction, second to component

  Vector<Real> plumeVal(1, m_params.thetaPlumeInflow);

  for (int dir =0; dir<SpaceDim; dir++)
  {
    customLoBC.push_back(Vector<int>(1, m_params.bcTypeTemperatureLo[dir]));
    customHiBC.push_back(Vector<int>(1, m_params.bcTypeTemperatureHi[dir]));

    customLoBCVal.push_back(Vector<Real>(1, m_params.bcValTemperatureLo[dir]));
    customHiBCVal.push_back(Vector<Real>(1, m_params.bcValTemperatureHi[dir]));

  }

  RefCountedPtr<AdvectDiffuseScalarBC>
  basicthetaBCFunction(new AdvectDiffuseScalarBC(true,
                                                 m_params, a_homogeneous, m_advVel,
                                                 m_dx,
                                                 plumeVal,
                                                 Interval(0,0),
                                                 customLoBC, customHiBC,
                                                 customLoBCVal, customHiBCVal));

  return BCHolder(basicthetaBCFunction);
}

void PhysBCUtil::applyFrameAdvectionBC(IntVect& bcTypeHi, IntVect& bcTypeLo) const
{

  if (m_params.nonDimVel > 0.0)
  {
    bcTypeHi[SpaceDim-1] = AdvectIBC::m_inflowOutflow;
  }

  if (m_params.nonDimVel < 0.0)
  {
    bcTypeLo[SpaceDim-1] = AdvectIBC::m_inflowOutflow;
  }

}

PhysIBC*
PhysBCUtil::scalarTraceH_IBC() const
{
  AdvectIBC* newIBCPtr = new AdvectIBC(1);
  RealVect bcValHi, bcValLo;
  IntVect bcTypeHi, bcTypeLo;

  for (int dir=0; dir<SpaceDim; dir++)
  {
    bcTypeHi[dir] = convertBCType(m_params.bcTypeEnthalpyHi[dir]);
    bcTypeLo[dir] = convertBCType(m_params.bcTypeEnthalpyLo[dir]);

    bcValLo[dir] = m_params.bcValEnthalpyLo[dir];
    bcValHi[dir] = m_params.bcValEnthalpyHi[dir];
  }

  applyFrameAdvectionBC(bcTypeHi, bcTypeLo);


  newIBCPtr->setBoundaryValues(bcValLo, bcTypeLo,  Side::Lo,
                               0); // component 0 - enthalpy
  newIBCPtr->setBoundaryValues(bcValHi, bcTypeHi,  Side::Hi,
                               0); // component 0 - enthalpy

  newIBCPtr->setAdvVel(m_advVel);


  Vector<Real> plumeVals;
  plumeVals.push_back(m_params.HPlumeInflow);
  //    plumeVals.push_back(m_params.ThetaPlumeInflow);
  newIBCPtr->setPlume(plumeVals, m_params.plumeBounds);


  return newIBCPtr;

}

PhysIBC*
PhysBCUtil::scalarTrace_IBC(Vector<int> a_bcTypeHi, Vector<int> a_bcTypeLo, RealVect a_bcValHi,RealVect a_bcValLo, Real a_plumeInflow) const
{
  AdvectIBC* newIBCPtr = new AdvectIBC(1);
  RealVect bcValHi, bcValLo;
  IntVect bcTypeHi, bcTypeLo;

  for (int dir=0; dir<SpaceDim; dir++)
  {
    bcTypeHi[dir] = convertBCType(a_bcTypeHi[dir]);
    bcTypeLo[dir] = convertBCType(a_bcTypeLo[dir]);

    bcValLo[dir] = a_bcValLo[dir];
    bcValHi[dir] = a_bcValHi[dir];
  }

  applyFrameAdvectionBC(bcTypeHi, bcTypeLo);

  newIBCPtr->setBoundaryValues(bcValLo, bcTypeLo,  Side::Lo,
                               0);
  newIBCPtr->setBoundaryValues(bcValHi, bcTypeHi,  Side::Hi,
                               0);

  newIBCPtr->setAdvVel(m_advVel);

  Vector<Real> plumeVals;
  plumeVals.push_back(a_plumeInflow);
  newIBCPtr->setPlume(plumeVals, m_params.plumeBounds);

  return newIBCPtr;
}

PhysIBC*
PhysBCUtil::scalarTraceC_IBC() const
{
  AdvectIBC* newIBCPtr = new AdvectIBC(1);
  RealVect bcValHi, bcValLo;
  IntVect bcTypeHi, bcTypeLo;

  for (int dir=0; dir<SpaceDim; dir++)
  {
    bcTypeHi[dir] = convertBCType(m_params.bcTypeBulkConcentrationHi[dir]);
    bcTypeLo[dir] = convertBCType(m_params.bcTypeBulkConcentrationLo[dir]);

    bcValLo[dir] = m_params.bcValBulkConcentrationLo[dir];
    bcValHi[dir] = m_params.bcValBulkConcentrationHi[dir];
  }

  applyFrameAdvectionBC(bcTypeHi, bcTypeLo);

  newIBCPtr->setBoundaryValues(bcValLo, bcTypeLo,  Side::Lo,
                               0);
  newIBCPtr->setBoundaryValues(bcValHi, bcTypeHi,  Side::Hi,
                               0);

  newIBCPtr->setAdvVel(m_advVel);


  Vector<Real> plumeVals;
  //   plumeVals.push_back(m_params.HPlumeInflow);
  plumeVals.push_back(m_params.ThetaPlumeInflow);
  newIBCPtr->setPlume(plumeVals, m_params.plumeBounds);


  return newIBCPtr;

}

PhysIBC*
PhysBCUtil::scalarTraceHC_IBC() const
{
  AdvectIBC* newIBCPtr = new AdvectIBC(2);

  // New N dimensional code. Just do enthalpy then bulk concentration
  RealVect bcValHi, bcValLo;
  IntVect bcTypeHi, bcTypeLo;

  for (int dir=0; dir<SpaceDim; dir++)
  {
    bcTypeHi[dir] = convertBCType(m_params.bcTypeEnthalpyHi[dir]);
    bcTypeLo[dir] = convertBCType(m_params.bcTypeEnthalpyLo[dir]);

    bcValLo[dir] = m_params.bcValEnthalpyLo[dir];
    bcValHi[dir] = m_params.bcValEnthalpyHi[dir];
  }

  applyFrameAdvectionBC(bcTypeHi, bcTypeLo);


  newIBCPtr->setBoundaryValues(bcValLo, bcTypeLo,  Side::Lo,
                               0); // component 0 - enthalpy
  newIBCPtr->setBoundaryValues(bcValHi, bcTypeHi,  Side::Hi,
                               0); // component 0 - enthalpy


  // Now bulk concentration
  for (int dir=0; dir<SpaceDim; dir++)
  {
    bcTypeHi[dir] = convertBCType(m_params.bcTypeBulkConcentrationHi[dir]);
    bcTypeLo[dir] = convertBCType(m_params.bcTypeBulkConcentrationLo[dir]);

    bcValLo[dir] = m_params.bcValBulkConcentrationLo[dir];
    bcValHi[dir] = m_params.bcValBulkConcentrationHi[dir];
  }

  applyFrameAdvectionBC(bcTypeHi, bcTypeLo);

  newIBCPtr->setBoundaryValues(bcValLo, bcTypeLo,  Side::Lo,
                               1); // component 1 - bulk concentration
  newIBCPtr->setBoundaryValues(bcValHi, bcTypeHi,  Side::Hi,
                               1); // component 1 - bulk concentration

  newIBCPtr->setAdvVel(m_advVel);


  Vector<Real> plumeVals;
  plumeVals.push_back(m_params.HPlumeInflow);
  plumeVals.push_back(m_params.ThetaPlumeInflow);
  newIBCPtr->setPlume(plumeVals, m_params.plumeBounds);

  return newIBCPtr;
}

PhysIBC*
PhysBCUtil::scalarTraceTSl_IBC() const
{
  AdvectIBC* newIBCPtr = new AdvectIBC(2);

  // New N dimensional code. Just do temperature then liquid concentration
  RealVect bcValHi, bcValLo;
  IntVect bcTypeHi, bcTypeLo;


  for (int dir=0; dir<SpaceDim; dir++)
  {
    bcTypeHi[dir] = convertBCType(m_params.bcTypeTemperatureHi[dir]);
    bcTypeLo[dir] = convertBCType(m_params.bcTypeTemperatureLo[dir]);

    bcValLo[dir] = m_params.bcValTemperatureLo[dir];
    bcValHi[dir] = m_params.bcValTemperatureHi[dir];
  }

  //  if (m_params.nonDimVel != 0.0)
  //  {
  //    bcTypeHi[SpaceDim-1] = AdvectIBC::m_inflowOutflow;
  //  }
  applyFrameAdvectionBC(bcTypeHi, bcTypeLo);

  newIBCPtr->setBoundaryValues(bcValLo, bcTypeLo,  Side::Lo,
                               0); // component 0 - temperature
  newIBCPtr->setBoundaryValues(bcValHi, bcTypeHi,  Side::Hi,
                               0); // component 0 - temperature


  // Now liquid concentration
  for (int dir=0; dir<SpaceDim; dir++)
  {
    bcTypeHi[dir] = convertBCType(m_params.bcTypeLiquidConcentrationHi[dir]);
    bcTypeLo[dir] = convertBCType(m_params.bcTypeLiquidConcentrationLo[dir]);
    bcValLo[dir] = m_params.bcValLiquidConcentrationLo[dir];
    bcValHi[dir] = m_params.bcValLiquidConcentrationHi[dir];
  }

  //  if (m_params.nonDimVel != 0.0)
  //  {
  //    bcTypeHi[SpaceDim-1] = AdvectIBC::m_inflowOutflow;
  //  }
  applyFrameAdvectionBC(bcTypeHi, bcTypeLo);

  newIBCPtr->setBoundaryValues(bcValLo, bcTypeLo,  Side::Lo,
                               1); // component 1 - Liquid concentration
  newIBCPtr->setBoundaryValues(bcValHi, bcTypeHi,  Side::Hi,
                               1); // component 1 - Liquid concentration

  newIBCPtr->setAdvVel(m_advVel);

  Vector<Real> plumeVals;
  plumeVals.push_back(m_params.thetaPlumeInflow);
  plumeVals.push_back(m_params.ThetaLPlumeInflow);
  newIBCPtr->setPlume(plumeVals, m_params.plumeBounds);



  return newIBCPtr;
}


// -------------------------------------------------------
/// this is a BC object used in the PatchGodunov stuff
PhysIBC*
PhysBCUtil::scalarTraceIBC(RealVect bcValLo, RealVect bcValHi,
                           IntVect bcTypeLo, IntVect bcTypeHi,
                           Real plumeVal, Vector<Real> plumeBounds) const
{
  AdvectIBC* newIBCPtr = new AdvectIBC;

  newIBCPtr->setBoundaryValues(bcValLo, bcTypeLo, Side::Lo);
  newIBCPtr->setBoundaryValues(bcValHi, bcTypeHi, Side::Hi);

  newIBCPtr->setAdvVel(m_advVel);

  Vector<Real> plumeVals;
  plumeVals.push_back(plumeVal);
  newIBCPtr->setPlume(plumeVals, m_params.plumeBounds);




  return newIBCPtr;
}


BCHolder PhysBCUtil::lambdaBC(bool a_homogeneous, bool a_scaleWithPorosity) const
{

  Vector<Vector<int> > customLoBC, customHiBC; // First index refers to direction, second to component
  Vector< Vector<Real> > customLoBCVal, customHiBCVal; // first index refers to direction, second to component

  customLoBC.resize(SpaceDim); customHiBC.resize(SpaceDim);
  customLoBCVal.resize(SpaceDim); customHiBCVal.resize(SpaceDim);

  for (int dir =0; dir<SpaceDim; dir++)
  {
    customLoBC[dir].resize(1, PhysBCUtil::Dirichlet);
    customHiBC[dir].resize(1, PhysBCUtil::Dirichlet);

    customLoBCVal[dir].resize(1, 1.0);
    customHiBCVal[dir].resize(1, 1.0);

    if (a_scaleWithPorosity)
    {
      customLoBCVal[dir][0] /= m_params.bcValPorosityLo[dir];
      customHiBCVal[dir][0] /= m_params.bcValPorosityHi[dir];
    }

  }

  RefCountedPtr<AdvectDiffuseScalarBC>
  basicScalarBCFunction(new AdvectDiffuseScalarBC(true,
                                                  m_params, a_homogeneous,
                                                  m_advVel, m_dx,
                                                  1.0, Interval(0,0),
                                                  customLoBC, customHiBC,
                                                  customLoBCVal,customHiBCVal
  ));

  return BCHolder(basicScalarBCFunction);
}


BCHolder PhysBCUtil::noFluxBC() const
{
  RefCountedPtr<NoFluxBCFunction>
  bc(new NoFluxBCFunction(true, m_params));

  return BCHolder(bc);
}



BCHolder PhysBCUtil::BasicPorosityFuncBC(bool a_homogeneous) const
{
  Vector<Vector<int> > customLoBC, customHiBC; // First index refers to direction, second to component
  Vector< Vector<Real> > customLoBCVal, customHiBCVal; // first index refers to direction, second to component
  Vector<Real> plumeVals;

  plumeVals.push_back(m_params.porosityPlume);

  for (int dir =0; dir<SpaceDim; dir++)
  {
    customLoBC.push_back(Vector<int>(1, m_params.bcTypePorosityLo[dir]));
    customHiBC.push_back(Vector<int>(1, m_params.bcTypePorosityHi[dir]));

    customLoBCVal.push_back(Vector<Real>(1, m_params.bcValPorosityLo[dir]));
    customHiBCVal.push_back(Vector<Real>(1, m_params.bcValPorosityHi[dir]));
  }

  RefCountedPtr<BasicPorosityPermeabilityBCFunction>
  basicPorosityBCFunction(new BasicPorosityPermeabilityBCFunction(true,
                                                                  m_params, a_homogeneous,
                                                                  m_advVel, m_dx,
                                                                  plumeVals,
                                                                  Interval(0,0),
                                                                  customLoBC, customHiBC,
                                                                  customLoBCVal, customHiBCVal));

  return BCHolder(basicPorosityBCFunction);
}

BCHolder PhysBCUtil::BasicPorosityFaceFuncBC(int a_comp, bool a_homogeneous) const
{
  Vector<Vector<int> > customLoBC, customHiBC; // First index refers to direction, second to component
  Vector< Vector<Real> > customLoBCVal, customHiBCVal; // first index refers to direction, second to component

  for (int dir =0; dir<SpaceDim; dir++)
  {
    customLoBC.push_back(Vector<int>(1, m_params.bcTypePorosityLo[dir]));
    customHiBC.push_back(Vector<int>(1, m_params.bcTypePorosityHi[dir]));

    customLoBCVal.push_back(Vector<Real>(1, m_params.bcValPorosityLo[dir]));
    customHiBCVal.push_back(Vector<Real>(1, m_params.bcValPorosityHi[dir]));
  }


  RefCountedPtr<BasicPorosityPermeabilityFaceBCFunction>
  basicPorosityBCFunction(new BasicPorosityPermeabilityFaceBCFunction(a_homogeneous,
                                                                      a_comp, m_params,
                                                                      m_params.porosityPlume,
                                                                      customLoBC, customHiBC,
                                                                      customLoBCVal, customHiBCVal));

  return BCHolder(basicPorosityBCFunction);
}

BCHolder PhysBCUtil::BasicPermeabilityFaceFuncBC(int a_comp, bool a_homogeneous) const
{
  Vector<Vector<int> > customLoBC, customHiBC; // First index refers to direction, second to component
  Vector< Vector<Real> > customLoBCVal, customHiBCVal; // first index refers to direction, second to component

  for (int dir =0; dir<SpaceDim; dir++)
  {
    customLoBC.push_back(Vector<int>(1, m_params.bcTypePermeabilityLo[dir]));
    customHiBC.push_back(Vector<int>(1, m_params.bcTypePermeabilityHi[dir]));

    customLoBCVal.push_back(Vector<Real>(1, m_params.bcValPermeabilityLo[dir]));
    customHiBCVal.push_back(Vector<Real>(1, m_params.bcValPermeabilityHi[dir]));

  }

  RefCountedPtr<BasicPorosityPermeabilityFaceBCFunction> basicPermeabilityFaceBCFunction
  (new BasicPorosityPermeabilityFaceBCFunction(a_homogeneous,
                                               a_comp, m_params,
                                               m_params.permeabilityPlume,
                                               customLoBC, customHiBC,
                                               customLoBCVal, customHiBCVal));

  return BCHolder(basicPermeabilityFaceBCFunction);
}

BCHolder PhysBCUtil::BasicPermeabilityFuncBC(bool a_homogeneous) const
{

  Vector<Vector<int> > customLoBC, customHiBC; // First index refers to direction, second to component
  Vector< Vector<Real> > customLoBCVal, customHiBCVal; // first index refers to direction, second to component

  Vector<Real> plumeVal(1, m_params.permeabilityPlume);


  for (int dir =0; dir<SpaceDim; dir++)
  {
    customLoBC.push_back(Vector<int>(1, m_params.bcTypePermeabilityLo[dir]));
    customHiBC.push_back(Vector<int>(1, m_params.bcTypePermeabilityHi[dir]));

    customLoBCVal.push_back(Vector<Real>(1, m_params.bcValPermeabilityLo[dir]));
    customHiBCVal.push_back(Vector<Real>(1, m_params.bcValPermeabilityHi[dir]));

  }

  RefCountedPtr<BasicPorosityPermeabilityBCFunction>
  basicPermeabilityBCFunction(new BasicPorosityPermeabilityBCFunction(true,
                                                                      m_params, a_homogeneous,
                                                                      m_advVel, m_dx,
                                                                      plumeVal,
                                                                      Interval(0,0),
                                                                      customLoBC, customHiBC,
                                                                      customLoBCVal, customHiBCVal
  ));

  return BCHolder(basicPermeabilityBCFunction);
}


BCHolder PhysBCUtil::BasicEnthalpyFuncBC(bool a_homogeneous,
                                         LevelData<FluxBox>* a_advVel) const
{
  Vector<Vector<int> > customLoBC, customHiBC; // First index refers to direction, second to component
  Vector< Vector<Real> > customLoBCVal, customHiBCVal; // first index refers to direction, second to component

  Vector<Real> plumeVal;
  plumeVal.push_back(m_params.HPlumeInflow);

  for (int dir =0; dir<SpaceDim; dir++)
  {
    customLoBC.push_back(Vector<int>(1, m_params.bcTypeEnthalpyLo[dir]));
    customHiBC.push_back(Vector<int>(1, m_params.bcTypeEnthalpyHi[dir]));

    customLoBCVal.push_back(Vector<Real>(1, m_params.bcValEnthalpyLo[dir]));
    customHiBCVal.push_back(Vector<Real>(1, m_params.bcValEnthalpyHi[dir]));

  }

  RefCountedPtr<AdvectDiffuseScalarBC>
  basicEnthalpyBCFunction(new AdvectDiffuseScalarBC(true,
                                                    m_params, a_homogeneous, m_advVel,
                                                    m_dx,
                                                    plumeVal,
                                                    Interval(0,0),
                                                    customLoBC, customHiBC,
                                                    customLoBCVal, customHiBCVal ));

  return BCHolder(basicEnthalpyBCFunction);
}
BCHolder PhysBCUtil::BasicLiquidusBC(bool a_homogeneous,
                                     LevelData<FluxBox>* a_advVel) const
{
  Vector<Vector<int> > customLoBC, customHiBC; // First index refers to direction, second to component
  Vector< Vector<Real> > customLoBCVal, customHiBCVal; // first index refers to direction, second to component

  Vector<Real> plumeVal;
  plumeVal.push_back(m_params.HLiquidusPlume);

  for (int dir =0; dir<SpaceDim; dir++)
  {
    customLoBC.push_back(Vector<int>(1, m_params.bcTypeEnthalpyLo[dir]));
    customHiBC.push_back(Vector<int>(1, m_params.bcTypeEnthalpyHi[dir]));

    customLoBCVal.push_back(Vector<Real>(1, m_params.bcValLiquidusLo[dir]));
    customHiBCVal.push_back(Vector<Real>(1, m_params.bcValLiquidusHi[dir]));

  }

  RefCountedPtr<AdvectDiffuseScalarBC>
  bcFunction(new AdvectDiffuseScalarBC(true,
                                       m_params, a_homogeneous, m_advVel,
                                       m_dx,
                                       plumeVal,
                                       Interval(0,0),
                                       customLoBC, customHiBC,
                                       customLoBCVal, customHiBCVal ));

  return BCHolder(bcFunction);
}

BCHolder PhysBCUtil::BasicEutecticBC(bool a_homogeneous,
                                     LevelData<FluxBox>* a_advVel) const
{
  Vector<Vector<int> > customLoBC, customHiBC; // First index refers to direction, second to component
  Vector< Vector<Real> > customLoBCVal, customHiBCVal; // first index refers to direction, second to component

  Vector<Real> plumeVal;
  plumeVal.push_back(m_params.HEutecticPlume);

  for (int dir =0; dir<SpaceDim; dir++)
  {
    customLoBC.push_back(Vector<int>(1, m_params.bcTypeEnthalpyLo[dir]));
    customHiBC.push_back(Vector<int>(1, m_params.bcTypeEnthalpyHi[dir]));

    customLoBCVal.push_back(Vector<Real>(1, m_params.bcValEutecticLo[dir]));
    customHiBCVal.push_back(Vector<Real>(1, m_params.bcValEutecticHi[dir]));

  }

  RefCountedPtr<AdvectDiffuseScalarBC>
  bcFunction(new AdvectDiffuseScalarBC(true,
                                       m_params, a_homogeneous, m_advVel,
                                       m_dx,
                                       plumeVal,
                                       Interval(0,0),
                                       customLoBC, customHiBC,
                                       customLoBCVal, customHiBCVal ));

  return BCHolder(bcFunction);
}

BCHolder PhysBCUtil::BasicSolidusBC(bool a_homogeneous,
                                    LevelData<FluxBox>* a_advVel) const
{
  Vector<Vector<int> > customLoBC, customHiBC; // First index refers to direction, second to component
  Vector< Vector<Real> > customLoBCVal, customHiBCVal; // first index refers to direction, second to component

  Vector<Real> plumeVal;
  plumeVal.push_back(m_params.HSolidusPlume);

  for (int dir =0; dir<SpaceDim; dir++)
  {
    customLoBC.push_back(Vector<int>(1, m_params.bcTypeEnthalpyLo[dir]));
    customHiBC.push_back(Vector<int>(1, m_params.bcTypeEnthalpyHi[dir]));

    customLoBCVal.push_back(Vector<Real>(1, m_params.bcValSolidusLo[dir]));
    customHiBCVal.push_back(Vector<Real>(1, m_params.bcValSolidusHi[dir]));

  }

  RefCountedPtr<AdvectDiffuseScalarBC>
  bcFunction(new AdvectDiffuseScalarBC(true,
                                       m_params, a_homogeneous, m_advVel,
                                       m_dx,
                                       plumeVal,
                                       Interval(0,0),
                                       customLoBC, customHiBC,
                                       customLoBCVal, customHiBCVal ));

  return BCHolder(bcFunction);
}



BCHolder PhysBCUtil::temperatureLiquidSalinityBC(bool a_homogeneous) const
{
  //temperature BC for comp 0, liquid salinity BC for comp 1
  Vector<Real> topVals, bottomVals, plumeVals;
  topVals.resize(2); bottomVals.resize(2); plumeVals.resize(2);
  //  topVals[0] = m_params.thetaTop;
  //  bottomVals[0] = m_params.thetaBottom;
  plumeVals[0] = m_params.thetaPlumeInflow;

  //  topVals[1] = m_params.ThetaLTop;
  //  bottomVals[1] = m_params.ThetaLBottom;
  plumeVals[1] = m_params.ThetaLPlumeInflow;

  Vector<Vector<int> > bcLo, bcHi;
  Vector<Vector<Real> > bcValLo, bcValHi;

  // Fill with defaults
  bcLo.resize(SpaceDim); bcHi.resize(SpaceDim);
  bcValLo.resize(SpaceDim); bcValHi.resize(SpaceDim);

  for (int dir = 0; dir<SpaceDim; dir++)
  {
    bcLo[dir].resize(2); bcHi[dir].resize(2);
    bcValLo[dir].resize(2); bcValHi[dir].resize(2);

    bcLo[dir][0] = m_params.bcTypeTemperatureLo[dir];
    bcHi[dir][0] = m_params.bcTypeTemperatureHi[dir];

    bcLo[dir][1] = m_params.bcTypeLiquidConcentrationLo[dir];
    bcHi[dir][1] = m_params.bcTypeLiquidConcentrationHi[dir];

    bcValLo[dir][0] = m_params.bcValTemperatureLo[dir];
    bcValHi[dir][0] = m_params.bcValTemperatureHi[dir];

    bcValLo[dir][1] = m_params.bcValLiquidConcentrationLo[dir];
    bcValHi[dir][1] = m_params.bcValLiquidConcentrationHi[dir];
  }


  RefCountedPtr<AdvectDiffuseScalarBC>
  bcFunc(new AdvectDiffuseScalarBC(true,
                                   m_params, a_homogeneous, m_advVel,
                                   m_dx,
                                   plumeVals,
                                   Interval(0, 1),
                                   bcLo,
                                   bcHi,
                                   bcValLo,
                                   bcValHi));

  return BCHolder(bcFunc);
}




BCHolder PhysBCUtil::enthalpySalinityBC(bool a_homogeneous) const
{
  //temperature BC for comp 0, liquid salinity BC for comp 1
  Vector<Real> topVals, bottomVals, plumeVals;
  topVals.resize(2); bottomVals.resize(2); plumeVals.resize(2);

  plumeVals[0] = m_params.HPlumeInflow;
  plumeVals[1] = m_params.ThetaPlumeInflow;

  Vector<Vector<int> > bcLo, bcHi;
  Vector<Vector<Real> > bcValLo, bcValHi;

  // Fill with defaults
  bcLo.resize(SpaceDim); bcHi.resize(SpaceDim);
  bcValLo.resize(SpaceDim); bcValHi.resize(SpaceDim);

  for (int dir = 0; dir<SpaceDim; dir++)
  {
    bcLo[dir].resize(2); bcHi[dir].resize(2);
    bcValLo[dir].resize(2); bcValHi[dir].resize(2);

    bcLo[dir][0] = m_params.bcTypeEnthalpyLo[dir];
    bcHi[dir][0] = m_params.bcTypeEnthalpyHi[dir];

    bcLo[dir][1] = m_params.bcTypeBulkConcentrationLo[dir];
    bcHi[dir][1] = m_params.bcTypeBulkConcentrationHi[dir];

    bcValLo[dir][0] = m_params.bcValEnthalpyLo[dir];
    bcValHi[dir][0] = m_params.bcValEnthalpyHi[dir];

    bcValLo[dir][1] = m_params.bcValBulkConcentrationLo[dir];
    bcValHi[dir][1] = m_params.bcValBulkConcentrationHi[dir];
  }


  RefCountedPtr<AdvectDiffuseScalarBC>
  bcFunc(new AdvectDiffuseScalarBC(true,
                                   m_params, a_homogeneous, m_advVel,
                                   m_dx,
                                   plumeVals,
                                   Interval(0, 1),
                                   bcLo, bcHi,
                                   bcValLo, bcValHi) );

  return BCHolder(bcFunc);
}

BCHolder PhysBCUtil::BasicScalarFuncBC(Real a_topVal, Real a_bottomVal, Real a_plumeVal, bool a_homogeneous) const
{
  Vector<Vector<int> > bcLo, bcHi;
  Vector<Vector<Real> > bcValLo, bcValHi;
  Vector<Real> plumeVal(1, a_plumeVal);

  int comp = 0; // just one component

  // Fill with defaults
  bcLo.resize(SpaceDim); bcHi.resize(SpaceDim);
  bcValLo.resize(SpaceDim, 0.0); bcValHi.resize(SpaceDim);

  for (int dir = 0; dir<SpaceDim; dir++)
  {
    bcLo[dir].resize(1); bcHi[dir].resize(1);
    bcValLo[dir].resize(1); bcValHi[dir].resize(1);

    if (dir == SpaceDim - 1)
    {
      bcLo[dir][comp] = Dirichlet;
      bcHi[dir][comp] = Dirichlet;
      bcValLo[dir][comp] = a_bottomVal;
      bcValHi[dir][comp] = a_topVal;
    }
    else
    {
      bcLo[dir][comp] = Neumann;
      bcHi[dir][comp] = Neumann;
      bcValLo[dir][comp] = 0.0;
      bcValHi[dir][comp] = 0.0;
    }
  }

  RefCountedPtr<AdvectDiffuseScalarBC>
  basicScalarBCFunction(new AdvectDiffuseScalarBC(true,
                                                  m_params, a_homogeneous,
                                                  m_advVel, m_dx,
                                                  plumeVal,
                                                  Interval(0,0),
                                                  bcLo, bcHi,
                                                  bcValLo, bcValHi));

  return BCHolder(basicScalarBCFunction);
}


BCHolder PhysBCUtil::ThetaFuncBC(bool a_homogeneous, LevelData<FluxBox>* a_advVel) const
{

  Vector<Vector<int> > bcLo, bcHi;
  Vector<Vector<Real> > bcValLo, bcValHi;

  int comp = 0; // just one component

  // Fill with defaults
  bcLo.resize(SpaceDim); bcHi.resize(SpaceDim);
  bcValLo.resize(SpaceDim, 0.0); bcValHi.resize(SpaceDim);

  Vector<Real> plumeVal;
  plumeVal.push_back(m_params.ThetaPlumeInflow);

  for (int dir = 0; dir<SpaceDim; dir++)
  {
    bcLo[dir].resize(1); bcHi[dir].resize(1);
    bcValLo[dir].resize(1); bcValHi[dir].resize(1);

    bcLo[dir][comp] = m_params.bcTypeBulkConcentrationLo[dir];
    bcHi[dir][comp] = m_params.bcTypeBulkConcentrationHi[dir];

    bcValLo[dir][comp] = m_params.bcValBulkConcentrationLo[dir];
    bcValHi[dir][comp] = m_params.bcValBulkConcentrationHi[dir];

  }


  RefCountedPtr<AdvectDiffuseScalarBC>
  basicThetaBCFunction(new AdvectDiffuseScalarBC(true,
                                                 m_params, a_homogeneous, m_advVel,
                                                 m_dx,
                                                 plumeVal,
                                                 Interval(0,0),
                                                 bcLo, bcHi,
                                                 bcValLo, bcValHi));

  return BCHolder(basicThetaBCFunction);
}


BCHolder PhysBCUtil::extrapVelFuncBC( int  a_comp,
                                      Interval ghostInterval) const
{
  RefCountedPtr<BasicExtrapBCFunction>
  extrapVelFuncBC(new BasicExtrapBCFunction(true, m_params, ghostInterval, a_comp));

  return BCHolder(extrapVelFuncBC);
}

BCHolder PhysBCUtil::extrapVelFuncInterior(bool a_isHomogeneous,
                                           int  a_comp,
                                           const Interval& a_interval,
                                           int n_ghost) const
{
  RefCountedPtr<BasicExtrapInteriorFunction>
  extrapVelFuncInterior(new BasicExtrapInteriorFunction(true, m_params, n_ghost));

  return BCHolder(extrapVelFuncInterior);
}




// ---------------------------------------------------------------
BCHolder PhysBCUtil::basicCCVelFuncBC(bool a_isHomogeneous,
                                      bool a_isViscous,
                                      int  a_comp,
                                      const Interval& a_interval) const
{

  RefCountedPtr<BasicCCVelBCFunction> basicCCVelBCFunction(new BasicCCVelBCFunction(a_isHomogeneous,
                                                                                    a_isViscous,
                                                                                    a_comp,
                                                                                    a_interval,
                                                                                    m_params));

  return BCHolder(basicCCVelBCFunction);
}





void  PhysBCUtil::setAdvVel(LevelData<FluxBox>* a_advVel)
{
  m_advVel = a_advVel;
}
LevelData<FluxBox>*  PhysBCUtil::getAdvVel()
{
  return m_advVel;
}

PhysBCUtil::PhysBCUtil() : m_dx(-1),
    m_defined(false),
    m_advVel(nullptr),
    m_time(0)
{ }

// ---------------------------------------------------------------
PhysBCUtil::PhysBCUtil(MushyLayerParams a_params, Real a_dx)
{
  //  pout() << "PhysBCUtil::PhysBCUtil" << endl;
  m_advVel = nullptr;

  // initialize to bogus values
  for (int idir=0; idir<SpaceDim; idir++)
  {
    m_loBC[idir] = bogusBC;
    m_hiBC[idir] = bogusBC;
  }

  m_dx = a_dx;

  // now initialize BC's
  setBCs();

  m_params = a_params;

  m_defined = true;

}

// ---------------------------------------------------------------
PhysBCUtil::~PhysBCUtil()
{
}



// ---------------------------------------------------------------
PhysBCUtil&
PhysBCUtil::operator= (const PhysBCUtil& rhs)
{
  m_loBC = rhs.m_loBC;
  m_hiBC = rhs.m_hiBC;
  m_advVel = rhs.m_advVel;
  return *this;
}

// ---------------------------------------------------------------
PhysBCUtil::PhysBCUtil(const PhysBCUtil& rhs) :   m_loBC(rhs.m_loBC),
    m_hiBC(rhs.m_hiBC),
    m_dx(rhs.m_dx),
    m_time(rhs.m_time),
    m_defined(false),
    m_advVel(nullptr)
{ }

// ---------------------------------------------------------------
void
PhysBCUtil::setBCs()
{
  //	ParmParse ppBC("physBC");
  //
  //	vector<int> tempBC(SpaceDim);
  //
  //	ppBC.getarr("lo", tempBC, 0, SpaceDim);
  //	for (int idir=0; idir<SpaceDim; idir++)
  //	{
  //		m_loBC[idir] = tempBC[idir];
  //	}
  //
  //	ppBC.getarr("hi", tempBC, 0, SpaceDim);
  //	for (int idir=0; idir<SpaceDim; idir++)
  //	{
  //		m_hiBC[idir] = tempBC[idir];
  //	}

  for (int idir=0; idir<SpaceDim; idir++)
  {
    m_hiBC[idir] = SolidWall;
    m_loBC[idir] = SolidWall;
  }
}

// ---------------------------------------------------------------
PhysBCUtil*
PhysBCUtil::newPhysBCUtil() const
{
  PhysBCUtil* newBCPtr = new PhysBCUtil(m_params, m_dx);
  return newBCPtr;
}

// ----------------------------------------------------------
// basic utility functions
void
PhysBCUtil::Dx(const Real a_dx)
{
  m_dx = a_dx;
}

int PhysBCUtil::convertBCType(const int a_implicitBC) const
{
  int a_explicitBC=0;
  Real temp=0.0;
  convertBCType(a_implicitBC, temp, a_explicitBC, temp);
  return a_explicitBC;
}

void PhysBCUtil::convertBCType(const int a_implicitBC,  const Real a_implicitVal,
                               int& a_explicitBC, Real& a_explicitVal) const
{
  // Default BCs
  //  int ibcType = AdvectIBC::m_dirichlet;

  switch (a_implicitBC)
  {
    case PhysBCUtil::Dirichlet:
      a_explicitBC = AdvectIBC::m_dirichlet;
      a_explicitVal = a_implicitVal;
      break;

    case PhysBCUtil::Neumann:
      a_explicitBC = AdvectIBC::m_neumann;
      a_explicitVal = 0;
      break;


    case PhysBCUtil::InflowOutflow:
      a_explicitBC = AdvectIBC::m_inflowOutflow;
      a_explicitVal = a_implicitVal;
      break;

    case PhysBCUtil::OnlyInflow:
      a_explicitBC = AdvectIBC::m_plumeInflow;
      a_explicitVal = a_implicitVal;
      break;
  }

}

// ---------------------------------------------------------------
Real
PhysBCUtil::Dx() const
{
  return m_dx;
}

bool PhysBCUtil::isDefined()
{
  return m_defined;
}
// ---------------------------------------------------------------
void
PhysBCUtil::Time(const Real a_time)
{
  m_time = a_time;

  this->updateTimeDependentBCs();
}

void PhysBCUtil::updateTimeDependentBCs()
{


  // Update BCs if necessary

  if (m_params.m_timeDependentBC == MushyLayerParams::m_sinusoid)
  {

    // For now, just the change the temperature at the top boundary if necessary
    int dir = SpaceDim - 1;

    // Make sure we're doing dirichlet
    CH_assert(m_params.bcTypeTemperatureHi[dir] == PhysBCUtil::Dirichlet);


    // dimensionless 0.1 mangitude is about 2 degrees C
//     m_params.sinusoidal_temperature_bc_timescale = 1;
//     m_params.sinusoidal_temperature_bc_amplitude = 0.5;
//     m_params.sinusoidal_temperature_bc_av = 0.5;
//     m_params.sinusoidal_temperature_bc_phase_diff = 0.0;
//    Real T = bcValTemperatureHi[dir] - amplitude*sin(2*M_PI*a_time/timescale);

    // Ensure T isn't less than freezing
//    T = max(T, 0.001);

    m_params.bcValTemperatureHi[dir] = m_params.sinusoidal_temperature_bc_av + m_params.sinusoidal_temperature_bc_amplitude
        *sin(2*M_PI*(m_params.sinusoidal_temperature_bc_phase_diff + m_time/m_params.sinusoidal_temperature_bc_timescale) );
//    bcValEnthalpyLo[dir] = stefan + T;

    //    pout() << "Set bottom temperature BC = " << T << endl;


    // m_BCtimescale = 1;
    //   m_time = -999;
    //   m_timeDependentBC = m_constant;)
  }

}

// ---------------------------------------------------------------
Real
PhysBCUtil::Time() const
{
  return m_time;
}




