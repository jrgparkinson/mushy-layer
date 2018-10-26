#ifndef _BCFUNCTIONS_H_
#define _BCFUNCTIONS_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"

void applyCorrectScalarBC(FArrayBox&      a_state,
                    LevelData<FluxBox>* a_advVel,
                    const Box&      a_valid,
                    bool            a_homogeneous,
                    Real            a_diriValue,
                    Real            a_inflowValue,
                    Real            a_neumValue,
                    int             a_dir,
                    Side::LoHiSide  a_side,
                    Real            a_dx,
                    int             a_order,
                    int             a_bcType,
                    Vector<Real>    a_plumeBounds,
                    Real fine_dx,
                    int a_comp = 0);

void ConstantDiriBC(FArrayBox&      a_state,
            const Box&      a_valid,
            bool            a_homogeneous,
            Real   			a_value,
            int             a_dir,
            Side::LoHiSide  a_side,
            int             a_order,
            int a_comp = 0);

void ConstantNeumBC(FArrayBox&      a_state,
            const Box&      a_valid,
            bool            a_homogeneous,
            Real   			a_value,
            int             a_dir,
            Side::LoHiSide  a_side,
            Real			a_dx,
            int a_comp = 0);

void ExtrapBC(FArrayBox&      a_state,
            const Box&      a_valid,
            int             a_dir,
            Side::LoHiSide  a_side,
            int             a_order,
            int a_comp = 0);

// Enforce dirichlet BCs with a_value if it's inflow,
// constant flux BCs if outflow
void InflowOutflowBC(FArrayBox&      a_state,
                    const LevelData<FluxBox>*      a_advVel,
            const Box&      a_valid,
            bool            a_homogeneous,
            Real                        a_value,
            int             a_dir,
            Side::LoHiSide  a_side,
            int             a_order,
            Real a_dx,
            Real fine_dx,
            int a_comp = 0);

void OnlyInflowBC(FArrayBox&      a_state,
                     const Box&      a_valid,
                     bool            a_homogeneous,
                     Real                        a_inflowValue,
                     Real                        a_otherValue,
                     Vector<Real> a_plumeBounds,
                     Real a_dx,
                     int             a_dir,
                     Side::LoHiSide  a_side,
                     int             a_order,
                     int a_comp = 0);

/// Inflow/outflow BCs for pressure
/** inflow - neumann
* outflow - dirichlet
*/
void PressureInflowOutflow(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            LevelData<FluxBox>* a_advVel,
            Real   a_diriValue,
            Real  a_neumValue,
            int             a_dir,
            Side::LoHiSide  a_side,
            int             a_order,
            Real a_fineDx,
            int a_comp = 0);

/// Plume inflow BCs
/**
 * Plume exists between a_plumeBounds[0] and a_plumeBounds[1].
 * Enforce neumann BCs for plume, dirichlet BCs for everywhere else
 */
void PressurePlumeInflow(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            Vector<Real> a_plumeBounds,
            Real   a_diriValue,
            Real  a_neumValue,
            int             a_dir,
            Side::LoHiSide  a_side,
            int             a_order,
            int a_comp = 0);


#include "NamespaceFooter.H"

#endif
