/*
 * ChomboSpline.h is a port of spline.h to match the Chombo interface,
 * primarily to include Real rather than double so we can also use floats
 * spline.h
 *
 * simple cubic spline interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014 Tino Kluge (ttk448 at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 *
 */


#ifndef _CHOMBOSPLINE_H_
#define _CHOMBOSPLINE_H_

#include "IntVect.H"
#include "RealVect.H"
//#include "MushyLayerParams.h"
//#include "LevelData.H"
//#include "FArrayBox.H"
//#include "BoxIterator.H"
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>

//#include "NamespaceHeader.H"



//void analyticSolnSolidificationNoFlow(LevelData<FArrayBox>& enthalpyAnalytic,
//                LevelData<FArrayBox>& bulkConcentrationAnalytic,
//                Real a_domainLength, Real a_dx,
//                MushyLayerParams a_params);

class band_matrix
{
private:
    std::vector< std::vector<Real> > m_upper;  // upper band
    std::vector< std::vector<Real> > m_lower;  // lower band
public:
    band_matrix() {};                             // constructor
    band_matrix(int dim, int n_u, int n_l);       // constructor
    ~band_matrix() {};                            // destructor
    void resize(int dim, int n_u, int n_l);      // init with dim,n_u,n_l
    int dim() const;                             // matrix dimension
    int num_upper() const
    {
        return m_upper.size()-1;
    }
    int num_lower() const
    {
        return m_lower.size()-1;
    }
    // access operator
    Real & operator () (int i, int j);            // write
    Real   operator () (int i, int j) const;      // read
    // we can store an additional diagonal (in m_lower)
    Real& saved_diag(int i);
    Real  saved_diag(int i) const;
    void lu_decompose();
    std::vector<Real> r_solve(const std::vector<Real>& b) const;
    std::vector<Real> l_solve(const std::vector<Real>& b) const;
    std::vector<Real> lu_solve(const std::vector<Real>& b,
                                 bool is_lu_decomposed=false);

};

// spline interpolation
class spline
{
public:
    enum bd_type {
        first_deriv = 1,
        second_deriv = 2
    };

private:
    std::vector<Real> m_x,m_y;            // x,y coordinates of points
    // interpolation parameters
    // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
    std::vector<Real> m_a,m_b,m_c;        // spline coefficients

    bd_type m_left, m_right;
    Real  m_left_value, m_right_value;
    bool    m_force_linear_extrapolation;
    Real  m_b0, m_c0;                     // for left extrapol

public:
    // set default boundary condition to be zero curvature at both ends
    spline(): m_left(second_deriv), m_right(second_deriv),
        m_left_value(0.0), m_right_value(0.0),
        m_force_linear_extrapolation(false), m_b0(0.0), m_c0(0.0)
    {
        ;
    }

    // optional, but if called it has to come be before set_points()
    void set_boundary(bd_type left, Real left_value,
                      bd_type right, Real right_value,
                      bool force_linear_extrapolation=false);
    void set_points(const std::vector<Real>& x,
                    const std::vector<Real>& y, bool cubic_spline=true);
    Real operator() (Real x) const;
};



#endif
