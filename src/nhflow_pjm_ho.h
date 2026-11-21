/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Architect: Hans Bihs
--------------------------------------------------------------------*/

#ifndef NHFLOW_PJM_HO_H_
#define NHFLOW_PJM_HO_H_

// Fourth-order horizontal operators of the explicit NHFLOW projection (A 510 2 / 3, A 520 1 / 2),
// shared by nhflow_pjm and nhflow_pjm_corr.
//
// The divergence at the pressure nodes and the gradient of the non-hydrostatic pressure at the cell centres are
// central differences over 2 dx on the collocated grid: (f_{i+1} - f_{i-1})/(2 dx). Their product has the symbol
// sin^2(k dx)/dx^2 instead of k^2, so the non-hydrostatic pressure of a wave is that of a longer wave and the wave
// runs too fast, while the hydrostatic part (WENO5 + HLL) has no such error: standing waves over a flat bed with
// 8 layers, kh 3.6 at 30 cells per wave length: period -0.27 %, kh 6.3 at 18 cells: -0.81 % (linear theory); on the
// Beji & Battjes bar the free third and fourth harmonics behind the bar ran 0.8 % and 2.4 % too fast.
// With (8 (f_{i+1} - f_{i-1}) - (f_{i+2} - f_{i-2}))/(12 dx) for both: +0.13 % and +0.26 %.
//
// The matrix of the Poisson equation stays the compact second-order one. It is the preconditioner of the
// projection: with the pressure increment (A 520 2) the divergence left by a step is taken out by the next ones.
//
// Fourth order where the two neighbours on both sides are regular water cells on a uniform Cartesian grid;
// second order elsewhere (next to walls and open boundaries, dry or shallow columns, solids) and on curvilinear,
// multi-block and refined grids, with membranes, thin bodies and floating bodies.

#include"lexer.h"
#include"fdm_nhf.h"
#include"increment.h"
#include<cmath>

namespace nhflow_pjm_ho
{

inline bool active(lexer *p, fdm_nhf *d)
{
    return p->geo_curv==0 && p->G1==0 && p->geo_mb==0 && p->X10==0 && d->MBETA==nullptr && d->thinbody==nullptr;
}

inline bool column(lexer *p, int a, int b)
{
    const int q = (a-p->imin)*p->jmax + (b-p->jmin);
    return p->wet[q]==1 && p->deep[q]==1;
}

inline bool cell(lexer *p, fdm_nhf *d, int a, int b, int c)
{
    const int q = (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + c-p->kmin;
    return p->flag4[q]>0 && !(d->SOLID[q]<0.0);   // SOLID: 0 everywhere without solids
}

// the stencil i-2 .. i+2 of the cell (i,j,k)
inline bool x(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    const int marge = increment::marge;
    const double dx = p->DXP[IP];
    
    if(fabs(p->DXP[IM1]-dx)>1.0e-8*dx || fabs(p->DXP[IM2]-dx)>1.0e-8*dx || fabs(p->DXP[IP1]-dx)>1.0e-8*dx)
    return false;
    
    for(int q=-2; q<=2; ++q)
    if(q!=0 && (!cell(p,d,i+q,j,k) || !column(p,i+q,j)))
    return false;
    
    return true;
}

// the stencil j-2 .. j+2 of the cell (i,j,k)
inline bool y(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    const int marge = increment::marge;
    const double dy = p->DYP[JP];
    
    if(p->j_dir==0)
    return false;
    
    if(fabs(p->DYP[JM1]-dy)>1.0e-8*dy || fabs(p->DYP[JM2]-dy)>1.0e-8*dy || fabs(p->DYP[JP1]-dy)>1.0e-8*dy)
    return false;
    
    for(int q=-2; q<=2; ++q)
    if(q!=0 && (!cell(p,d,i,j+q,k) || !column(p,i,j+q)))
    return false;
    
    return true;
}

// d/dx and d/dy of the layer mean (P_k + P_{k+1})/2 of a nodal field at the centre of the cell (i,j,k)
inline double gradx(lexer *p, double *P, int i, int j, int k)
{
    const int marge = increment::marge;
    auto m = [&](int q) { const int n = (i+q-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin; return 0.5*(P[n]+P[n+1]); };
    
    return (8.0*(m(1)-m(-1)) - (m(2)-m(-2)))/(6.0*(p->DXP[IP]+p->DXP[IM1]));
}

inline double grady(lexer *p, double *P, int i, int j, int k)
{
    const int marge = increment::marge;
    auto m = [&](int q) { const int n = (i-p->imin)*p->jmax*p->kmaxF + (j+q-p->jmin)*p->kmaxF + k-p->kmin; return 0.5*(P[n]+P[n+1]); };
    
    return (8.0*(m(1)-m(-1)) - (m(2)-m(-2)))/(6.0*(p->DYP[JP]+p->DYP[JM1]));
}

// d/dx and d/dy at the pressure node k of the column (i,j) of a cell field: the node value is (1-fac) f_k + fac f_{k-1},
// below the bed node 0 (nhflow_pjm::rhs)
inline double divx(lexer *p, double *U, int i, int j, int k, double fac)
{
    const int marge = increment::marge;
    auto m = [&](int q) { const int n = (i+q-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin; return (1.0-fac)*U[n] + (k==0 ? 0.0 : fac*U[n-1]); };
    
    return (8.0*(m(1)-m(-1)) - (m(2)-m(-2)))/(6.0*(p->DXP[IP]+p->DXP[IM1]));
}

inline double divy(lexer *p, double *V, int i, int j, int k, double fac)
{
    const int marge = increment::marge;
    auto m = [&](int q) { const int n = (i-p->imin)*p->jmax*p->kmax + (j+q-p->jmin)*p->kmax + k-p->kmin; return (1.0-fac)*V[n] + (k==0 ? 0.0 : fac*V[n-1]); };
    
    return (8.0*(m(1)-m(-1)) - (m(2)-m(-2)))/(6.0*(p->DYP[JP]+p->DYP[JM1]));
}

}

#endif
