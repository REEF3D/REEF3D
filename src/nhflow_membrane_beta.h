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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef NHFLOW_MEMBRANE_BETA_H_
#define NHFLOW_MEMBRANE_BETA_H_

// Membrane mobility (X 330) in the NHFLOW pressure Poisson equation.
//
// The membrane forcing is integrated implicitly, (1 + a K H)(u - u_m) = (u* - u_m) - a/rho grad(q),
// so in the smeared membrane layer the pressure acts with the reduced mobility
// beta = 1/(1 + a K H) (d->MBETA, cell centred). The same factor multiplies the velocity
// correction in nhflow_pjm (A 520 1, required for membranes).
//
// Row n belongs to the pressure node (i,j,k) at the bottom face of cell k, between cell k-1
// (below) and cell k (above):
//   node value         bF = min(beta(k-1), beta(k))
//   horizontal faces   harmonic mean of the node values
//   vertical faces     beta of the cell between the nodes
//   sigma cross terms  bF (both the sigxx first-derivative term and the lagged cross terms,
//                      so that a constant pressure stays in the null space)
//
// rhs0 is the right-hand side before the cross terms of this row were added.

#include"lexer.h"
#include"fdm_nhf.h"
#include"increment.h"

// mobility of the velocity correction in cell (i,j,k)
//  vertical:   cell value, the correction uses the compact difference of the pressure nodes k, k+1,
//              exactly the vertical face of the Poisson matrix
//  horizontal: minimum over the cell and its horizontal neighbours. The horizontal correction uses
//              the wide (2 dx) collocated gradient, while the matrix couples neighbours through the
//              small harmonic face value. With the cell value, a cell next to the membrane layer is
//              corrected with mobility 1 by pressure differences across the layer that the matrix
//              hardly sees, the projection over-corrects (eigenvalues of L_wide L_compact^-1 > 2)
//              and in 3D the divergence error grows from stage to stage.
#define MBETAVAL  (d->MBETA!=nullptr ? d->MBETA[IJK] : 1.0)
#define MBETAHVAL (d->MBETA!=nullptr ? nhflow_membrane_hmin(p,d,i,j,k) : 1.0)

// mobility of the faces i+1/2 and j+1/2 for the depth-jump dissipation of the HLL continuity flux.
// The dissipation c/2 (D_n - D_s) is a hydrostatic, pressure-driven mass flux. It is suppressed
//  - in the membrane layer with the mobility beta: otherwise it transports mass through the
//    membrane even when the velocity in the layer vanishes,
//  - below the bag floor, wherever the prescribed static pressure acts (chi = d->MCHI > 0):
//    there the column depth carries the inner water level, which does not drive the flow under
//    the floor. Without this the dissipation drains the bag through the gap below its edge.
#define MBETAFACEX (d->MBETA!=nullptr ? nhflow_membrane_face(d,IJK,Ip1JK) : 1.0)
#define MBETAFACEY (d->MBETA!=nullptr ? nhflow_membrane_face(d,IJK,IJp1K) : 1.0)

inline double nhflow_membrane_face(fdm_nhf *d, int q0, int q1)
{
    return (d->MCHI[q0]>0.0 || d->MCHI[q1]>0.0) ? 0.0 : MIN(d->MBETA[q0],d->MBETA[q1]);
}

inline double nhflow_membrane_hmin(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    const double *B = d->MBETA;
    double b = MIN(B[IJK],MIN(B[Ip1JK],B[Im1JK]));
    
    if(p->j_dir==1)
    b = MIN(b,MIN(B[IJp1K],B[IJm1K]));
    
    return b;
}

inline void nhflow_membrane_row(lexer *p, fdm_nhf *d, int i, int j, int k, int n, double ct, double cb, double sxx, double rhs0)
{
    // ct, cb: vertical Laplacian coefficients of the row (M.t, M.b without the sigxx term)
    // sxx:    sigxx first-derivative coefficient (M.t = ct - sxx, M.b = cb + sxx)
    const int marge = increment::marge;
    const double *B = d->MBETA;
    
    const double bc  = B[IJK];
    const double bcm = k>0 ? B[IJKm1] : bc;
    const double bF  = MIN(bc,bcm);
    
    const double bFn = k>0 ? MIN(B[Ip1JK],B[Ip1JKm1]) : B[Ip1JK];
    const double bFs = k>0 ? MIN(B[Im1JK],B[Im1JKm1]) : B[Im1JK];
    const double bFw = p->j_dir==1 ? (k>0 ? MIN(B[IJp1K],B[IJp1Km1]) : B[IJp1K]) : 1.0;
    const double bFe = p->j_dir==1 ? (k>0 ? MIN(B[IJm1K],B[IJm1Km1]) : B[IJm1K]) : 1.0;
    
    if(bF==1.0 && bcm==1.0 && bFn==1.0 && bFs==1.0 && bFw==1.0 && bFe==1.0)
    return;
    
    auto harm = [](double a, double b) {return 2.0*a*b/(a+b);};
    
    const double bn = harm(bF,bFn);
    const double bs = harm(bF,bFs);
    const double bw = harm(bF,bFw);
    const double be = harm(bF,bFe);
    const double bt = bc;
    const double bb = bcm;
    
    d->M.n[n] *= bn;
    d->M.s[n] *= bs;
    d->M.w[n] *= bw;
    d->M.e[n] *= be;
    
    d->M.t[n] = bt*ct - bF*sxx;
    d->M.b[n] = bb*cb + bF*sxx;
    
    d->M.p[n] = -(d->M.n[n] + d->M.s[n] + d->M.w[n] + d->M.e[n]) - bt*ct - bb*cb;
    
    d->rhsvec.V[n] = rhs0 + bF*(d->rhsvec.V[n] - rhs0);
}

#endif
