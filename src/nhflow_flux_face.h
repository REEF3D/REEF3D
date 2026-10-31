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

#ifndef NHFLOW_FLUX_FACE_H_
#define NHFLOW_FLUX_FACE_H_

// Physical face fluxes of the NHFLOW momentum and continuity equations,
// evaluated at one face (i,j,k). Single source for nhflow_flux_build_f
// (which stores them in d->Fs, d->Fn, d->Fe, d->Fw) and for the fused
// nhflow_HLL sweeps (which evaluate them on the fly).
//
// s/n: left/right state at the x-face i+1/2, e/w: at the y-face j+1/2.
//
// The face functions are noinline on purpose: g++ contracts a*b+c into FMA
// across inlined expressions (-ffp-contract=fast is its C++ default), so
// inlining them into the HLL formula would change the last bits of the
// fluxes. As separate calls they round exactly as the stored d->Fs/Fn/Fe/Fw
// did before, and the fused HLL sweep stays bitwise identical. (Measured:
// no speed difference between inline and noinline.)

#include"lexer.h"
#include"fdm_nhf.h"
#include"slice.h"
#include"vrans_definitions.h"
#include<cmath>

namespace nhflow_face
{
// porosity at the free surface (VRANS, d->PORS), otherwise of the surface cell, averaged to the x- / y-face: storage coefficient of the
// continuity equation n_s deta/dt + div(Q) = 0 -> the HLL/HLLC jump term is n_s*(D_R - D_L)
inline double por_surf_x(lexer *p, fdm_nhf *d, int i, int j)
{
    const int k = p->knoz-1;
    
    if(p->B200==1)
    return 0.5*(d->PORS(i,j) + d->PORS(i+1,j));
    
    return 0.5*(d->POR[IJK] + d->POR[Ip1JK]);
}

inline double por_surf_y(lexer *p, fdm_nhf *d, int i, int j)
{
    const int k = p->knoz-1;
    
    if(p->B200==1)
    return 0.5*(d->PORS(i,j) + d->PORS(i,j+1));
    
    return 0.5*(d->POR[IJK] + d->POR[IJp1K]);
}

// VRANS (B 200 1, B 206 1): HLL jump term of the momentum equations.
// Every momentum tendency is scaled by the inertia factor C = CPORNH (vrans_definitions.h), i.e. the
// system is M q_t + f_x = 0 with M = 1/C for the momentum q = UH. The HLL flux of that system has the
// jump M (q_R - q_L) = (q_R - q_L)/C (as the continuity has n_s (D_R - D_L)); with the jump q_R - q_L
// the numerical dissipation of the porous momentum was C times (0.17-0.27 in an armour layer) that
// of open water. Face value of C as in the advective flux (PORVALNH1, PORVALNH2); 1 in open water.
inline double jinv_x(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    if(p->B200!=1 || p->B206==0)
    return 1.0;
    
    return 1.0/CPORNHval(PORVALNH1);
}

inline double jinv_y(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    if(p->B200!=1 || p->B206==0)
    return 1.0;
    
    return 1.0/CPORNHval(PORVALNH2);
}

// VRANS (B 200 1, B 206 1): wave speeds of the porous system with M = diag(n_s, 1/C)
//      n_s deta/dt + dQ/dx = 0,      dQ/dt + C d/dx (Q^2/(n^2 D) + g D^2/2) = ...
// eigenvalues of M^-1 df/dq, u = Q/D (Darcy velocity), n the face porosity of the advective flux:
//      lambda = a +- b,   a = C u/n^2,   b^2 = C^2 u^2/n^4 + (C/n_s)(g D - u^2/n^2)
// n = n_s = C = 1 gives u +- sqrt(g D). In an armour layer (n = 0.3, B 260 0.34) the advective speed
// is about 3 u (B 265 0) or 1.9 u (B 265 1) instead of u.
inline void por_speed(lexer *p, double n, double ns, double u, double D, double &a, double &b)
{
    const double C = CPORNHval(n);
    a = C*u/(n*n);
    const double b2 = a*a + (C/ns)*(9.81*MAX(D,0.0) - u*u/(n*n));
    b = sqrt(MAX(b2,0.0));
}

// porosity (VRANS, Darcy velocities): the advective flux UH*U/n^2 of the face i+1/2 (j+1/2) is evaluated
// with the porosity of that face, PORVALNH1 (PORVALNH2), for the left and the right state alike: equal
// states give equal fluxes, and the porosity gradient acts at the face where it is

// U momentum
__attribute__((noinline)) inline double U_s(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return d->UHs[IJK]*d->Us[IJK]/(PORVALNH1*PORVALNH1)
            + 0.5*fabs(p->W22)*d->ETAs(i,j)*d->ETAs(i,j)
            + fabs(p->W22)*d->ETAs(i,j)*d->dfx(i,j);
}

__attribute__((noinline)) inline double U_n(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return d->UHn[IJK]*d->Un[IJK]/(PORVALNH1*PORVALNH1)
            + 0.5*fabs(p->W22)*d->ETAn(i,j)*d->ETAn(i,j)
            + fabs(p->W22)*d->ETAn(i,j)*d->dfx(i,j);
}

__attribute__((noinline)) inline double U_e(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return d->Ve[IJK]*d->UHe[IJK]/(PORVALNH2*PORVALNH2);
}

__attribute__((noinline)) inline double U_w(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return d->Vw[IJK]*d->UHw[IJK]/(PORVALNH2*PORVALNH2);
}

// V momentum
__attribute__((noinline)) inline double V_s(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return d->Us[IJK]*d->VHs[IJK]/(PORVALNH1*PORVALNH1);
}

__attribute__((noinline)) inline double V_n(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return d->Un[IJK]*d->VHn[IJK]/(PORVALNH1*PORVALNH1);
}

__attribute__((noinline)) inline double V_e(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return d->VHe[IJK]*d->Ve[IJK]/(PORVALNH2*PORVALNH2)
            + 0.5*fabs(p->W22)*d->ETAe(i,j)*d->ETAe(i,j)
            + fabs(p->W22)*d->ETAe(i,j)*d->dfy(i,j);
}

__attribute__((noinline)) inline double V_w(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return d->VHw[IJK]*d->Vw[IJK]/(PORVALNH2*PORVALNH2)
            + 0.5*fabs(p->W22)*d->ETAw(i,j)*d->ETAw(i,j)
            + fabs(p->W22)*d->ETAw(i,j)*d->dfy(i,j);
}

// W momentum
__attribute__((noinline)) inline double W_s(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return d->Us[IJK]*d->WHs[IJK]/(PORVALNH1*PORVALNH1);
}

__attribute__((noinline)) inline double W_n(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return d->Un[IJK]*d->WHn[IJK]/(PORVALNH1*PORVALNH1);
}

__attribute__((noinline)) inline double W_e(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return d->Ve[IJK]*d->WHe[IJK]/(PORVALNH2*PORVALNH2);
}

__attribute__((noinline)) inline double W_w(lexer *p, fdm_nhf *d, int i, int j, int k)
{
    return d->Vw[IJK]*d->WHw[IJK]/(PORVALNH2*PORVALNH2);
}

// vertical (sigma) flux, upwinded with omega: Fz = omega*Fb (omega>=0) or omega*Ft (omega<0)
inline void zflux(lexer *p, fdm_nhf *d, const double *Fb, const double *Ft)
{
    int i,j,k;
    
    WLOOP
    {
    if(d->omegaF[FIJKp1]>=0.0)
    d->Fz[IJK] = (d->omegaF[FIJKp1]*(Fb[IJK]))/(PORVALNH3*PORVALNH3);
    
    if(d->omegaF[FIJKp1]<0.0)
    d->Fz[IJK] = (d->omegaF[FIJKp1]*(Ft[IJK]))/(PORVALNH3*PORVALNH3);
    }
}
}

#endif
