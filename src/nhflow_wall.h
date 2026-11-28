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

#ifndef NHFLOW_WALL_H_
#define NHFLOW_WALL_H_

// Walls of a wet NHFLOW cell, shared by the momentum wall friction (nhflow_bcmom, A519) and the
// turbulence wall functions (nhflow_kepsilon_bc, nhflow_komega_bc, nhflow_strain::wallf_update),
// so that both see the same walls.
//
//   face 0 x-, 1 x+, 2 y-, 3 y+, 4 bed, 5 top
//
// - x: a solid neighbour (DF<0) or a boundary ghost that is not an in/outflow cell (IO 0). The domain
//   ends in x are open boundaries (inflow, outflow, wave generation and absorption): never a wall.
// - y: as x, only in 3D. The domain sides in y are walls.
// - bed: k==0, or a solid / dry cell below.
// - top: a solid cell above (e.g. under a fixed structure), not the free surface.
//
// Momentum (nhflow_bcmom): A519 1 bed friction; A519 2 friction at every wall face.
// Turbulence: the bed always has a wall function (bed shear from A519 1/2 at the wall-function bed A518 1, or
// from the diffusion at the resolved no-slip bed A518 2, which has no bed friction). Side and top walls only with A519 2; otherwise they are slip walls in the momentum
// equations and the turbulence sees them with zero gradient, without a wall function.
// The faces at fixed solids (direct forcing, DF < 0) and, with A 541, at dry land above the water level of a deep
// cell (an abutment or a quay as a raised bed, A 580 0; not the shoreline of a beach) have a wall function for the turbulence whatever A519 (nhflow_turb_solid):
// the velocity there is zero, and without the wall function the jump to it in the strain made the production
// (double abutment, k-omega, A 519 1: k up to 0.18 m2/s2, nu_t 0.008 m2/s at the abutments).
// B11 0 switches the turbulence wall functions off (and with them the wall cells in WALLF).
// Sharp thin bodies (nhflow_thinbody, X 330 'mobility sharp'): the faces at the body are walls with a wall
// function (ks = B57), and with a turbulence model (A 560 > 0) wall friction in the momentum equations.

#include"lexer.h"
#include"fdm_nhf.h"
#include"increment.h"
#include"nhflow_thinbody.h"
#include<cmath>

inline void nhflow_wall_faces(lexer *p, int i, int j, int k, int *w)
{
    w[0] = (p->DF[Im1JK]<0 || (p->flag4[Im1JK]<0 && p->IO[Im1JK]==0)) && i+p->origin_i!=0;
    w[1] = (p->DF[Ip1JK]<0 || (p->flag4[Ip1JK]<0 && p->IO[Ip1JK]==0)) && i+p->origin_i!=p->gknox-1;
    w[2] = p->j_dir==1 && (p->DF[IJm1K]<0 || (p->flag4[IJm1K]<0 && p->IO[IJm1K]==0));
    w[3] = p->j_dir==1 && (p->DF[IJp1K]<0 || (p->flag4[IJp1K]<0 && p->IO[IJp1K]==0));
    w[4] = k==0 || p->flag4[IJKm1]<0 || p->DF[IJKm1]<0;
    w[5] = k!=p->knoz-1 && (p->flag4[IJKp1]<0 || p->DF[IJKp1]<0);
}

// horizontal faces at fixed solids (DF < 0) and, with A 541, at dry land above the water level of a deep cell: walls
// for the turbulence (face 0 x-, 1 x+, 2 y-, 3 y+)
inline void nhflow_turb_solid(lexer *p, fdm_nhf *d, int i, int j, int k, int *w)
{
    const int marge = increment::marge;   // for IP, JP

    // dry land is a wall for a column deep enough for A 541 (nhflow_drywall.h): not at the shoreline of a beach
    const double level = d->bed(i,j) + d->WL(i,j);
    const bool deep = p->A541>0.0 && d->WL(i,j) > p->A541*MAX(p->DXN[IP], p->j_dir==1 ? p->DYN[JP] : 0.0);

    auto land = [&](int a, int b)
    {
        const int q = (i+a-p->imin)*p->jmax + (j+b-p->jmin);
        return deep && p->flagslice4[q]>0 && p->wet[q]==0 && d->bed(i+a,j+b) > level;
    };

    w[0] = p->DF[Im1JK]<0 || land(-1,0);
    w[1] = p->DF[Ip1JK]<0 || land(1,0);
    w[2] = p->j_dir==1 && (p->DF[IJm1K]<0 || land(0,-1));
    w[3] = p->j_dir==1 && (p->DF[IJp1K]<0 || land(0,1));
}

// nearest wall with a turbulence wall function: distance of the cell centre, roughness and the velocity
// magnitude tangential to that wall. Returns 0 if the cell has none.
inline int nhflow_turb_wall(lexer *p, fdm_nhf *d, int i, int j, int k, double &dist, double &ks, double &ut)
{
    if(p->B11!=1 || p->DF[IJK]<0)
    return 0;

    const int marge = increment::marge;   // for IP, JP, KP
    int w[6];
    nhflow_wall_faces(p,i,j,k,w);

    const double U = d->U[IJK];
    const double V = d->V[IJK];
    const double W = d->W[IJK];
    const double dz = 0.5*p->DZN[KP]*d->WL(i,j);
    int found=0;

    dist = 1.0e20;

    if(w[4])
    {
    dist = dz;
    ks = d->ks(i,j);
    ut = sqrt(U*U + V*V);
    found = 1;
    }

    if(p->A519==2)
    {
        if((w[0] || w[1]) && 0.5*p->DXN[IP]<dist)
        {
        dist = 0.5*p->DXN[IP];
        ks = p->B57;
        ut = sqrt(V*V + W*W);
        found = 1;
        }

        if((w[2] || w[3]) && 0.5*p->DYN[JP]<dist)
        {
        dist = 0.5*p->DYN[JP];
        ks = p->B57;
        ut = sqrt(U*U + W*W);
        found = 1;
        }

        if(w[5] && dz<dist)
        {
        dist = dz;
        ks = p->B57;
        ut = sqrt(U*U + V*V);
        found = 1;
        }
    }

    // fixed solids and dry land above the water level: walls with a wall function whatever A519
    {
        int sw[4];
        nhflow_turb_solid(p,d,i,j,k,sw);

        if((sw[0] || sw[1]) && 0.5*p->DXN[IP]<dist)
        {
        dist = 0.5*p->DXN[IP];
        ks = p->B57;
        ut = sqrt(V*V + W*W);
        found = 1;
        }

        if((sw[2] || sw[3]) && 0.5*p->DYN[JP]<dist)
        {
        dist = 0.5*p->DYN[JP];
        ks = p->B57;
        ut = sqrt(U*U + W*W);
        found = 1;
        }
    }

    // sharp thin bodies (X 330 'mobility sharp'): the faces at the body are walls with a wall function
    if(d->thinbody!=nullptr)
    {
        int tw[6];
        d->thinbody->wall_faces(p,i,j,k,tw);
        
        if((tw[0] || tw[1]) && 0.5*p->DXN[IP]<dist)
        {
        dist = 0.5*p->DXN[IP];
        ks = p->B57;
        ut = sqrt(V*V + W*W);
        found = 1;
        }
        
        if((tw[2] || tw[3]) && 0.5*p->DYN[JP]<dist)
        {
        dist = 0.5*p->DYN[JP];
        ks = p->B57;
        ut = sqrt(U*U + W*W);
        found = 1;
        }
        
        if((tw[4] || tw[5]) && dz<dist)
        {
        dist = dz;
        ks = p->B57;
        ut = sqrt(U*U + V*V);
        found = 1;
        }
    }

    if(found && ks<=0.0)
    ks = 0.0001;   // same clamp as CFD roughness::ks_val

    return found;
}

// Implicit dissipation of the wall-function k equation, cmu^3/4 k^1/2 u+/y k = cw k^1/2 k (nhflow_kepsilon_bc,
// nhflow_komega_bc): the k^1/2 of the coefficient is not that of the last step but the root k* of the cell's own
// equation with the neighbours at their last values,
//     a k* + cw k*^3/2 = R,   a = M.p (the row without the wall term), R = rhs - sum_nb M.nb k_nb,
// with the row as the solver gets it (after bckin_matrix). With k^1/2 of the last step, a cell whose k was ~0 got a
// dissipation of ~0 against the full production tau u/y in one time step: k = (k_eq^3/2 / k_old^1/2), k_eq the
// equilibrium value u*^2/cmu^1/2. The long steps of the semi-implicit momentum (dt 20-60 s on the tidal flat of
// nhflow_2d_si_tidal_flat_kepsilon_elm) gave k = 10-150 m2/s2 in the bed cell of the column at the wet-dry front for
// one step (k_old 1e-5, k_eq 0.6), and through the vertical diffusion in the whole column. With k*: k <= k_eq from
// the wall terms alone. In a steady state k* is the k of the last step, the solution is that of the lagged form.
// Returns k* (R <= 0: k of the last step, the row is a sink anyway).
inline double nhflow_wall_kstar(lexer *p, fdm_nhf *d, const double *KIN, int n, int i, int j, int k, double cw)
{
    const int marge = increment::marge;   // for the IJK macros

    const double R = d->rhsvec.V[n] - d->M.n[n]*KIN[Ip1JK] - d->M.s[n]*KIN[Im1JK]
                                     - d->M.w[n]*KIN[IJp1K] - d->M.e[n]*KIN[IJm1K]
                                     - d->M.t[n]*KIN[IJKp1] - d->M.b[n]*KIN[IJKm1];
    const double a = d->M.p[n]>0.0 ? d->M.p[n] : 0.0;

    if(R<=0.0 || cw<=0.0)
    return KIN[IJK]>0.0 ? KIN[IJK] : 0.0;

    // f(s) = cw s^3 + a s^2 - R, s = k^1/2: increasing and convex for s >= 0; Newton from an upper bound of the root
    // (min of the roots of the two terms alone) converges monotonically from above
    double s = std::cbrt(R/cw);

    if(a>0.0)
    s = std::min(s, std::sqrt(R/a));

    for(int it=0; it<30; ++it)
    {
        const double fs = (cw*s + a)*s*s - R;
        const double df = (3.0*cw*s + 2.0*a)*s;

        if(df<=0.0)
        break;

        const double ds = fs/df;
        s -= ds;

        if(fabs(ds)<=1.0e-10*s)
        break;
    }

    return s>0.0 ? s*s : 0.0;
}

#endif
