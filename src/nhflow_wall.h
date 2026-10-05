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
// Turbulence: the bed always has a wall function (bed shear from A519 1/2 or from the no-slip bed in the
// diffusion). Side and top walls only with A519 2; otherwise they are slip walls in the momentum
// equations and the turbulence sees them with zero gradient, without a wall function.
// B11 0 switches the turbulence wall functions off (and with them the wall cells in WALLF).

#include"lexer.h"
#include"fdm_nhf.h"
#include"increment.h"
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

    if(found && ks<=0.0)
    ks = 0.0001;   // same clamp as CFD roughness::ks_val

    return found;
}

#endif
