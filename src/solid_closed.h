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

#ifndef SOLID_CLOSED_H_
#define SOLID_CLOSED_H_

#include"lexer.h"
#include"fdm.h"

// D 32: the faces of the CFD projection closed inside the direct-forcing solids and the bed.
//
// The direct forcing sets the predictor to 0 in a solid, but the projection solved through it:
// the correction -dt/rho grad(p) moved the fluid across the solid again, a Darcy medium of
// permeability dt/rho.  With the density of the air (rho 1.2) a solid is no barrier for the air
// (a baffle 4 cells thick, the bays under a deck: no pressure difference held).  A closed face is
// Neumann in the pressure correction, its velocity 0 and not corrected.
//
// D 32 1: closed where the centres of both cells lie in the solid (the forcing at the surface of
//         the solid unchanged; a solid needs 2 cells across to be closed)
// D 32 2: closed where the centre of one of the cells lies in the solid (thin solids of 1 cell
//         closed as well; the wall moves to the cell face next to the solid)
//
// D 1 2 (instead of D 32 and of the forcing at these solids): the open fraction of every face from
//         the distance to the solid (solid_aperture.h): the wall where the solid cuts the faces;
//         closed where the fraction is 0
//
// Only the static solids (S, solid level set) and the bed of the direct forcing (topoforcing);
// not the 6DOF bodies (their forced velocity inside is the rigid-body motion).

// the solids are walls of the projection
inline bool solid_walls(lexer *p)
{
    return (p->D32>0 || p->D1==2) && (p->solidread>0 || p->topoforcing>0);
}

inline bool solid_closed_inside(lexer *p, fdm *a, int ii, int jj, int kk)
{
    return (p->solidread>0 && a->solid(ii,jj,kk)<0.0) || (p->topoforcing>0 && a->topo(ii,jj,kk)<0.0);
}

// the face between cell (ii,jj,kk) and cell (ii+di,jj+dj,kk+dk)
inline bool solid_face_closed(lexer *p, fdm *a, int ii, int jj, int kk, int di, int dj, int dk)
{
    if(!solid_walls(p))
    return false;

    // D 1 2: the face of the lower cell holds the open fraction (at the position of u, v, w)
    if(p->D1==2)
    {
        if(di!=0)
        return a->apx(di>0 ? ii : ii-1,jj,kk)==0.0;
        if(dj!=0)
        return a->apy(ii,dj>0 ? jj : jj-1,kk)==0.0;
        return a->apz(ii,jj,dk>0 ? kk : kk-1)==0.0;
    }

    const bool A = solid_closed_inside(p,a,ii,jj,kk);
    const bool B = solid_closed_inside(p,a,ii+di,jj+dj,kk+dk);

    return (p->D32==1) ? (A && B) : (A || B);
}

// D 1 2: a cell with all its faces closed (deep in a solid; the faces within the arrays)
inline bool solid_cell_closed(lexer *p, fdm *a, int ii, int jj, int kk)
{
    if(ii>p->imin && a->apx(ii-1,jj,kk)!=0.0) return false;
    if(ii<p->imin+p->imax-1 && a->apx(ii,jj,kk)!=0.0) return false;
    if(p->j_dir==1)
    {
    if(jj>p->jmin && a->apy(ii,jj-1,kk)!=0.0) return false;
    if(jj<p->jmin+p->jmax-1 && a->apy(ii,jj,kk)!=0.0) return false;
    }
    if(kk>p->kmin && a->apz(ii,jj,kk-1)!=0.0) return false;
    if(kk<p->kmin+p->kmax-1 && a->apz(ii,jj,kk)!=0.0) return false;
    return true;
}

// the open fraction of the face between cell (ii,jj,kk) and cell (ii+di,jj+dj,kk+dk): D 32 1/2 0 or 1
inline double solid_face_aperture(lexer *p, fdm *a, int ii, int jj, int kk, int di, int dj, int dk)
{
    if(!solid_walls(p))
    return 1.0;

    if(p->D1==2)
    {
        if(di!=0)
        return a->apx(di>0 ? ii : ii-1,jj,kk);
        if(dj!=0)
        return a->apy(ii,dj>0 ? jj : jj-1,kk);
        return a->apz(ii,jj,dk>0 ? kk : kk-1);
    }

    return solid_face_closed(p,a,ii,jj,kk,di,dj,dk) ? 0.0 : 1.0;
}

#endif
