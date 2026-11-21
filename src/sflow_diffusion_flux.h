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

#ifndef SFLOW_DIFFUSION_FLUX_H_
#define SFLOW_DIFFUSION_FLUX_H_

// The horizontal diffusion of SFLOW in flux form (sflow_ediff, sflow_idiff; as NHFLOW hdiffusion_brk):
//   d(hu)/dt = d/dx (2 (h nu)_f du/dx) + d/dy ((h nu)_f (du/dy + dv/dx)),   y-momentum likewise,
// with (h nu)_f of a face between two wet fluid cells = the mean of the two viscosities times the smaller of the two
// depths. The faces cancel in the sum over the cells, so the diffusion has no net force. The velocity form
// h_P div(nu grad u) is not conservative: at a breaking front (A 246) u and h fall together, the term
// -nu du/dx dh/dx is a net force against the waves, and with the viscosity of the cell alone a jump of nu
// between two cells (the edge of a breaking zone) gives the two cells different stresses at their common face.
// h_f <= h_P keeps the explicit limit that of nu. A face to a dry cell carries no stress (free slip at the shoreline);
// a face to a wall or an open boundary (a ghost cell) keeps the cell's own h nu with the ghost value.

#include"lexer.h"
#include"fdm2D.h"
#include"slice.h"
#include"increment.h"

namespace sflow_dflux
{
    inline int q(lexer *p, int a, int c)
    {
        return (a-p->imin)*p->jmax + (c-p->jmin);
    }

    // the viscosity of cell (a,c): molecular, eddy and breaking (A 246 2: A 250 in the breaking cells, A 246 4: vb)
    inline double nu(lexer *p, fdm2D *b, int a, int c)
    {
        double v = p->W2 + b->eddyv(a,c);

        if(p->A246==2 && b->breaking(a,c)==1)
        v += p->A250;

        if(p->A246==4)
        v += b->vb(a,c);

        return v;
    }

    inline bool fluid(lexer *p, int a, int c)
    {
        const int k = q(p,a,c);
        return p->flagslice4[k]>0 && p->wet[k]==1;
    }

    // (h nu) of the face between cell (i,j) and its neighbour (a,c): 0 towards a dry cell, the cell's own h nu
    // towards a ghost cell (wall, open boundary), else mean(nu) min(h)
    inline double hnu(lexer *p, fdm2D *b, slice &WL, int i, int j, int a, int c)
    {
        const int k = q(p,a,c);
        const double hP = WL(i,j)>0.0 ? WL(i,j) : 0.0;

        if(p->flagslice4[k]<0)
        return nu(p,b,i,j)*hP;

        if(p->wet[k]==0)
        return 0.0;

        const double hN = WL(a,c)>0.0 ? WL(a,c) : 0.0;

        return 0.5*(nu(p,b,i,j) + nu(p,b,a,c))*(hP<hN ? hP : hN);
    }

    // the cross terms in flux form, d/dy ((h nu)_f dv/dx) of the x-momentum (comp 1) and d/dx ((h nu)_f du/dy) of
    // the y-momentum (comp 2), f = V or U, with the derivative at a face the mean of the central differences of its
    // two cells; 0 unless the eight neighbours are wet fluid cells (no cross term next to dry cells and walls)
    inline double cross(lexer *p, fdm2D *b, slice &WL, slice &f, int i, int j, int comp)
    {
        if(p->j_dir==0)
        return 0.0;

        for(int a=i-1; a<=i+1; ++a)
        for(int c=j-1; c<=j+1; ++c)
        if(!fluid(p,a,c))
        return 0.0;

        const int m = increment::marge;

        if(comp==1)
        {
        const double sx = p->DXP[i+m] + p->DXP[i-1+m];
        const double gn = 0.5*((f(i+1,j) - f(i-1,j)) + (f(i+1,j+1) - f(i-1,j+1)))/sx;
        const double gs = 0.5*((f(i+1,j) - f(i-1,j)) + (f(i+1,j-1) - f(i-1,j-1)))/sx;

        return (hnu(p,b,WL,i,j,i,j+1)*gn - hnu(p,b,WL,i,j,i,j-1)*gs)/p->DYN[j+m];
        }

        const double sy = p->DYP[j+m] + p->DYP[j-1+m];
        const double ge = 0.5*((f(i,j+1) - f(i,j-1)) + (f(i+1,j+1) - f(i+1,j-1)))/sy;
        const double gw = 0.5*((f(i,j+1) - f(i,j-1)) + (f(i-1,j+1) - f(i-1,j-1)))/sy;

        return (hnu(p,b,WL,i,j,i+1,j)*ge - hnu(p,b,WL,i,j,i-1,j)*gw)/p->DXN[i+m];
    }
}

#endif
