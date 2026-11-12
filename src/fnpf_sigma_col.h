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

#ifndef FNPF_SIGMA_COL_H_
#define FNPF_SIGMA_COL_H_

#include"lexer.h"
#include"fdm_fnpf.h"
#include"fnpf_memory.h"
#include<cmath>

// The sigma metrics of FNPF are affine in sigma within a column: sigx, sigy and sigxx follow from
// the depth W and the bed and surface slopes and curvatures of the column (fsfdisc), and
// sig = ZN. Instead of three 3D arrays (and sig) they are evaluated where they are read: the Laplace
// assembly, the velocities and the bed boundary condition. sigma_update fills the 3D arrays only
// for the other readers (fnpf_sigma_3D), with the same functions, so the values are the same bit
// for bit, and they are those of the old arrays: the column quantities do not change between
// sigma_update and these readers (fsfdisc and wetdry, which write them, come before sigma_update
// in every RK stage).
struct fnpf_sigma_col
{
    double W,W2,bxW,exW,byW,eyW,tbx,tex,tby,tey,sbex,bex,sbey,bey;
    
    // the column (i,j); W as WLVL in fnpf_sigma_update.cpp
    inline void set(lexer *p, fdm_fnpf *c, int i, int j)
    {
        W   = (fabs(c->WL(i,j))>1.0e-20?c->WL(i,j):1.0e-5);
        W2  = W*W;
        const double bx  = c->Bx(i,j),  by  = c->By(i,j);
        const double ex  = c->Ex(i,j),  ey  = c->Ey(i,j);
        const double bxx = c->Bxx(i,j), byy = c->Byy(i,j);
        const double exx = c->Exx(i,j), eyy = c->Eyy(i,j);
        bxW = bx/W; exW = ex/W;
        byW = by/W; eyW = ey/W;
        tbx = bxx - bx*bx/W;
        tex = exx - ex*ex/W;
        tby = byy - by*by/W;
        tey = eyy - ey*ey/W;
        sbex = bx + ex; bex = bx*ex;
        sbey = by + ey; bey = by*ey;
    }
    
    inline double sx(double s) const
    {
        const double sgx = (1.0 - s)*bxW - s*exW;
        return sgx;
    }
    
    inline double sy(double s) const
    {
        const double sgy = (1.0 - s)*byW - s*eyW;
        return sgy;
    }
    
    inline double sxx(double s, double sgx, double sgy) const
    {
        const double v = ((1.0 - s)/W)*tbx
                       - (s/W)*tex
                       - (sgx/W)*sbex
                       - ((1.0 - 2.0*s)/W2)*bex
                       + ((1.0 - s)/W)*tby
                       - (s/W)*tey
                       - (sgy/W)*sbey
                       - ((1.0 - 2.0*s)/W2)*bey;
        return v;
    }
};

#endif
