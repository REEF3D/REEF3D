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

#ifndef FNPF_DDWENO_F_NUG_H_
#define FNPF_DDWENO_F_NUG_H_

#include "increment.h"
#include "weno_nug_func.h"

class field;
class slice;
class lexer;

using namespace std;

class fnpf_ddweno_f_nug : public weno_nug_func
{
public:
    fnpf_ddweno_f_nug(lexer*);
    ~fnpf_ddweno_f_nug();

    // field
    double ddwenox(field&, double);
    double ddwenoy(field&, double);
    double ddwenoz(field&, double);

    // slice, wet-dry aware: zero gradient unless every cell of the stencil is wet
    double dswenox(slice&, double);
    double dswenoy(slice&, double);

    // Upwinded WENO5 gradient at (i,j) from weno_nug_func::dsdiffx/dsdiffy output, bit-identical to dswenox/dswenoy.
    inline double dswenox_dq(slice &dq, double uvel)
    {
        if(uvel>0.0)
        {
            if(p->wet[Im3J]>0 && p->wet[Im2J]>0 && p->wet[Im1J]>0 && p->wet[IJ]>0 && p->wet[Ip1J]>0 && p->wet[Ip2J]>0)
            {
                q1 = dq(i-3,j);
                q2 = dq(i-2,j);
                q3 = dq(i-1,j);
                q4 = dq(i,j);
                q5 = dq(i+1,j);
                return weno_min_x();
            }
            return 0.0;
        }
        else if(uvel<0.0)
        {
            if(p->wet[Im2J]>0 && p->wet[Im1J]>0 && p->wet[IJ]>0 && p->wet[Ip1J]>0 && p->wet[Ip2J]>0 && p->wet[Ip3J]>0)
            {
                q1 = dq(i-2,j);
                q2 = dq(i-1,j);
                q3 = dq(i,j);
                q4 = dq(i+1,j);
                q5 = dq(i+2,j);
                return weno_max_x();
            }
            return 0.0;
        }
        else
            return 0.0;
    }

    inline double dswenoy_dq(slice &dq, double vvel)
    {
        if(vvel>0.0)
        {
            if(p->wet[IJm3]>0 && p->wet[IJm2]>0 && p->wet[IJm1]>0 && p->wet[IJ]>0 && p->wet[IJp1]>0 && p->wet[IJp2]>0)
            {
                q1 = dq(i,j-3);
                q2 = dq(i,j-2);
                q3 = dq(i,j-1);
                q4 = dq(i,j);
                q5 = dq(i,j+1);
                return weno_min_y();
            }
            return 0.0;
        }
        else if(vvel<0.0)
        {
            if(p->wet[IJm2]>0 && p->wet[IJm1]>0 && p->wet[IJ]>0 && p->wet[IJp1]>0 && p->wet[IJp2]>0 && p->wet[IJp3]>0)
            {
                q1 = dq(i,j-2);
                q2 = dq(i,j-1);
                q3 = dq(i,j);
                q4 = dq(i,j+1);
                q5 = dq(i,j+2);
                return weno_max_y();
            }
            return 0.0;
        }
        else
            return 0.0;
    }

    // Wet-dry aware upwinded gradient with two options (see A315, A316):
    //   zerosym  = 1: symmetric gradient at zero speed instead of 0
    //   fallback = 1: first-order gradient from wet neighbours when the WENO5
    //                 stencil is not fully wet, instead of 0
    // zerosym = fallback = 0 is identical to dswenox_dq/dswenoy_dq.
    inline double dswenox_dq_wd(slice &dq, double uvel, int zerosym, int fallback)
    {
        const bool wl = p->wet[Im3J]>0 && p->wet[Im2J]>0 && p->wet[Im1J]>0 && p->wet[IJ]>0 && p->wet[Ip1J]>0 && p->wet[Ip2J]>0;
        const bool wr = p->wet[Im2J]>0 && p->wet[Im1J]>0 && p->wet[IJ]>0 && p->wet[Ip1J]>0 && p->wet[Ip2J]>0 && p->wet[Ip3J]>0;

        if(uvel>0.0)
        {
            if(wl)
                return wmin_x(dq);
        }
        else if(uvel<0.0)
        {
            if(wr)
                return wmax_x(dq);
        }
        else
        {
            if(zerosym==0)
                return 0.0;

            if(wl && wr)
                return 0.5*(wmin_x(dq) + wmax_x(dq));
            if(wl)
                return wmin_x(dq);
            if(wr)
                return wmax_x(dq);
        }

        if(fallback==1)
            return fallback_x(dq,uvel);

        return 0.0;
    }

    inline double dswenoy_dq_wd(slice &dq, double vvel, int zerosym, int fallback)
    {
        const bool wl = p->wet[IJm3]>0 && p->wet[IJm2]>0 && p->wet[IJm1]>0 && p->wet[IJ]>0 && p->wet[IJp1]>0 && p->wet[IJp2]>0;
        const bool wr = p->wet[IJm2]>0 && p->wet[IJm1]>0 && p->wet[IJ]>0 && p->wet[IJp1]>0 && p->wet[IJp2]>0 && p->wet[IJp3]>0;

        if(vvel>0.0)
        {
            if(wl)
                return wmin_y(dq);
        }
        else if(vvel<0.0)
        {
            if(wr)
                return wmax_y(dq);
        }
        else
        {
            if(zerosym==0)
                return 0.0;

            if(wl && wr)
                return 0.5*(wmin_y(dq) + wmax_y(dq));
            if(wl)
                return wmin_y(dq);
            if(wr)
                return wmax_y(dq);
        }

        if(fallback==1)
            return fallback_y(dq,vvel);

        return 0.0;
    }


private:

    // WENO5 reconstructions from the face differences, no wet-dry checks
    inline double wmin_x(slice &dq)
    {
        q1 = dq(i-3,j);
        q2 = dq(i-2,j);
        q3 = dq(i-1,j);
        q4 = dq(i,j);
        q5 = dq(i+1,j);
        return weno_min_x();
    }
    inline double wmax_x(slice &dq)
    {
        q1 = dq(i-2,j);
        q2 = dq(i-1,j);
        q3 = dq(i,j);
        q4 = dq(i+1,j);
        q5 = dq(i+2,j);
        return weno_max_x();
    }
    inline double wmin_y(slice &dq)
    {
        q1 = dq(i,j-3);
        q2 = dq(i,j-2);
        q3 = dq(i,j-1);
        q4 = dq(i,j);
        q5 = dq(i,j+1);
        return weno_min_y();
    }
    inline double wmax_y(slice &dq)
    {
        q1 = dq(i,j-2);
        q2 = dq(i,j-1);
        q3 = dq(i,j);
        q4 = dq(i,j+1);
        q5 = dq(i,j+2);
        return weno_max_y();
    }

    // First-order one-sided gradient from the wet faces next to (i,j):
    // upwind face if wet, else the other wet face, else 0.
    // dq(i-1) is the face i-1/2, dq(i) the face i+1/2.
    inline double fallback_x(slice &dq, double uvel)
    {
        const bool bw = p->wet[Im1J]>0 && p->wet[IJ]>0;
        const bool fw = p->wet[IJ]>0 && p->wet[Ip1J]>0;

        if(bw && fw)
        {
            if(uvel>0.0)
                return dq(i-1,j);
            if(uvel<0.0)
                return dq(i,j);
            return 0.5*(dq(i-1,j) + dq(i,j));
        }
        if(bw)
            return dq(i-1,j);
        if(fw)
            return dq(i,j);

        return 0.0;
    }
    inline double fallback_y(slice &dq, double vvel)
    {
        const bool bw = p->wet[IJm1]>0 && p->wet[IJ]>0;
        const bool fw = p->wet[IJ]>0 && p->wet[IJp1]>0;

        if(bw && fw)
        {
            if(vvel>0.0)
                return dq(i,j-1);
            if(vvel<0.0)
                return dq(i,j);
            return 0.5*(dq(i,j-1) + dq(i,j));
        }
        if(bw)
            return dq(i,j-1);
        if(fw)
            return dq(i,j);

        return 0.0;
    }

    // wet-dry aware stencils, all q zero unless every cell of the stencil is wet
    inline void isqmin(slice&);
    inline void jsqmin(slice&);
    inline void isqmax(slice&);
    inline void jsqmax(slice&);

    lexer *p;
};

#endif
