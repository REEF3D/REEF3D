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

    // Face divided differences dq(ii,j) = (f(ii+1,j)-f(ii,j))/DXP[ii] over every face the
    // WENO5 stencils of the interior cells touch, so that dswenox_dq() needs no divisions.
    void dsdiffx(slice&, slice&);
    void dsdiffy(slice&, slice&);

    // Upwinded WENO5 gradient at (i,j) from dsdiffx/dsdiffy output, bit-identical to dswenox/dswenoy.
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


private:
    // wet-dry aware stencils, all q zero unless every cell of the stencil is wet
    inline void isqmin(slice&);
    inline void jsqmin(slice&);
    inline void isqmax(slice&);
    inline void jsqmax(slice&);

    // WENO5 reconstruction from q1..q5
    inline double weno_min_x()
    {
        is_min_x();
        weight_min_x();

        return w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
             + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
             + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
    }
    inline double weno_max_x()
    {
        is_max_x();
        weight_max_x();

        return w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
             + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
             + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
    }

    inline double weno_min_y()
    {
        is_min_y();
        weight_min_y();

        return w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
             + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
             + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
    }
    inline double weno_max_y()
    {
        is_max_y();
        weight_max_y();

        return w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
             + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
             + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
    }

    inline double weno_min_z()
    {
        is_min_z();
        weight_min_z();

        return w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
             + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
             + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
    }
    inline double weno_max_z()
    {
        is_max_z();
        weight_max_z();

        return w1z*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
             + w2z*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
             + w3z*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
    }

    lexer *p;
};

#endif
