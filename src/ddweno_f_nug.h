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

#ifndef DDWENO_F_NUG_H_
#define DDWENO_F_NUG_H_

#include "weno_nug_func.h"

class fdm;
class field;
class lexer;
class slice;

class ddweno_f_nug : public weno_nug_func
{
public:
    ddweno_f_nug(lexer*);
    ~ddweno_f_nug();

    double ddwenox(field&, double);
    double ddwenoy(field&, double);
    double ddwenoz(field&, double);

    double dswenox(slice&, double);
    double dswenoy(slice&, double);

    // Upwinded WENO5 gradient at (i,j) from weno_nug_func::dsdiffx/dsdiffy output, bit-identical to dswenox/dswenoy.
    // Zero velocity returns 0, like fnpf_weno5::sx/sy.
    inline double dswenox_dq(slice &dq, double uvel)
    {
        if(uvel>0.0)
        {
            q1 = dq(i-3,j);
            q2 = dq(i-2,j);
            q3 = dq(i-1,j);
            q4 = dq(i,j);
            q5 = dq(i+1,j);
            return weno_min_x();
        }
        else if(uvel<0.0)
        {
            q1 = dq(i-2,j);
            q2 = dq(i-1,j);
            q3 = dq(i,j);
            q4 = dq(i+1,j);
            q5 = dq(i+2,j);
            return weno_max_x();
        }
        else
            return 0.0;
    }

    inline double dswenoy_dq(slice &dq, double vvel)
    {
        if(vvel>0.0)
        {
            q1 = dq(i,j-3);
            q2 = dq(i,j-2);
            q3 = dq(i,j-1);
            q4 = dq(i,j);
            q5 = dq(i,j+1);
            return weno_min_y();
        }
        else if(vvel<0.0)
        {
            q1 = dq(i,j-2);
            q2 = dq(i,j-1);
            q3 = dq(i,j);
            q4 = dq(i,j+1);
            q5 = dq(i,j+2);
            return weno_max_y();
        }
        else
            return 0.0;
    }

protected:
    double grad;
};

#endif
