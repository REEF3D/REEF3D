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

#ifndef FNPF_HIRES_H_
#define FNPF_HIRES_H_

#include "fnpf_convection.h"
#include "increment.h"
#include "lexer.h"
#include "slice.h"

#include <cmath>

class fnpf_hires final : public fnpf_convection, public increment
{
public:
    fnpf_hires() = default;
    virtual ~fnpf_hires() = default;

    inline double fx(lexer*, field&, double, double) override final {return 0.0;};
    inline double fy(lexer*, field&, double, double) override final {return 0.0;};
    inline double fz(lexer*, field&, double, double) override final {return 0.0;};

    inline double sx(lexer *p, slice &f, double ivel) override final
    {
        const double dfdx_plus = (f(i+1,j) - f(i,j))/p->DXP[IP];
        const double dfdx_min  = (f(i,j) - f(i-1,j))/p->DXP[IM1];

        return limiter(dfdx_plus,dfdx_min);
    }

    inline double sy(lexer *p, slice &f, double ivel) override final
    {
        const double dfdy_plus = (f(i,j+1) - f(i,j))/p->DYP[JP];
        const double dfdy_min  = (f(i,j) - f(i,j-1))/p->DYP[JM1];

        return limiter(dfdy_plus,dfdy_min);
    }

    inline double sz(lexer*, double*) override final {return 0.0;};

private:

    inline double limiter(double v1, double v2)
    {
        double denom = fabs(v1) + fabs(v2);

        denom = fabs(denom)>1.0e-10?denom:1.0e10;

        return (v1*fabs(v2) + fabs(v1)*v2)/denom;
    }
};

#endif
