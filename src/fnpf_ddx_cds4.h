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

#ifndef FNPF_DDX_CDS4_H_
#define FNPF_DDX_CDS4_H_

#include "fnpf_ddx.h"
#include "lexer.h"
#include "increment.h"
#include "slice.h"

class fnpf_ddx_cds4 final : public fnpf_ddx, public increment
{
public:
    fnpf_ddx_cds4() = default;
    virtual ~fnpf_ddx_cds4() = default;

    inline double sxx(lexer *p, slice &f) override final
    {
        const double X0 = -0.5*p->XP[IP2] + 13.0*p->XP[IP1] - 13.0*p->XP[IM1] + 0.5*p->XP[IM2];
        const double X1 = (-p->XP[IP3] + 27.0*p->XP[IP2] -27.0*p->XP[IP1] + p->XP[IP])*X0;
        const double X2 = (-p->XP[IP2] + 27.0*p->XP[IP1] -27.0*p->XP[IP] + p->XP[IM1])*X0;
        const double X3 = (-p->XP[IP1] + 27.0*p->XP[IP] -27.0*p->XP[IM1] + p->XP[IM2])*X0;
        const double X4 = (-p->XP[IP] + 27.0*p->XP[IM1] -27.0*p->XP[IM2] + p->XP[IM3])*X0;

        return (      f(i+3,j) - 27.0*f(i+2,j) + 27.0*f(i+1,j) - f(i,j))/X1
             + (-27.0*f(i+2,j) + 729.0*f(i+1,j) - 729.0*f(i,j) + 27.0*f(i-1,j))/X2
             + ( 27.0*f(i+1,j) - 729.0*f(i,j) + 729.0*f(i-1,j) - 27.0*f(i-2,j))/X3
             + (-f(i,j) + 27.0*f(i-1,j) - 27.0*f(i-2,j) + f(i-3,j))/X4;
    }

    inline double syy(lexer *p, slice &f) override final
    {
        const double Y0 = -0.5*p->YP[JP2] + 13.0*p->YP[JP1] - 13.0*p->YP[JM1] + 0.5*p->YP[JM2];
        const double Y1 = (-p->YP[JP3] + 27.0*p->YP[JP2] -27.0*p->YP[JP1] + p->YP[JP])*Y0;
        const double Y2 = (-p->YP[JP2] + 27.0*p->YP[JP1] -27.0*p->YP[JP] + p->YP[JM1])*Y0;
        const double Y3 = (-p->YP[JP1] + 27.0*p->YP[JP] -27.0*p->YP[JM1] + p->YP[JM2])*Y0;
        const double Y4 = (-p->YP[JP] + 27.0*p->YP[JM1] -27.0*p->YP[JM2] + p->YP[JM3])*Y0;

        return (      f(i,j+3) - 27.0*f(i,j+2) + 27.0*f(i,j+1) - f(i,j))/Y1
             + (-27.0*f(i,j+2) + 729.0*f(i,j+1) - 729.0*f(i,j) + 27.0*f(i,j-1))/Y2
             + ( 27.0*f(i,j+1) - 729.0*f(i,j) + 729.0*f(i,j-1) - 27.0*f(i,j-2))/Y3
             + (-f(i,j) + 27.0*f(i,j-1) - 27.0*f(i,j-2) + f(i,j-3))/Y4;
    }
};

#endif
