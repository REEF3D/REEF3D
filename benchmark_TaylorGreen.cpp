/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"benchmark_TaylorGreen.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"

benchmark_TaylorGreen::benchmark_TaylorGreen(lexer *p, fdm *a)
{
    double L = 1000;
    double U = 10;
    
    double x,y,z;
    
    ULOOP
    {
        x = p->pos1_x();
        y = p->pos1_y();
        z = p->pos1_z();
        a->u(i,j,k) = U*(sin(x/L)*cos(y/L)*cos(z/L));
    }   
    
    VLOOP
    {
        x = p->pos2_x();
        y = p->pos2_y();
        z = p->pos2_z();
        a->v(i,j,k) = U*(cos(x/L)*sin(y/L)*cos(z/L));
    }
    
    WLOOP
    a->w(i,j,k) = 0.0;
}

benchmark_TaylorGreen::~benchmark_TaylorGreen()
{
}

void benchmark_TaylorGreen::start(lexer* p, fdm *a, ghostcell *pgc, convection *pconvec )
{
}
