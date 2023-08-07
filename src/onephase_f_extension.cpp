/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the B117, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"onephase_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"

void onephase_f::uvel(lexer *p, fdm *a, ghostcell *pgc, field &u)
{
    UAIRLOOP
    urk1(i,j,k) = u(i,j,k) + dt*ddwenox(a->phi,a->u(i,j,k))*ddwenox(u,a->u(i,j,k));
    
    pgc->start1(p,urk1,gcval_u);

    UAIRLOOP
    u(i,j,k) = 0.5*u(i,j,k) + 0.5*urk1(i,j,k) + 0.5*dt*ddwenox(a->phi,a->u(i,j,k))*ddwenox(urk1,a->u(i,j,k));
}

void onephase_f::vvel(lexer *p, fdm *a, ghostcell*, field&)
{
}

void onephase_f::wvel(lexer *p, fdm *a, ghostcell *pgc, field &w)
{
    WAIRLOOP
    wrk1(i,j,k) = w(i,j,k) + dt*ddwenox(a->phi,a->w(i,j,k))*ddwenox(w,a->w(i,j,k));
    
    pgc->start1(p,wrk1,gcval_w);

    WAIRLOOP
    w(i,j,k) = 0.5*w(i,j,k) + 0.5*wrk1(i,j,k) + 0.5*dt*ddwenox(a->phi,a->w(i,j,k))*ddwenox(wrk1,a->w(i,j,k));
}