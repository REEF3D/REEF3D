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
    AIRLOOP
    {
    nx = (a->phi(i+1,j,k)-a->phi(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]);
    ny = (a->phi(i,j+1,k)-a->phi(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP]);
    nz = (a->phi(i,j,k+1)-a->phi(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]);  

        dnorm = sqrt(nx*nx + ny*ny + nz*nz);
        
        nx/=(dnorm>1.0e-8?dnorm:1.0e10);
        ny/=(dnorm>1.0e-8?dnorm:1.0e10);
        nz/=(dnorm>1.0e-8?dnorm:1.0e10);
        
    xphi(i,j,k) = -nx;
    yphi(i,j,k) = -ny;
    zphi(i,j,k) = -nz;
    }
    
    for(int qn=0;qn<3;++qn)
    {
    UAIRLOOP
    urk1(i,j,k) = u(i,j,k) - dt*(xphi(i,j,k)*ddwenox(u,a->u(i,j,k)) + yphi(i,j,k)*ddwenoy(u,a->u(i,j,k)) + zphi(i,j,k)*ddwenoz(u,a->u(i,j,k)) );
    
    pgc->start1(p,urk1,gcval_u);

    UAIRLOOP
    u(i,j,k) = 0.5*u(i,j,k) + 0.5*urk1(i,j,k) - 0.5*dt*(xphi(i,j,k)*ddwenox(urk1,a->u(i,j,k)) + yphi(i,j,k)*ddwenoy(urk1,a->u(i,j,k)) + zphi(i,j,k)*ddwenoz(urk1,a->u(i,j,k)) );
    }
}

void onephase_f::vvel(lexer *p, fdm *a, ghostcell*, field&)
{
}

void onephase_f::wvel(lexer *p, fdm *a, ghostcell *pgc, field &w)
{
    for(int qn=0;qn<3;++qn)
    {
    WAIRLOOP
    wrk1(i,j,k) = w(i,j,k) - dt*(xphi(i,j,k)*ddwenox(w,a->w(i,j,k)) + yphi(i,j,k)*ddwenoy(w,a->w(i,j,k)) + zphi(i,j,k)*ddwenoz(w,a->w(i,j,k)) );
    
    pgc->start3(p,wrk1,gcval_w);

    WAIRLOOP
    w(i,j,k) = 0.5*w(i,j,k) + 0.5*wrk1(i,j,k) - 0.5*dt*(xphi(i,j,k)*ddwenox(wrk1,a->w(i,j,k)) + yphi(i,j,k)*ddwenoy(wrk1,a->w(i,j,k)) + zphi(i,j,k)*ddwenoz(wrk1,a->w(i,j,k)) );
    }
}