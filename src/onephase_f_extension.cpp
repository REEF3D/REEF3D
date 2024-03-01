/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
    FLUIDLOOP
    {
    nx0 = -(a->phi(i+1,j,k) - a->phi(i-1,j,k))/(p->DXP[IP]+p->DXP[IM1]);
    ny0 = -(a->phi(i,j+1,k) - a->phi(i,j-1,k))/(p->DYP[JP]+p->DYP[JM1]);
    nz0 = -(a->phi(i,j,k+1) - a->phi(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]);
    
    nx = ddwenox(a->phi,nx0);
    ny = ddwenoy(a->phi,ny0);
    nz = ddwenoz(a->phi,nz0);

        dnorm = sqrt(nx*nx + ny*ny + nz*nz);
        
        nx/=(dnorm>1.0e-8?dnorm:1.0e10);
        ny/=(dnorm>1.0e-8?dnorm:1.0e10);
        nz/=(dnorm>1.0e-8?dnorm:1.0e10);
        
    //a->test(i,j,k) = -nx;
    xphi(i,j,k) = -nx;
    yphi(i,j,k) = -ny;
    zphi(i,j,k) = -nz;
    }
    
    for(int qn=0;qn<2;++qn)
    {
    UFLUIDLOOP
    urk1(i,j,k) = u(i,j,k) 
                - dt*(xphi(i,j,k)*ddwenox(u,xphi(i,j,k)) 
                    + yphi(i,j,k)*ddwenoy(u,yphi(i,j,k)) 
                    + zphi(i,j,k)*ddwenoz(u,zphi(i,j,k)));
    
    pgc->start1(p,urk1,gcval_u);

    UFLUIDLOOP
    uf(i,j,k) = 0.5*u(i,j,k) + 0.5*urk1(i,j,k) 
             - 0.5*dt*(xphi(i,j,k)*ddwenox(urk1,xphi(i,j,k)) 
                     + yphi(i,j,k)*ddwenoy(urk1,yphi(i,j,k)) 
                     + zphi(i,j,k)*ddwenoz(urk1,zphi(i,j,k)));
                     
    UAIRLOOP
    u(i,j,k)=uf(i,j,k);
    }
}

void onephase_f::vvel(lexer *p, fdm *a, ghostcell*, field&)
{
}

void onephase_f::wvel(lexer *p, fdm *a, ghostcell *pgc, field &w)
{
    for(int qn=0;qn<2;++qn)
    {
    WFLUIDLOOP
   a->test(i,j,k) =   wrk1(i,j,k) = w(i,j,k) 
                - dt*(xphi(i,j,k)*ddwenox(w,xphi(i,j,k)) 
                    + yphi(i,j,k)*ddwenoy(w,yphi(i,j,k)) 
                    + zphi(i,j,k)*ddwenoz(w,zphi(i,j,k)));
    
    pgc->start3(p,wrk1,gcval_w);

    WFLUIDLOOP
    wf(i,j,k) = 0.5*w(i,j,k) + 0.5*wrk1(i,j,k) 
             - 0.5*dt*(xphi(i,j,k)*ddwenox(wrk1,xphi(i,j,k)) 
                     + yphi(i,j,k)*ddwenoy(wrk1,yphi(i,j,k)) 
                     + zphi(i,j,k)*ddwenoz(wrk1,zphi(i,j,k)));
    
    WAIRLOOP
    w(i,j,k)=wf(i,j,k);
    }
}