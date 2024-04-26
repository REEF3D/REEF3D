/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"momentum_forcing.h"
#include"6DOF.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"FSI.h"

momentum_forcing::momentum_forcing(lexer* p)
{
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
}

momentum_forcing::~momentum_forcing()
{
}

void momentum_forcing::momentum_forcing_start(fdm* a, lexer* p, ghostcell *pgc, sixdof* p6dof, vrans* pvrans, vector<net*>& pnet, fsi* pfsi,
                                              field &u, field &v, field &w, field &fx, field &fy, field &fz, int iter, double alpha, bool final)
{
    if(p->X10==1 || p->Z10>0 || p->G3==1)
    {
	starttime=pgc->timer();
    
        // Forcing
        ULOOP
        fx(i,j,k) = 0.0;
       
        VLOOP
        fy(i,j,k) = 0.0;
      
        WLOOP
        fz(i,j,k) = 0.0;
        
        pgc->start1(p,fx,10);
        pgc->start2(p,fy,11);
        pgc->start3(p,fz,12); 
         

        if(p->G3==1)
        pgc->solid_forcing(p,a,alpha,u,v,w,fx,fy,fz);         
        
        p6dof->start_twoway(p,a,pgc,pvrans,pnet,iter,u,v,w,fx,fy,fz,final);
        
        pfsi->forcing(p,a,pgc,alpha,u,v,w,fx,fy,fz,final);
 
        ULOOP
        {
        u(i,j,k) += alpha*p->dt*CPOR1*fx(i,j,k);
        
        if(p->count<10)
        a->maxF = MAX(fabs(alpha*CPOR1*fx(i,j,k)), a->maxF);
        
        p->fbmax = MAX(fabs(alpha*CPOR1*fx(i,j,k)), p->fbmax);
        }
        
        VLOOP
        {
        v(i,j,k) += alpha*p->dt*CPOR2*fy(i,j,k);
        
        if(p->count<10)
        a->maxG = MAX(fabs(alpha*CPOR2*fy(i,j,k)), a->maxG);
        
        p->fbmax = MAX(fabs(alpha*CPOR2*fy(i,j,k)), p->fbmax);
        }
        
        WLOOP
        {
        w(i,j,k) += alpha*p->dt*CPOR3*fz(i,j,k);
        
        if(p->count<10)
        a->maxH = MAX(fabs(alpha*CPOR3*fz(i,j,k)), a->maxH);
        
        p->fbmax = MAX(fabs(alpha*CPOR3*fz(i,j,k)), p->fbmax);
        }

        pgc->start1(p,u,gcval_u);
        pgc->start2(p,v,gcval_v);
        pgc->start3(p,w,gcval_w);
        
        p->fbtime+=pgc->timer()-starttime;
    }
}


