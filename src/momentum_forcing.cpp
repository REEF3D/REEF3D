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

#include"momentum_forcing.h"
#include<algorithm>
#include"6DOF.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"FSI.h"
#include"rodtree_coupling.h"
#include"fem_coupling.h"
#include"dem.h"
#include"sediment.h"

momentum_forcing::momentum_forcing(lexer* p) : prodtree(nullptr), pfem(nullptr)
{
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
}

momentum_forcing::~momentum_forcing()
{
    delete prodtree;
    delete pfem;
}

void momentum_forcing::momentum_forcing_start(fdm* a, lexer* p, ghostcell *pgc, sixdof* p6dof, fsi* pfsi,
                                              field &u, field &v, field &w, field &fx, field &fy, field &fz, int iter, double alpha, bool final)
{

	starttime=pgc->timer();
    
        // Forcing: zero everywhere, halos included (only the interior is applied to u/v/w below,
        // the forcing modules that spread into fx/fy/fz update the halos themselves)
        std::fill(fx.V, fx.V + p->imax*p->jmax*p->kmax, 0.0);
        std::fill(fy.V, fy.V + p->imax*p->jmax*p->kmax, 0.0);
        std::fill(fz.V, fz.V + p->imax*p->jmax*p->kmax, 0.0);
         
        pgc->solid_forcing(p,a,alpha,u,v,w,fx,fy,fz);         
        
        p6dof->start_cfd(p,a,pgc,iter,u,v,w,fx,fy,fz,final);
        
        pfsi->forcing(p,a,pgc,alpha,u,v,w,fx,fy,fz,final);
        
        // flexible rod trees: sample, spread reaction, advance structure (final stage)
        if(p->Z20>0)
        {
            if(prodtree==nullptr)
            prodtree = new rodtree_coupling(p,pgc);
            
            prodtree->start_cfd(p,a,pgc,alpha,u,v,w,fx,fy,fz,final);
        }

        // FEM solid structures: surface direct forcing, fluid loads, advance solid (final stage)
        if(p->Z30>0)
        {
            if(pfem==nullptr)
            pfem = new fem_coupling(p,pgc);

            pfem->start_cfd(p,a,pgc,alpha,u,v,w,fx,fy,fz,final);
        }
 
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
        
        p->fbtime+=pgc->timer()-starttime;

        // DEM: resolved particle forcing and unresolved momentum source
        if(pdem!=nullptr)
        pdem->forcing_cfd(p,a,pgc,iter,alpha,u,v,w,final);
        
        // particle sediment (CPM): two-way coupling
        if(psed!=nullptr)
        psed->forcing_cfd(p,a,pgc,alpha,u,v,w);
        
        
    // ghostcell update: the flags depend on the geometry only (solid, topo, fb), rebuild when it changed
    if(pgc->geometry_changed(p,a))
    {
    pgc->solid_forcing_flag_update(p,a);
    pgc->gcdf_update(p,a);
    pgc->gcb_velflagio(p,a);
    }
}


