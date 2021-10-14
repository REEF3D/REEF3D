/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_exner.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"reinitopo.h"
#include"bedconc.h"
#include"topo_relax.h"
#include"sediment_fou.h"
#include"sediment_cds.h"
#include"sediment_wenoflux.h"
#include"sediment_weno_hj.h"
#include<math.h>

sediment_exner::sediment_exner(lexer* p, fdm *a, ghostcell* pgc, turbulence *pturb) :  bedshear(p,pturb), q0(p),dqx0(p),dqy0(p)
{
	if(p->S50==1)
	gcval_topo=151;

	if(p->S50==2)
	gcval_topo=152;

	if(p->S50==3)
	gcval_topo=153;
	
	if(p->S50==4)
	gcval_topo=154;
    
    
    rhosed=p->S22;
    rhowat=p->W1;
    g=9.81;
    d50=p->S20;
    ws=1.1*(rhosed/rhowat-1.0)*g*d50*d50;
    Ls = p->S20;
    
    
    pcb = new bedconc(p, pturb);
    
    prelax = new topo_relax(p);
    
    if(p->S32==1)
    pdx = new sediment_fou(p);
    
    if(p->S32==2)
    pdx = new sediment_cds(p);
    
    if(p->S32==4)
    pdx = new sediment_wenoflux(p);
    
    if(p->S32==5)
    pdx = new sediment_weno_hj(p);
}

sediment_exner::~sediment_exner()
{
}

void sediment_exner::start(fdm* a,lexer* p, convection* pconvec, ghostcell* pgc,reinitopo* preto, sediment_fdm *s)
{   
    //non_equillibrium_solve(p,a,pgc); 
   
    SLICELOOP4
    {
		topovel(p,a,pgc,vx,vy,vz);
        dqx0(i,j) = vx;
        dqy0(i,j) = vy;
		s->vz(i,j) = vz;
	}
    
	pgc->gcsl_start4(p,s->vz,1);
	
    timestep(p,a,pgc,s);
    
    
	
	SLICELOOP4
    a->bedzh(i,j) += p->dtsed*s->vz(i,j);

	pgc->gcsl_start4(p,a->bedzh,1);
}


















