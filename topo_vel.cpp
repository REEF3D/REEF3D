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

#include"topo_vel.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bedconc.h"
#include"topo_relax.h"
#include"fnpf_weno.h"

topo_vel::topo_vel(lexer* p, turbulence *pturb): ccipol(p), norm_vec(p), dx(p->dx),epsi(1.6*p->dx)
{
    rhosed=p->S22;
    rhowat=p->W1;
    g=9.81;
    d50=p->S20;
    ws=1.1*(rhosed/rhowat-1.0)*g*d50*d50;

    pcb = new bedconc(p, pturb);
    
    prelax = new topo_relax(p);
    
    pdx = new fnpf_weno(p);
}

topo_vel::~topo_vel()
{
}

void topo_vel::topovel(lexer* p,fdm* a, ghostcell *pgc, double& vx, double& vy, double& vz)
{
	double uvel,vvel,u_abs;
	double signx,signy;
	double dqx,dqy;
	
	vx=0.0;
	vy=0.0;
	vz=0.0;
	 
	if(p->pos_x()>=p->S71 && p->pos_x()<=p->S72)
	{						
        pip=1;
        uvel=0.5*(a->P(i,j)+a->P(i-1,j));
        pip=0;

        pip=2;
        vvel=0.5*(a->Q(i,j)+a->Q(i,j-1));
        pip=0;
		
		u_abs = sqrt(uvel*uvel + vvel*vvel);
		signx=fabs(u_abs)>1.0e-10?uvel/fabs(u_abs):0.0;
		signy=fabs(u_abs)>1.0e-10?vvel/fabs(u_abs):0.0;
	
		
	dqx=dqy=0.0;

    //dqx = (a->bedload(i+1,j)-a->bedload(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);
    //dqy = (a->bedload(i,j+1)-a->bedload(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);
    
    dqx = pdx->sx(p,a->bedload,signx);
    dqy = pdx->sy(p,a->bedload,signy);
		
	// Exner equations
    vz =  -prelax->rf(p,a,pgc)*(1.0/(1.0-p->S24))*(dqx*signx + dqy*signy) + ws*(a->conc(i,j,k) - pcb->cbed(p,a,pgc,a->topo)); 
	}
}

void topo_vel::rkio_update(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
	int n;
		for(n=0;n<p->gcout4a_count;++n)
        {
        i=p->gcout4a[n][0];
        j=p->gcout4a[n][1];
        k=p->gcout4a[n][2];

        f(i+1,j,k)=a->topo(i+1,j,k);
        f(i+2,j,k)=a->topo(i+2,j,k);
        f(i+3,j,k)=a->topo(i+3,j,k);
        }
		
		for(n=0;n<p->gcin4a_count;++n)
        {
        i=p->gcin4a[n][0];
        j=p->gcin4a[n][1];
        k=p->gcin4a[n][2];

        f(i-1,j,k)=a->topo(i-1,j,k);
        f(i-2,j,k)=a->topo(i-2,j,k);
        f(i-3,j,k)=a->topo(i-3,j,k);
        }
}

