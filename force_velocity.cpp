/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"force.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void force::velocity(lexer* p, fdm *a, ghostcell *pgc)
{
    double timestep = p->simtime/double(MAX(p->count,1));
	int kvert = ke-ks;
	
	double uvel2,vvel2,uvel1,vvel1,uvel0,vvel0,zcoor;
	double wsf_height0,wsf_height1,wsf_height2;
	double zloc;
	double dh = 0.05*p->dx;
	
	uvel0=vvel0=uvel1=vvel1=uvel2=vvel2=0.0;
    
    i=0;
    j=0;
	
	wsf_height0 = wave_h(p,pgc,xm,0.0,0.0);
	wsf_height1 = wave_h(p,pgc,xm,0.0,0.0);
	wsf_height2 = wave_h(p,pgc,xm,0.0,0.0);
	
	
	kvert = p->conv((wsf_height0-p->originz)/dh)+1;
	
	int count0=0;
	int count1=0;
	int count2=0;
	
	for(int kk=0;kk<kvert;++kk)
	{
		zcoor = -(p->phimean-(double(kk) * dh + p->originz)-0.5*dh);
		zloc  =  (double(kk)+0.5)*dh+p->originz;
		
		p->simtime-=p->dt;
		if(zloc<wsf_height0)
		{
		uvel0 += wave_u(p,pgc,xm,0.0,zcoor);
		vvel0 += wave_v(p,pgc,xm,0.0,zcoor);
		++count0;
		}
		p->simtime+=p->dt;
		if(zloc<wsf_height1)
		{
		uvel1 += wave_u(p,pgc,xm,0.0,zcoor);
		vvel1 += wave_v(p,pgc,xm,0.0,zcoor);
		++count1;
		}
		p->simtime+=p->dt;
		if(zloc<wsf_height2)
		{
		uvel2 += wave_u(p,pgc,xm,0.0,zcoor);
		vvel2 += wave_v(p,pgc,xm,0.0,zcoor);
		++count2;
		}
		p->simtime-=p->dt;
	}
	
	
	uvel0 = uvel0 /double(count0);
	vvel0 = vvel0 /double(count0);
	
	uvel1 = uvel1 /double(count1);
	vvel1 = vvel1 /double(count1);
	
	uvel2 = uvel2 /double(count2);
	vvel2 = vvel2 /double(count2);
	
	uinf = uvel2;
	vinf = vvel2;
	
	//uinf_dt = (3.0*uvel2 - 4.0*uvel1 + uvel0)/(2.0*timestep);
	//vinf_dt = (3.0*vvel2 - 4.0*vvel1 + vvel0)/(2.0*timestep);
	
	uinf_dt = (uvel2-uvel1)/p->dt;
	vinf_dt = (vvel2-vvel1)/p->dt;
}



