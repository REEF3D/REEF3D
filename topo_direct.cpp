/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"topo_direct.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinitopo.h"
#include"slice4.h"

topo_direct::topo_direct(lexer* p, fdm *a, ghostcell* pgc, turbulence *pturb) : topo_vel(p,pturb)
{
	if(p->S50==1)
	gcval_topo=151;

	if(p->S50==2)
	gcval_topo=152;

	if(p->S50==3)
	gcval_topo=153;
	
	if(p->S50==4)
	gcval_topo=154;
}

topo_direct::~topo_direct()
{
}

void topo_direct::start(fdm* a,lexer* p, discrete* pdisc, ghostcell* pgc,reinitopo* preto)
{
    starttime=pgc->timer();
	
	slice4 dh(p),dh1(p),dh2(p);
	
	double maxdh;
    
	if(p->S31==1)
    SLICELOOP4
    {
		topovel(p,a,pgc,vx,vy,vz);
		
		dh(i,j) = vz;
	}
    
    if(p->S31==2)
    SLICELOOP4
    {
		topovel_xy_cds(p,a,pgc,vx,vy,vz);
		
		dh(i,j) = vz;
	}
    
    LOOP
    a->test(i,j,k) = dh(i,j);
    
    pgc->start4(p,a->test,1);
	
	pgc->gcsl_start4(p,dh,1);
	
	if(p->S39==2)
	SLICELOOP4
	dh(i,j) = 0.5*(3.0*dh(i,j) - dh1(i,j));

	if(p->S39==3)
	SLICELOOP4
	dh(i,j) = (1.0/12.0)*(23.0*dh(i,j) - 16.0*dh1(i,j) + 5.0*dh2(i,j));

	
	SLICELOOP4
    {
    dh2(i,j)=dh1(i,j);
    dh1(i,j)=dh(i,j);
    }
	
	maxdh=0.0;
	SLICELOOP4
	maxdh = MAX(fabs(dh(i,j)),maxdh);	

	
	double localmaxdh = maxdh;
	maxdh=pgc->globalmax(maxdh);
	
	if(p->S15==0)
    p->dtsed=MIN(p->S13, (p->S14*p->DXM)/(fabs(maxdh)>1.0e-15?maxdh:1.0e-15));

    if(p->S15==1)
    p->dtsed=MIN(p->dt, (p->S14*p->DXM)/(fabs(maxdh)>1.0e-15?maxdh:1.0e-15));
    
    if(p->S15==2)
    p->dtsed=p->S13;

    p->dtsed=pgc->timesync(p->dtsed);
	
	if(p->mpirank==0)
	cout<<p->mpirank<<" maxdh: "<<setprecision(4)<<maxdh<<" dtsed: "<<setprecision(4)<<p->dtsed<<endl;
	
	SLICELOOP4
    a->bedzh(i,j) += p->dtsed*dh(i,j);


	pgc->gcsl_start4(p,a->bedzh,1);
    pgc->dgcslpol(p,a->bedzh,p->dgcsl4,p->dgcsl4_count,14);
    a->bedzh.ggcpol(p);

	p->topotime=pgc->timer()-starttime;
}




