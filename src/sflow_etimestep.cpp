/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"sflow_etimestep.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

#define HXIJ (fabs(b->hx(i,j))>wd_criterion?b->hx(i,j):wd_criterion)
#define HYIJ (fabs(b->hy(i,j))>wd_criterion?b->hy(i,j):wd_criterion)

sflow_etimestep::sflow_etimestep(lexer *p, fdm2D *b)
{
    wd_criterion=0.00005;
    
    if(p->A244==1)
    wd_criterion=p->A244_val;
    
    if(p->A245==1)
    wd_criterion=p->A245_val*p->DXM;
}

sflow_etimestep::~sflow_etimestep()
{
}

void sflow_etimestep::start(lexer *p, fdm2D* b, ghostcell* pgc)
{	
	p->umax=p->vmax=p->viscmax=0.0;
	p->dt_old=p->dt;
	double depthmax=-10.0;


// maximum velocities


	SLICELOOP1
	p->umax=MAX(p->umax,fabs(b->P(i,j)));

	p->umax=pgc->globalmax(p->umax);


	SLICELOOP2
	p->vmax=MAX(p->vmax,fabs(b->Q(i,j)));

	p->vmax=pgc->globalmax(p->vmax);

	SLICELOOP4
	depthmax=MAX(depthmax,b->depth(i,j));
    
    depthmax=MAX(depthmax,0.00001);
	
	depthmax=pgc->globalmax(depthmax);
    
    
	
    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(3)<<p->umax<<" \t utime: "<<p->utime<<endl;
	cout<<"vmax: "<<setprecision(3)<<p->vmax<<" \t vtime: "<<p->vtime<<endl;
	cout<<"fsftime: "<<p->lsmtime<<endl;
    }

    velmax=MAX(p->umax,p->vmax);


// 
    
	cu=2.0/((p->umax/p->DXM)+sqrt(4.0*fabs(b->maxF)/p->DXM));
	cv=2.0/((p->vmax/p->DXM)+sqrt(4.0*fabs(b->maxG)/p->DXM));
    
	if(p->A219==1)
    {
    cu=MIN(cu,2.0/((p->umax+sqrt(9.81*depthmax))/p->DXM));
	cv=MIN(cv,2.0/((p->vmax+sqrt(9.81*depthmax))/p->DXM));
    }
    
    if(p->A219==2)
    {
    cu=p->DXM/(fabs(p->umax)>1.0e-20?p->umax:1.0);
	cv=p->DXM/(fabs(p->vmax)>1.0e-20?p->vmax:1.0);
    }
    
    if(p->A219==3)
    {
    cu=2.0/((p->umax+sqrt(9.81*depthmax))/p->DXM);
	cv=2.0/((p->vmax+sqrt(9.81*depthmax))/p->DXM);
    }
    

	p->dt=p->N47*MIN(cu,cv);
	p->dt=pgc->timesync(p->dt);


	b->maxF=0.0;
	b->maxG=0.0;
}


void sflow_etimestep::ini(lexer *p, fdm2D* b, ghostcell* pgc)
{	
	
    p->umax=p->W10;
    
	SLICELOOP1
	p->umax=MAX(p->umax,fabs(b->P(i,j)));
    

	p->umax=pgc->globalmax(p->umax);


	SLICELOOP2
	p->vmax=MAX(p->vmax,fabs(b->Q(i,j)));

	p->vmax=pgc->globalmax(p->vmax);
	
	p->umax=MAX(p->umax,p->vmax);
	p->umax=MAX(p->umax,5.0);
	

	
	cu=2.0/((p->umax/p->DXM));
	
	
	
	p->dt=p->N47*cu;
	p->dt=pgc->timesync(p->dt);

	p->dt_old=p->dt;

	b->maxF=0.0;
	b->maxG=0.0;
}
