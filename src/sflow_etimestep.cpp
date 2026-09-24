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

#include"sflow_etimestep.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"


sflow_etimestep::sflow_etimestep(lexer *p, fdm2D *b)
{
    wd_criterion=0.00005;
    
    wd_criterion=p->A244;
}

sflow_etimestep::~sflow_etimestep()
{
}

void sflow_etimestep::start(lexer *p, fdm2D* b, ghostcell* pgc)
{	
    const double g = fabs(p->W22);
    double dx,dy,c,cmin;
    
	p->umax=p->vmax=p->viscmax=0.0;
	p->dt_old=p->dt;

// maximum velocities
	SLICELOOP4
    WETDRY
    {
	p->umax=MAX(p->umax,fabs(b->U(i,j)));
	p->vmax=MAX(p->vmax,fabs(b->V(i,j)));
    }

	p->umax=pgc->globalmax(p->umax);
	p->vmax=pgc->globalmax(p->vmax);
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(3)<<p->umax<<" \t utime: "<<p->utime<<endl;
	cout<<"vmax: "<<setprecision(3)<<p->vmax<<" \t vtime: "<<p->vtime<<endl;
	cout<<"fsftime: "<<p->lsmtime<<endl;
    }

// CFL: dt = N47 * 2*min( dx/(|u|+c), dy/(|v|+c) ),  c = sqrt(g h)
// (factor 2 as in the previous SFLOW time step definition, N 47 0.2 -> Courant number 0.4)
    cmin=1.0e20;
    
    SLICELOOP4
    WETDRY
    {
    c = sqrt(g*MAX(b->WL(i,j),wd_criterion));
    
    dx = p->DXN[IP]/(fabs(b->U(i,j)) + c);
    cmin = MIN(cmin,dx);
    
    if(p->j_dir==1)
    {
    dy = p->DYN[JP]/(fabs(b->V(i,j)) + c);
    cmin = MIN(cmin,dy);
    }
    
        // advection-only limit (A 219 2)
        if(p->A219==2)
        {
        dx = p->DXN[IP]/(fabs(b->U(i,j))>1.0e-20?fabs(b->U(i,j)):1.0e-20);
        cmin = MIN(cmin,dx);
        }
    }
    
    cmin = pgc->globalmin(cmin);
    
    if(cmin>1.0e19)
    cmin = p->DXM/sqrt(g*MAX(p->wd,wd_criterion));

	p->dt=p->N47*2.0*cmin;
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
    
    // include the shallow water wave speed
    cu=MIN(cu,p->DXM/(p->umax + sqrt(fabs(p->W22)*MAX(MAX(p->F60,p->wd),1.0))));
	
	
	
	p->dt=p->N47*cu;
	p->dt=pgc->timesync(p->dt);

	p->dt_old=p->dt;

	b->maxF=0.0;
	b->maxG=0.0;
}
