/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"fnpf_timestep.h"
#include<iomanip>
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

fnpf_timestep::fnpf_timestep(lexer* p):epsi(1.0e-3),maxtimestep(p->N49),c0_orig(p->N47)
{
}

fnpf_timestep::~fnpf_timestep()
{
}

void fnpf_timestep::start(fdm_fnpf *c, lexer *p,ghostcell *pgc)
{
    double depthmax=0.0;


    p->umax=p->vmax=p->wmax=p->viscmax=irsm=jrsm=krsm=0.0;
    p->epsmax=p->kinmax=p->pressmax=0.0;
	p->dt_old=p->dt;

	p->umax=p->vmax=p->wmax=p->viscmax=0.0;
	sqd=1.0/(p->dx*p->dx);

// maximum velocities

    SLICELOOP4
    FPWDCHECK
	depthmax=MAX(depthmax,c->depth(i,j));
	
	depthmax=pgc->globalmax(depthmax);

	LOOP
    FPWDCHECK
	p->umax=MAX(p->umax,fabs(c->u(i,j,k)));

	p->umax=pgc->globalmax(p->umax);


	LOOP
    FPWDCHECK
	p->vmax=MAX(p->vmax,fabs(c->v(i,j,k)));

	p->vmax=pgc->globalmax(p->vmax);


	LOOP
    FPWDCHECK
	p->wmax=MAX(p->wmax,fabs(c->w(i,j,k)));

	p->wmax=pgc->globalmax(p->wmax);
	

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(3)<<p->umax<<endl;
	cout<<"vmax: "<<setprecision(3)<<p->vmax<<endl;
	cout<<"wmax: "<<setprecision(3)<<p->wmax<<endl;
    }
	
	p->umax=MAX(p->umax,p->ufbmax);
	p->vmax=MAX(p->vmax,p->vfbmax);
	p->wmax=MAX(p->wmax,p->wfbmax);


    cu=cv=cw=1.0e10;
    
    
    SLICELOOP4
    FPWDCHECK
    p->wmax = MAX(fabs(c->Fz(i,j)),p->wmax);
    
    p->wmax=pgc->globalmax(p->wmax);
    
    //p->wmax = MAX3(p->wmax,p->umax,p->vmax);
        
    FLOOP
    FPWDCHECK
    {
    if(p->y_dir==1)
    dx = MIN(p->DXN[IP],p->DYN[JP]);
    
    if(p->y_dir==0)
    dx = p->DXN[IP];
    
    cu = MIN(cu, 1.0/((fabs(MAX(p->umax, 1.0*sqrt(9.81*depthmax)))/dx)));
    cv = MIN(cv, 1.0/((fabs(MAX(p->vmax, 1.0*sqrt(9.81*depthmax)))/dx)));
    //cw = MIN(cw, 1.0/((fabs(p->wmax)/dx)));
    }

    cw = MIN3(cu,cv,cw);
    
   	p->dt=p->N47*cw;
    
    //p->dt=MIN(1.1*p->dt_old,p->dt);

    
	p->dt=pgc->timesync(p->dt);

	p->dt=MIN(p->dt,maxtimestep);
	p->turbtimestep=p->dt;

}

void fnpf_timestep::ini(fdm_fnpf* c, lexer* p,ghostcell* pgc)
{
	p->umax=p->vmax=p->wmax=p->viscmax=-1e19;
	
    LOOP
	p->umax=MAX(p->umax,fabs(c->u(i,j,k)));

	p->umax=pgc->globalmax(p->umax);


	LOOP
	p->vmax=MAX(p->vmax,fabs(c->v(i,j,k)));

	p->vmax=pgc->globalmax(p->vmax);


	LOOP
	p->wmax=MAX(p->wmax,fabs(c->w(i,j,k)));

	p->wmax=pgc->globalmax(p->wmax);
	
	p->umax=MAX(p->umax,2.0*p->ufbmax);
	p->umax=MAX(p->umax,2.0*p->vfbmax);
	p->umax=MAX(p->umax,2.0*p->wfbmax);
    
    p->umax=MAX(p->umax,2.0*p->X210_u);
	p->umax=MAX(p->umax,2.0*p->X210_v);
	p->umax=MAX(p->umax,2.0*p->X210_w);

	p->dt=p->dx/(p->umax+epsi);
    

    p->umax+=10.0;

	cu= + 2.0/((p->umax/p->dx)+sqrt(pow(p->umax/p->dx,2.0)+(4.0*sqrt(fabs(c->gi) + fabs(c->gj) +fabs(c->gk)))/p->dx));


	p->dt=p->N47*cu;

    
	p->dt=pgc->timesync(p->dt);
	p->dt_old=p->dt;

	p->maxkappa=0.0;
}

