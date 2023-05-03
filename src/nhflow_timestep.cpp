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

#include"nhflow_timestep.h"
#include<iomanip>
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

nhflow_timestep::nhflow_timestep(lexer* p):epsi(1.0e-3),maxtimestep(p->N49),c0_orig(p->N47)
{
}

nhflow_timestep::~nhflow_timestep()
{
}

void nhflow_timestep::start(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    double depthmax=0.0;


    p->umax=p->vmax=p->wmax=p->viscmax=irsm=jrsm=krsm=0.0;
    p->epsmax=p->kinmax=p->pressmax=0.0;
	p->dt_old=p->dt;

	p->umax=p->vmax=p->wmax=p->viscmax=0.0;
	sqd=1.0/(p->DXM*p->DXM);

// maximum velocities

    SLICELOOP4
	depthmax=MAX(depthmax,d->depth(i,j));
	
	depthmax=pgc->globalmax(depthmax);

	LOOP
	p->umax=MAX(p->umax,fabs(d->U[IJK]));

	p->umax=pgc->globalmax(p->umax);


	LOOP
	p->vmax=MAX(p->vmax,fabs(d->U[IJK]));

	p->vmax=pgc->globalmax(p->vmax);


	LOOP
	p->wmax=MAX(p->wmax,fabs(d->W[IJK]));

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
    
    
    // visc
    /*SLICELOOP4
	p->viscmax=MAX(p->viscmax, d->vb(i,j));

	p->viscmax=pgc->globalmax(p->viscmax);
    
    if(p->mpirank==0 && (p->count%p->P12==0) && p->viscmax>0.0)
	cout<<"viscmax: "<<p->viscmax<<endl;*/
    

    LOOP
    {
    if(p->j_dir==1 && p->knoy>1)
    dx = MIN(p->DXN[IP],p->DYN[JP]);
    
    if(p->j_dir==0 || p->knoy==1)
    dx = p->DXN[IP];
    
    cu = MIN(cu, 1.0/((fabs(MAX(p->umax, sqrt(9.81*depthmax)))/dx)));
    
    if(p->j_dir==1 )
    cv = MIN(cv, 1.0/((fabs(MAX(p->vmax, sqrt(9.81*depthmax)))/dx)));
    
    }

    if(p->j_dir==1 )
    cu = MIN(cu,cv);
    
   	p->dt=p->N47*cu;
    
	p->dt=pgc->timesync(p->dt);
    
    if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"dt: "<<p->dt<<endl;

    
    if (p->N48==0) 
    p->dt=maxtimestep;
    
    else
	p->dt=MIN(p->dt,maxtimestep);
    
    
}

void nhflow_timestep::ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    double depthmax;
    
	p->umax = p->vmax = p->wmax = -1e19;
    depthmax = -1e19;
    
    cu=cv=1.0e10;
    
    
    SLICELOOP4
	depthmax=MAX(depthmax,p->wd - d->bed(i,j));
	
	depthmax=pgc->globalmax(depthmax);
	
    LOOP
	p->umax=MAX(p->umax,fabs(d->U[IJK]));

	p->umax=pgc->globalmax(p->umax);


	LOOP
	p->vmax=MAX(p->vmax,fabs(d->U[IJK]));

	p->vmax=pgc->globalmax(p->vmax);


	LOOP
	p->wmax=MAX(p->wmax,fabs(d->W[IJK]));

	p->wmax=pgc->globalmax(p->wmax);
    
	
	p->umax=MAX(p->umax,2.0*p->ufbmax);
	p->umax=MAX(p->umax,2.0*p->vfbmax);
	p->umax=MAX(p->umax,2.0*p->wfbmax);
    
    p->umax=MAX(p->umax,2.0*p->X210_u);
	p->umax=MAX(p->umax,2.0*p->X210_v);
	p->umax=MAX(p->umax,2.0*p->X210_w);

	p->dt=p->DXM/(p->umax+epsi);
    

    p->umax+=10.0;


    LOOP
    {
    if(p->j_dir==1 && p->knoy>1)
    dx = MIN(p->DXN[IP],p->DYN[JP]);
    
    if(p->j_dir==0 || p->knoy==1)
    dx = p->DXN[IP];
    
    //if(p->mpirank==0)
    //cout<<"dx: "<<dx<<" p->DXN[IP]: "<<p->DXN[IP]<<" p->DYN[JP]: "<<p->DYN[JP]<<endl;
    
    cu = MIN(cu, 1.0/((fabs(MAX(p->umax, sqrt(9.81*depthmax)))/dx)));
    cv = MIN(cv, 1.0/((fabs(MAX(p->vmax, sqrt(9.81*depthmax)))/dx)));
    }

	cu = MIN(cu,cv);
    
   	p->dt=p->N47*cu;
    
	//p->dt=pgc->timesync(p->dt);
    p->dt=pgc->globalmin(p->dt);
	p->dt_old=p->dt;

    
    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(3)<<p->umax<<endl;
	cout<<"vmax: "<<setprecision(3)<<p->vmax<<endl;
	cout<<"wmax: "<<setprecision(3)<<p->wmax<<endl;
    cout<<"dmax: "<<setprecision(3)<<depthmax<<endl;
    }
}

