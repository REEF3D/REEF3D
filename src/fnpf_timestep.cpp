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
	sqd=1.0/(p->DXM*p->DXM);

// maximum velocities

    FFILOOP4
    FPWDCHECK
	depthmax=MAX(depthmax,c->depth(i,j));
	
	depthmax=pgc->globalmax(depthmax);

	FFILOOP4
    FPWDCHECK
    {
	p->umax=MAX(p->umax,fabs(c->U[FIJK]));
    }

	p->umax=pgc->globalmax(p->umax);


	FFILOOP4
    FPWDCHECK
	p->vmax=MAX(p->vmax,fabs(c->V[FIJK]));

	p->vmax=pgc->globalmax(p->vmax);


	FFILOOP4
    FPWDCHECK
	p->wmax=MAX(p->wmax,fabs(c->W[FIJK]));

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
    
    
    FFILOOP4
    FPWDCHECK
    p->wmax = MAX(fabs(c->Fz(i,j)),p->wmax);
    
    p->wmax=pgc->globalmax(p->wmax);
    
    
    // visc
    SLICELOOP4
	p->viscmax=MAX(p->viscmax, c->vb(i,j));

	p->viscmax=pgc->globalmax(p->viscmax);
    
    if(p->mpirank==0 && (p->count%p->P12==0) && p->viscmax>0.0)
	cout<<"viscmax: "<<p->viscmax<<endl;
    

    FLOOP
    FPWDCHECK
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

void fnpf_timestep::ini(fdm_fnpf* c, lexer* p,ghostcell* pgc)
{
    double depthmax;
    
	p->umax = p->vmax = p->wmax = -1e19;
    depthmax = -1e19;
    
    cu=cv=1.0e10;
    
    
    FFILOOP4
    FPWDCHECK
	depthmax=MAX(depthmax,p->wd - c->bed(i,j));
	
	depthmax=pgc->globalmax(depthmax);
	
    k=p->knoz;
    SLICELOOP4
	p->umax=MAX(p->umax,fabs((c->Fifsf[Ip1J]-c->Fifsf[Im1J])/(p->DXP[IP]+p->DXP[IM1])));

	p->umax=pgc->globalmax(p->umax);


	k=p->knoz;
    SLICELOOP4
	p->vmax=MAX(p->vmax,fabs((c->Fifsf[IJp1]-c->Fifsf[IJm1])/(p->DYP[JP]+p->DYP[JM1])));

	p->vmax=pgc->globalmax(p->vmax);


	FLOOP
	p->wmax=MAX(p->wmax,fabs(c->W[FIJK]));

	p->wmax=pgc->globalmax(p->wmax);
	
	p->umax=MAX(p->umax,2.0*p->ufbmax);
	p->umax=MAX(p->umax,2.0*p->vfbmax);
	p->umax=MAX(p->umax,2.0*p->wfbmax);
    
    p->umax=MAX(p->umax,2.0*p->X210_u);
	p->umax=MAX(p->umax,2.0*p->X210_v);
	p->umax=MAX(p->umax,2.0*p->X210_w);

	p->dt=p->DXM/(p->umax+epsi);
    

    p->umax+=10.0;

 
    FLOOP
    FPWDCHECK
    {
    if(p->j_dir==1 && p->knoy>1)
    dx = MIN(p->DXN[IP],p->DYN[JP]);
    
    if(p->j_dir==0 || p->knoy==1)
    dx = p->DXN[IP];
    
    
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

