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
    
    // k
    k=p->knoz;

// maximum velocities
    SLICELOOP4
	depthmax=MAX(depthmax,c->depth(i,j));

	SLICELOOP4
    if(i+p->origin_i>=5 || p->j_dir==1)
    {
	// the surface node (k = knoz): of the 3D velocities, or the surface ones (fnpf_memory.h)
	p->umax=MAX(p->umax,fabs(c->U!=nullptr ? c->U[FIJK] : c->Us(i,j)));
    p->vmax=MAX(p->vmax,fabs(c->V!=nullptr ? c->V[FIJK] : c->Vs(i,j)));
    p->wmax=MAX(p->wmax,fabs(c->W!=nullptr ? c->W[FIJK] : c->Ws(i,j)));
    }
    
    SLICELOOP4
	p->viscmax=MAX(p->viscmax, c->vb(i,j));
    
    // one MPI_Allreduce instead of five
    double red[5] = {depthmax, p->umax, p->vmax, p->wmax, p->viscmax};
    pgc->globalmax(red,5);
    depthmax  = red[0];
    p->umax   = red[1];
    p->vmax   = red[2];
    p->wmax   = red[3];
    p->viscmax= red[4];

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(3)<<p->umax<<endl;
	cout<<"vmax: "<<setprecision(3)<<p->vmax<<endl;
	cout<<"wmax: "<<setprecision(3)<<p->wmax<<endl;
    }
	
    cu=cv=cw=1.0e10;
    
    
    

    if(p->mpirank==0 && (p->count%p->P12==0) && p->viscmax>0.0)
	cout<<"viscmax: "<<p->viscmax<<endl;
    

    SLICELOOP4
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
    
    // minimum over all ranks (was the value of rank 0, i.e. of its cells only)
    cu = pgc->globalmin(cu);
    
    // no water and no velocity: the time step would be N47*1e10 and the run would end silently
    if(cu>=1.0e10)
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! FNPF: no time step, the domain holds no water (max water depth "<<depthmax<<" m, umax "<<p->umax<<" m/s) !!!"<<endl<<endl;
        
        pgc->final(true);
    }

   	p->dt=p->N47*cu;
    
    if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"dt: "<<p->dt<<endl;

    
    if (p->N48==0) 
    p->dt=maxtimestep;
    
    else
	p->dt=MIN(p->dt,maxtimestep);
    
    
}

void fnpf_timestep::ini(fdm_fnpf* c, lexer* p,ghostcell* pgc)
{
    double depthmax,wdt;
    
	p->umax = p->vmax = p->wmax = -1e19;
    depthmax = -1e19;
    
    cu=cv=1.0e10;
    
    
    SLICELOOP4
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


	// (only printed; without the 3D velocities the surface ones, fnpf_memory.h)
	if(c->W!=nullptr)
	FLOOP
	p->wmax=MAX(p->wmax,fabs(c->W[FIJK]));
	
	if(c->W==nullptr)
	SLICELOOP4
	p->wmax=MAX(p->wmax,fabs(c->Ws(i,j)));

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
    
    if(p->B94==0)
	wdt=p->phimean;

	if(p->B94==1)
	wdt=p->B94_wdt;
    
    // all columns dry: no wet cell set the time step above, it would be N47*1e10 and the run
    // would end after one step without a message (e.g. a bed level above the still water level
    // everywhere from a broken geo.dat interpolation)
    int nwet=0, ncol=0;
    double bedmin=1.0e20, bedmax=-1.0e20;
    
    SLICELOOP4
    {
    ++ncol;
    
    if(p->wet[IJ]>0)
    ++nwet;
    
    bedmin=MIN(bedmin,c->bed(i,j));
    bedmax=MAX(bedmax,c->bed(i,j));
    }
    
    nwet=pgc->globalisum(nwet);
    ncol=pgc->globalisum(ncol);
    bedmin=pgc->globalmin(bedmin);
    bedmax=pgc->globalmax(bedmax);
    
    if(nwet==0)
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! FNPF: all "<<ncol<<" columns are dry: bed level "<<bedmin<<" to "<<bedmax
            <<" m, still water level F 60 = "<<p->phimean<<" m, wetting criterion A 344 = "<<p->A344<<" m."
            <<" Check the bed level of the grid (DIVEMesh geo.dat, G 51/G 52, G 23) and F 60 !!!"<<endl<<endl;
        
        pgc->final(true);
    }
    
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(3)<<p->umax<<endl;
	cout<<"vmax: "<<setprecision(3)<<p->vmax<<endl;
	cout<<"wmax: "<<setprecision(3)<<p->wmax<<endl;
    cout<<"depth_max: "<<setprecision(3)<<depthmax<<" |  wave depth: "<<wdt<<endl;
    }
}

