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

#include"pftimestep.h"
#include<iomanip>
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"

pftimestep::pftimestep(lexer* p):epsi(1.0e-19),maxtimestep(p->N49),c0_orig(p->N47)
{
}

pftimestep::~pftimestep()
{
}

void pftimestep::start(fdm *a, lexer *p,ghostcell *pgc, turbulence *pturb)
{
    p->umax=p->vmax=p->wmax=p->viscmax=irsm=jrsm=krsm=0.0;
    p->epsmax=p->kinmax=p->pressmax=0.0;
	p->dt_old=p->dt;
    depthmax=0.0;

	p->umax=p->vmax=p->wmax=p->viscmax=0.0;
	sqd=1.0/(p->DXM*p->DXM);
    
    FILOOP4
	depthmax=MAX(depthmax,a->depth(i,j));
	
	depthmax=pgc->globalmax(depthmax);

// maximum velocities

	ULOOP
	p->umax=MAX(p->umax,fabs(a->u(i,j,k)));

	p->umax=pgc->globalmax(p->umax);


	VLOOP
	p->vmax=MAX(p->vmax,fabs(a->v(i,j,k)));

	p->vmax=pgc->globalmax(p->vmax);


	WLOOP
	p->wmax=MAX(p->wmax,fabs(a->w(i,j,k)));

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
    
    if(p->A10==4)
    {
    ULOOP
    {
    dx = MIN3(p->DXP[IP],p->DYN[JP],p->DZN[KP]);

	cu = MIN(cu, 2.0/((fabs(p->umax)/dx)));
    }

    VLOOP
    {
    dx = MIN3(p->DXN[IP],p->DYP[JP],p->DZN[KP]);
    
	cv = MIN(cv, 2.0/((fabs(p->vmax)/dx)));
    }
    
    WLOOP
    {
    dx = MIN3(p->DXN[IP],p->DYN[JP],p->DZP[KP]);

	cw = MIN(cw, 2.0/((fabs(p->wmax)/dx)));
    }
    }
    
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
    
    if(p->j_dir==1)
   	cu = MIN(cu,cv);
    
    
   	p->dt=p->N47*cu;
    
	p->dt=pgc->timesync(p->dt);

	p->dt=MIN(p->dt,maxtimestep);
    
    //cout<<p->mpirank<<" p->umax: "<<p->umax<<" depthmax: "<<depthmax<<" cu: "<<cu<<" cv: "<<cv<<" dt: "<<p->dt<<endl;
     
	
	a->maxF=0.0;
	a->maxG=0.0;
	a->maxH=0.0;
}

void pftimestep::ini(fdm* a, lexer* p,ghostcell* pgc)
{
    cu=cv=1.0e10;
    p->umax=p->vmax=0.0;
    
    
    depthmax=0.0;
    
    FILOOP4
	depthmax=MAX(depthmax,a->depth(i,j));
	
	depthmax=pgc->globalmax(depthmax);
    
    
    
	p->umax=p->vmax=p->wmax=p->viscmax=-1e19;
    

    SLICELOOP4
	p->umax=MAX(p->umax,fabs((a->Fifsf(i+1,j)-a->Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1])));

	p->umax=pgc->globalmax(p->umax);


	k=p->knoz;
    SLICELOOP4
	p->vmax=MAX(p->vmax,fabs((a->Fifsf(i+1,j+1)-a->Fifsf(i,j-1))/(p->DYP[JP]+p->DYP[JM1])));
	
	p->umax=MAX(p->umax,2.0*p->ufbmax);
	p->umax=MAX(p->umax,2.0*p->vfbmax);
	p->umax=MAX(p->umax,2.0*p->wfbmax);
    
    p->umax=MAX(p->umax,2.0*p->X210_u);
	p->umax=MAX(p->umax,2.0*p->X210_v);
	p->umax=MAX(p->umax,2.0*p->X210_w);
    
    p->umax+=10.0;

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
    
    
    
    if(p->j_dir==1) 
   	cu = MIN(cu,cv);
    
   	p->dt=p->N47*cu;
    
    //cout<<p->mpirank<<" p->umax: "<<p->umax<<" cu: "<<cu<<" cv: "<<cv<<" dt: "<<p->dt<<endl;
    
	p->dt=pgc->timesync(p->dt);

	p->dt=MIN(p->dt,maxtimestep);
	p->dt_old=p->dt;
    
    

}

double pftimestep::min(double val1,double val2,double val3)
{
	double mini;

	mini=val1;

	if(mini>val2)
	mini=val2;

	if(mini>val3)
	mini=val3;

	if(mini<0.0)
	mini=0.0;

	return mini;
}

double pftimestep::min(double val1,double val2)
{
	double mini;

	mini=val1;

	if(mini>val2)
	mini=val2;

	if(mini<0.0)
	mini=0.0;

	return mini;
}

double pftimestep::max(double val1,double val2,double val3)
{
	double maxi;

	maxi=val1;

	if(maxi<val2)
	maxi=val2;

	if(maxi<val3)
	maxi=val3;

	if(maxi<0.0)
	maxi=0.0;

	return maxi;
}

double pftimestep::max(double val1,double val2)
{
	double maxi;

	maxi=val1;

	if(maxi<val2)
	maxi=val2;

	if(maxi<0.0)
	maxi=0.0;

	return maxi;
}

