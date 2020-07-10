/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"poisson_f.h"
#include"lexer.h"
#include"fdm.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"

poisson_f::poisson_f(lexer *p, heat *&pheat, concentration *&pconc) 
{
    if((p->F80==0||p->A10==4) && p->H10==0 && p->W30==0)
	pd = new density_f(p);
	
	if(p->F80==0 && p->H10==0 && p->W30==1)
	pd = new density_comp(p);
	
	if(p->F80==0 && p->H10>0)
	pd = new density_heat(p,pheat);
	
	if(p->F80==0 && p->C10>0)
	pd = new density_conc(p,pconc);
    
    if(p->F80>0 && p->H10==0 && p->W30==0)
	pd = new density_vof(p);

}

poisson_f::~poisson_f()
{
}
/*
void poisson_f::start(lexer* p, fdm *a, field &press)
{	
   
	n=0;
    LOOP
	{
	a->M.p[n]  =   (CPOR1*PORVAL1)/(pd->roface(p,a,1,0,0)*p->DRM*p->DRM)*p->x_dir*p->DRDXP[IP]*p->DRDXN[IP]
                + (CPOR1m*PORVAL1m)/(pd->roface(p,a,-1,0,0)*p->DRM*p->DRM)*p->x_dir*p->DRDXP[IM1]*p->DRDXN[IP]
                
                + (CPOR2*PORVAL2)/(pd->roface(p,a,0,1,0)*p->DSM*p->DSM)*p->y_dir*p->DSDYP[JP]*p->DSDYN[JP]
                + (CPOR2m*PORVAL2m)/(pd->roface(p,a,0,-1,0)*p->DSM*p->DSM)*p->y_dir*p->DSDYP[JM1]*p->DSDYN[JP]
                
                + (CPOR3*PORVAL3)/(pd->roface(p,a,0,0,1)*p->DTM*p->DTM)*p->z_dir*p->DTDZP[KP]*p->DTDZN[KP]
                + (CPOR3m*PORVAL3m)/(pd->roface(p,a,0,0,-1)*p->DTM*p->DTM)*p->z_dir*p->DTDZP[KM1]*p->DTDZN[KP];


   	a->M.n[n] = -(CPOR1*PORVAL1)/(pd->roface(p,a,1,0,0)*p->DRM*p->DRM)*p->x_dir*p->DRDXP[IP]*p->DRDXN[IP]
               -(CPOR1*PORVAL1)/(pd->roface(p,a,1,0,0)*2.0*p->DRM)*p->x_dir*p->DDRDDXP[IP];
                
	a->M.s[n] = -(CPOR1m*PORVAL1m)/(pd->roface(p,a,-1,0,0)*p->DRM*p->DRM)*p->x_dir*p->DRDXP[IM1]*p->DRDXN[IP]
               +(CPOR1m*PORVAL1m)/(pd->roface(p,a,-1,0,0)*2.0*p->DRM)*p->x_dir*p->DDRDDXP[IP];
               

	a->M.w[n] = -(CPOR2*PORVAL2)/(pd->roface(p,a,0,1,0)*p->DSM*p->DSM)*p->y_dir*p->DSDYP[JP]*p->DSDYN[JP]
               +(CPOR2*PORVAL2)/(pd->roface(p,a,0,1,0)*2.0*p->DSM)*p->y_dir*p->DDSDDYN[JP];
    
	a->M.e[n] = -(CPOR2m*PORVAL2m)/(pd->roface(p,a,0,-1,0)*p->DSM*p->DSM)*p->y_dir*p->DSDYP[JM1]*p->DSDYN[JP]
               -(CPOR2m*PORVAL2m)/(pd->roface(p,a,0,-1,0)*2.0*p->DSM)*p->y_dir*p->DDSDDYN[JP];
    

	a->M.t[n] = -(CPOR3*PORVAL3)/(pd->roface(p,a,0,0,1)*p->DTM*p->DTM)*p->z_dir*p->DTDZP[KP]*p->DTDZN[KP]
                -(CPOR3*PORVAL3)/(pd->roface(p,a,0,0,1)*2.0*p->DTM)*p->z_dir*p->DDTDDZN[KP];
    
	a->M.b[n] = -(CPOR3m*PORVAL3m)/(pd->roface(p,a,0,0,-1)*p->DTM*p->DTM)*p->z_dir*p->DTDZP[KM1]*p->DTDZN[KP]
               +(CPOR3m*PORVAL3m)/(pd->roface(p,a,0,0,-1)*2.0*p->DTM)*p->z_dir*p->DDTDDZN[KP];
	
	++n;
	}
    
    n=0;
	LOOP
	{
        
		if(p->flag4[Im1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.s[n]*press(i-1,j,k);
		a->M.s[n] = 0.0;
		}
		
		if(p->flag4[Ip1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.n[n]*press(i+1,j,k);
		a->M.n[n] = 0.0;
		}
		
		if(p->flag4[IJm1K]<0)
		{
		a->rhsvec.V[n] -= a->M.e[n]*press(i,j-1,k);
		a->M.e[n] = 0.0;
		}
		
		if(p->flag4[IJp1K]<0)
		{
		a->rhsvec.V[n] -= a->M.w[n]*press(i,j+1,k);
		a->M.w[n] = 0.0;
		}
		
		if(p->flag4[IJKm1]<0)
		{
		a->rhsvec.V[n] -= a->M.b[n]*press(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag4[IJKp1]<0)
		{
		a->rhsvec.V[n] -= a->M.t[n]*press(i,j,k+1);
		a->M.t[n] = 0.0;
		}
	++n;
	}
}
*/



void poisson_f::start(lexer* p, fdm *a, field &press)
{	
   
	n=0;
    LOOP
	{
	a->M.p[n]  =   (CPOR1*PORVAL1)/(pd->roface(p,a,1,0,0)*p->DXP[IP]*p->DXN[IP])*p->x_dir
                + (CPOR1m*PORVAL1m)/(pd->roface(p,a,-1,0,0)*p->DXP[IM1]*p->DXN[IP])*p->x_dir
                
                + (CPOR2*PORVAL2)/(pd->roface(p,a,0,1,0)*p->DYP[JP]*p->DYN[JP])*p->y_dir
                + (CPOR2m*PORVAL2m)/(pd->roface(p,a,0,-1,0)*p->DYP[JM1]*p->DYN[JP])*p->y_dir
                
                + (CPOR3*PORVAL3)/(pd->roface(p,a,0,0,1)*p->DZP[KP]*p->DZN[KP])*p->z_dir
                + (CPOR3m*PORVAL3m)/(pd->roface(p,a,0,0,-1)*p->DZP[KM1]*p->DZN[KP])*p->z_dir;


   	a->M.n[n] = -(CPOR1*PORVAL1)/(pd->roface(p,a,1,0,0)*p->DXP[IP]*p->DXN[IP])*p->x_dir;
	a->M.s[n] = -(CPOR1m*PORVAL1m)/(pd->roface(p,a,-1,0,0)*p->DXP[IM1]*p->DXN[IP])*p->x_dir;

	a->M.w[n] = -(CPOR2*PORVAL2)/(pd->roface(p,a,0,1,0)*p->DYP[JP]*p->DYN[JP])*p->y_dir;
	a->M.e[n] = -(CPOR2m*PORVAL2m)/(pd->roface(p,a,0,-1,0)*p->DYP[JM1]*p->DYN[JP])*p->y_dir;

	a->M.t[n] = -(CPOR3*PORVAL3)/(pd->roface(p,a,0,0,1)*p->DZP[KP]*p->DZN[KP])*p->z_dir;
	a->M.b[n] = -(CPOR3m*PORVAL3m)/(pd->roface(p,a,0,0,-1)*p->DZP[KM1]*p->DZN[KP])*p->z_dir;
	
	++n;
	}
    
    n=0;
	LOOP
	{
        
		if(p->flag4[Im1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.s[n]*press(i-1,j,k);
		a->M.s[n] = 0.0;
		}
		
		if(p->flag4[Ip1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.n[n]*press(i+1,j,k);
		a->M.n[n] = 0.0;
		}
		
		if(p->flag4[IJm1K]<0)
		{
		a->rhsvec.V[n] -= a->M.e[n]*press(i,j-1,k);
		a->M.e[n] = 0.0;
		}
		
		if(p->flag4[IJp1K]<0)
		{
		a->rhsvec.V[n] -= a->M.w[n]*press(i,j+1,k);
		a->M.w[n] = 0.0;
		}
		
		if(p->flag4[IJKm1]<0)
		{
		a->rhsvec.V[n] -= a->M.b[n]*press(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag4[IJKp1]<0)
		{
		a->rhsvec.V[n] -= a->M.t[n]*press(i,j,k+1);
		a->M.t[n] = 0.0;
		}
	++n;
	}
}


