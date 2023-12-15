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

#include"poisson_pcorr.h"
#include"lexer.h"
#include"fdm.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_df.h"
#include"density_sf.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"
#include"density_rheo.h"

poisson_pcorr::poisson_pcorr(lexer *p, heat *&pheat, concentration *&pconc) 
{
    if((p->F80==0||p->A10==55) && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0 && p->X10==0)
	pd = new density_f(p);
    
    if((p->F80==0||p->A10==55) && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0 && p->X10==1)  
	pd = new density_df(p);
    
	if(p->F80==0 && p->H10==0 && p->W30==1  && p->F300==0 && p->W90==0)
	pd = new density_comp(p);
	
	if(p->F80==0 && p->H10>0 && p->F300==0 && p->W90==0)
	pd = new density_heat(p,pheat);
	
	if(p->F80==0 && p->C10>0 && p->F300==0 && p->W90==0)
	pd = new density_conc(p,pconc);
    
    if(p->F80>0 && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0)
	pd = new density_vof(p);
    
    if(p->F30>0 && p->H10==0 && p->W30==0  && p->F300==0 && p->W90>0)
    pd = new density_rheo(p);
    
    if(p->F300>=1)
    pd = new density_rheo(p);
    
    if(p->G3==1)  
	pd = new density_sf(p);
}

poisson_pcorr::~poisson_pcorr()
{
}

void poisson_pcorr::start(lexer* p, fdm *a, field &press)
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
        // inflow
		if(p->flag4[Im1JK]<0 && (i+p->origin_i>0 || p->periodic1==0))
		{
		a->rhsvec.V[n] -= a->M.s[n]*press(i-1,j,k);
		a->M.s[n] = 0.0;
		}
        
        /*
        if(p->flag4[Im1JK]<0 && (i+p->origin_i>0 || p->periodic1==0) && p->BC[Im1JK]==1)
		{
             pval=a->press(i,j,k);
             
		//a->rhsvec.V[n] -= a->M.s[n]*(-a->press(i,j,k)+pval);
        a->rhsvec.V[n] -= a->M.s[n]*press(i-1,j,k);
		a->M.s[n] = 0.0;
		}*/
		
        // outflow
		if(p->flag4[Ip1JK]<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0) && (p->BC[Ip1JK]!=2 || p->B60!=1))
		{
		a->rhsvec.V[n] -= a->M.n[n]*press(i+1,j,k);
		a->M.n[n] = 0.0;
		}
        
         if(p->flag4[Ip1JK]<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0) && (p->BC[Ip1JK]==2 && p->B60==1))
		{
             if(p->B77==1)
             {
             if(p->F50==2 || p->F50==3)
             pval=(p->fsfout - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
             
             if(p->F50==1 || p->F50==4)
             pval=a->press(i,j,k);
             }
             
             if(p->B77==2)
             pval=0.0;
        
		a->rhsvec.V[n] -= a->M.n[n]*(-a->press(i,j,k)+pval);
		a->M.n[n] = 0.0;
		}
        
    
		if(p->flag4[IJm1K]<0 && (j+p->origin_j>0 || p->periodic2==0))
		{
		a->rhsvec.V[n] -= a->M.e[n]*press(i,j-1,k);
		a->M.e[n] = 0.0;
		}
		
		if(p->flag4[IJp1K]<0 && (j+p->origin_j<p->gknoy-1 || p->periodic2==0))
		{
		a->rhsvec.V[n] -= a->M.w[n]*press(i,j+1,k);
		a->M.w[n] = 0.0;
		}
		
		if(p->flag4[IJKm1]<0 && (k+p->origin_k>0 || p->periodic3==0))
		{
		a->rhsvec.V[n] -= a->M.b[n]*press(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag4[IJKp1]<0 && (k+p->origin_k<p->gknoz-1 || p->periodic3==0))
		{
		a->rhsvec.V[n] -= a->M.t[n]*press(i,j,k+1);
		a->M.t[n] = 0.0;
		}
	++n;
	}
    
  
}
