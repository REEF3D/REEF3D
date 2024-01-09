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

#include"poisson_nse.h"
#include"lexer.h"
#include"fdm.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_df.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"

poisson_nse::poisson_nse(lexer * p, heat *&pheat, concentration *&pconc)
{
    if((p->F80==0||p->A10==55) && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0 && p->X10==0)
	pd = new density_f(p);
    
    if((p->F80==0||p->A10==55) && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0 && p->X10==1)  
	pd = new density_df(p);
	
	if(p->F80==0 && p->H10==0 && p->W30==1)
	pd = new density_comp(p);
	
	if(p->F80==0 && p->H10>0)
	pd = new density_heat(p,pheat);
	
	if(p->F80==0 && p->C10>0)
	pd = new density_conc(p,pconc);
    
    if(p->F80>0 && p->H10==0 && p->W30==0)
	pd = new density_vof(p);
}

poisson_nse::~poisson_nse()
{
}

void poisson_nse::start(lexer* p, fdm *a, field &press)
{
	n=0;
    FLUIDLOOP
	{
        if(p->flag4[IJK]>0)
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
        }
        
        if(p->flag4[IJK]<0)
        {
        a->M.p[n] = 1.0;


        a->M.n[n] = 0.0;
        a->M.s[n] = 0.0;

        a->M.w[n] = 0.0;
        a->M.e[n] = 0.0;

        a->M.t[n] = 0.0;
        a->M.b[n] = 0.0;
        }
	
	++n;
	}
     
    n=0;
	FLUIDLOOP
	{
        if(p->flag4[IJK]>0)
        {
            // Solid boundaries
            if(p->flag4[Im1JK]<AIR)
            {
            a->rhsvec.V[n] -= a->M.s[n]*press(i-1,j,k);
            a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<AIR)
            {
            a->rhsvec.V[n] -= a->M.n[n]*press(i+1,j,k);
            a->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<AIR)
            {
            a->rhsvec.V[n] -= a->M.e[n]*press(i,j-1,k);
            a->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<AIR)
            {
            a->rhsvec.V[n] -= a->M.w[n]*press(i,j+1,k);
            a->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<AIR)
            {
            a->rhsvec.V[n] -= a->M.b[n]*press(i,j,k-1);
            a->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<AIR)
            {
            a->rhsvec.V[n] -= a->M.t[n]*press(i,j,k+1);
            a->M.t[n] = 0.0;
            }
            
            // FSFBC
            if(p->flag4[Im1JK]==AIR)
            {
                if(p->D37==1)
                {
                a->rhsvec.V[n] -= a->M.s[n]*press(i-1,j,k);
                a->M.s[n] = 0.0;
                }
                
                if(p->D37==2)
                {
                teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i-1,j,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DXN[IM1]/(fabs(a->phi(i,j,k-1))+fabs(a->phi(i,j,k)));
                
                //cout<<" Teta: "<<teta<<" a->phi(i,j,k): "<<a->phi(i,j,k)<<" a->phi(i-1,j,k): "<<a->phi(i-1,j,k)<<endl;
                
                a->M.p[n] -= (CPOR1m*PORVAL1m)/(pd->roface(p,a,-1,0,0)*p->DXP[IM1]*p->DXN[IP])*p->x_dir;
                a->M.p[n] += (CPOR1m*PORVAL1m)/(pd->roface(p,a,-1,0,0)*teta*p->DXP[IM1]*p->DXN[IP])*p->x_dir;
                           
                a->M.s[n] = 0.0;
                }
            }
            
            if(p->flag4[Ip1JK]==AIR)
            {
                if(p->D37==1)
                {
                a->rhsvec.V[n] -= a->M.n[n]*press(i+1,j,k);
                a->M.n[n] = 0.0;
                }
                
                if(p->D37==2)
                {
                teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i+1,j,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DXN[IP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k)));
                
                //cout<<" Teta: "<<teta<<" a->phi(i,j,k): "<<a->phi(i,j,k)<<" a->phi(i+1,j,k): "<<a->phi(i+1,j,k)<<endl;
                
                a->M.p[n] -= (CPOR1*PORVAL1)/(pd->roface(p,a,1,0,0)*p->DXP[IP]*p->DXN[IP])*p->x_dir;
                a->M.p[n] += (CPOR1*PORVAL1)/(pd->roface(p,a,1,0,0)*teta*p->DXP[IP]*p->DXN[IP])*p->x_dir;
                           
                a->M.n[n] = 0.0;
                }
            }
            
            if(p->flag4[IJm1K]==AIR)
            {
                if(p->D37==1)
                {
                a->rhsvec.V[n] -= a->M.e[n]*press(i,j-1,k);
                a->M.e[n] = 0.0;
                }
            }
            
            if(p->flag4[IJp1K]==AIR)
            {
                if(p->D37==1)
                {
                a->rhsvec.V[n] -= a->M.w[n]*press(i,j+1,k);
                a->M.w[n] = 0.0;
                }
            }
            
            if(p->flag4[IJKm1]==AIR)
            {
                if(p->D37==1)
                {
                a->rhsvec.V[n] -= a->M.b[n]*press(i,j,k-1);
                a->M.b[n] = 0.0;
                }
            }
            
            if(p->flag4[IJKp1]==AIR)
            {
                if(p->D37==1)
                {
                a->rhsvec.V[n] -= a->M.t[n]*press(i,j,k+1);
                a->M.t[n] = 0.0;
                }
            

                if(p->D37==2)
                {
                teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k)));
                
                //cout<<" Teta: "<<teta<<" a->phi(i,j,k): "<<a->phi(i,j,k)<<" a->phi(i+1,j,k): "<<a->phi(i+1,j,k)<<endl;
                
                a->M.p[n] -= (CPOR3*PORVAL3)/(pd->roface(p,a,0,0,1)*p->DZP[KP]*p->DZN[KP])*p->z_dir;
                a->M.p[n] += (CPOR3*PORVAL3)/(pd->roface(p,a,0,0,1)*teta*p->DZP[KP]*p->DZN[KP])*p->z_dir;
                           
                a->M.t[n] = 0.0;
                }
            }
            
        }

	++n;
	}
}

