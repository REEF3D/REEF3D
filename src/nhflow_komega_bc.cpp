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

#include"nhflow_komega_bc.h"
#include"fdm_nhf.h"
#include"lexer.h"
 
nhflow_komega_bc::nhflow_komega_bc(lexer *p) : roughness(p)
{
    kappa=0.4;
}

nhflow_komega_bc::~nhflow_komega_bc()
{
}

void nhflow_komega_bc::bckomega_start(lexer *p, fdm_nhf *d, double *KIN, double *EPS, int gcval)
{
	if(gcval==20)
    wall_law_kin(p,d,KIN,EPS);
        
// ----------------- 

	if(gcval==30)
	wall_law_omega(p,d,KIN,EPS);
}

void nhflow_komega_bc::wall_law_kin(lexer *p, fdm_nhf *d, double *KIN, double *EPS)
{
    double uvel,vvel,wvel;
    double zval;
    int check;
    
    count=0;
    LOOP
    {
            check=0;
        
            if(k==0)
            check=1;
        
            if((p->flag4[Im1JK]<0 || p->DF[Im1JK]<0) && i+p->origin_i != 0)
            {
            dist = 0.5*p->DXN[IP];
            check=1;
            }

            if((p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0) && i+p->origin_i != p->gknox-1)
            {
            dist = 0.5*p->DXN[IP];
            check=1;
            }

            if((p->flag4[IJm1K]<0 || p->DF[IJm1K]<0) && p->j_dir==1)
            {
            dist = 0.5*p->DYN[JP];
            check=1;
            }
                
            if((p->flag4[IJp1K]<0 || p->DF[IJp1K]<0) && p->j_dir==1)
            {
            dist = 0.5*p->DYN[JP];
            check=1;
            }
                
            if(p->flag4[IJKm1]<0 || p->DF[IJKm1]<0)
            {
            dist = 0.5*p->DZN[KP]*d->WL(i,j);
            check=1;
            }

            if((p->flag4[IJKp1]<0 || p->DF[IJKp1]<0) && k!=p->knoz-1)
            {
            dist = 0.5*p->DZN[KP]*d->WL(i,j);
            check=1;
            }
        
        
        if(check==1)
        {
            ks=p->B50;
        

            uvel=d->U[IJK];
            vvel=d->V[IJK];
            wvel=d->W[IJK];

            
            if(k==0 && p->S10>0)
            ks = p->S20*p->S21;

            u_abs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);

            if(30.0*dist<ks)
            dist=ks/30.0;
            
            uplus = (1.0/kappa)*MAX(0.01,log(30.0*(dist/ks)));

            tau=(u_abs*u_abs)/pow((uplus>0.0?uplus:(1.0e20)),2.0);
        
        d->M.p[count] += (pow(p->cmu,0.75)*pow(fabs(KIN[IJK]),0.5)*uplus)/dist;
        d->rhsvec.V[count] += (tau*u_abs)/dist;
        }
    ++count;
    }
}

void nhflow_komega_bc::wall_law_omega(lexer *p, fdm_nhf *d, double *KIN, double *EPS)
{
    int check=0;
    
    
    count=0;
    LOOP
    {
        check=0;
        
            if(k==0)
            check=1;
        
            if((p->flag4[Im1JK]<0 || p->DF[Im1JK]<0) && i+p->origin_i != 0)
            {
            dist = 0.5*p->DXN[IP];
            check=1;
            }

            if((p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0) && i+p->origin_i != p->gknox-1)
            {
            dist = 0.5*p->DXN[IP];
            check=1;
            }

            if((p->flag4[IJm1K]<0 || p->DF[IJm1K]<0) && p->j_dir==1)
            {
            dist = 0.5*p->DYN[JP];
            check=1;
            }
                
            if((p->flag4[IJp1K]<0 || p->DF[IJp1K]<0) && p->j_dir==1)
            {
            dist = 0.5*p->DYN[JP];
            check=1;
            }
                
            if(p->flag4[IJKm1]<0 || p->DF[IJKm1]<0)
            {
            dist = 0.5*p->DZN[KP]*d->WL(i,j);
            check=1;
            }

            if((p->flag4[IJKp1]<0 || p->DF[IJKp1]<0) && k!=p->knoz-1)
            {
            dist = 0.5*p->DZN[KP]*d->WL(i,j);
            check=1;
            }
    
        if(check==1)
        {
        eps_star = pow((KIN[IJK]>(0.0)?(KIN[IJK]):(0.0)),0.5) / (0.4*dist*pow(p->cmu, 0.25));

        EPS[IJK] = eps_star;
        }
    }
}

void nhflow_komega_bc::bckin_matrix(lexer *p, fdm_nhf *d, double *KIN, double *EPS)
{
	int q;
    int inflow=0;
    int outflow=0;
    
    if(p->B98>=3 || p->B60==1)
    inflow=1;
    
    if(p->B99>=3)
    outflow=1;
    
    if(p->B60==1)
    outflow=2;
        
        n=0;
        LOOP
        {
            if((p->flag4[Im1JK]<0 || p->DF[Im1JK]<0))// && inflow==0)
            {
            d->rhsvec.V[n] -= d->M.s[n]*KIN[Im1JK];
            d->M.s[n] = 0.0;
            }
            
            if((p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0))// && outflow==0)
            {
            d->rhsvec.V[n] -= d->M.n[n]*KIN[Ip1JK];
            d->M.n[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJm1K]<0 || p->DF[IJm1K]<0)
            {
            d->rhsvec.V[n] -= d->M.e[n]*KIN[IJm1K];
            d->M.e[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJp1K]<0 || p->DF[IJp1K]<0)
            {
            d->rhsvec.V[n] -= d->M.w[n]*KIN[IJp1K];
            d->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0 || p->DF[IJKm1]<0)
            {
            d->rhsvec.V[n] -= d->M.b[n]*KIN[IJKm1];
            d->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0 || p->DF[IJKp1]<0)
            {
            d->rhsvec.V[n] -= d->M.t[n]*KIN[IJKp1];
            d->M.t[n] = 0.0;
            }

        ++n;
        }
    
    // turn off inside direct forcing body
        n=0;
        LOOP
        {
            if(p->DF[IJK]==0)
            {
            KIN[IJK] = 0.0;
            d->EV[IJK] = 0.0;
            
            d->M.p[n]  =   1.0;

            d->M.n[n] = 0.0;
            d->M.s[n] = 0.0;

            d->M.w[n] = 0.0;
            d->M.e[n] = 0.0;

            d->M.t[n] = 0.0;
            d->M.b[n] = 0.0;
            
            d->rhsvec.V[n] = 0.0;
            }
            ++n;
        }
}


void nhflow_komega_bc::bcomega_matrix(lexer *p, fdm_nhf *d, double *KIN, double *EPS)
{
    
	int q;
    int inflow=0;
    int outflow=0;
    
    if(p->B98>=3 || p->B60==1)
    inflow=1;
    
    if(p->B99>=3)
    outflow=1;
    
    if(p->B60==1)
    outflow=2;
    
        n=0;
        LOOP
        {
            if((p->flag4[Im1JK]<0 || p->DF[Im1JK]<0))// && inflow==0)
            {
            d->rhsvec.V[n] -= d->M.s[n]*EPS[Im1JK];
            d->M.s[n] = 0.0;
            }
            
            if((p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0))// && outflow==0)
            {
            d->rhsvec.V[n] -= d->M.n[n]*EPS[Ip1JK];
            d->M.n[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJm1K]<0 || p->DF[IJm1K]<0)
            {
            d->rhsvec.V[n] -= d->M.e[n]*EPS[IJm1K];
            d->M.e[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJp1K]<0 || p->DF[IJp1K]<0)
            {
            d->rhsvec.V[n] -= d->M.w[n]*EPS[IJp1K];
            d->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0 || p->DF[IJKm1]<0)
            {
            d->rhsvec.V[n] -= d->M.b[n]*EPS[IJKm1];
            d->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0 || p->DF[IJKp1]<0)
            {
            d->rhsvec.V[n] -= d->M.t[n]*EPS[IJKp1];
            d->M.t[n] = 0.0;
            }

        ++n;
        }
    
    
    // turn off inside direct forcing body
        n=0;
        LOOP
        {
            if(p->DF[IJK]==0)
            {
            EPS[IJK] = 0.0;
            
            
            
            d->M.p[n]  =   1.0;

            d->M.n[n] = 0.0;
            d->M.s[n] = 0.0;

            d->M.w[n] = 0.0;
            d->M.e[n] = 0.0;

            d->M.t[n] = 0.0;
            d->M.b[n] = 0.0;
            
            d->rhsvec.V[n] = 0.0;
            }
            ++n;
        }
}
