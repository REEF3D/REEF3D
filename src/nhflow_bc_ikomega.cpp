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

#include"nhflow_bc_ikomega.h"
#include"fdm_nhf.h"
#include"lexer.h"
 
nhflow_bc_ikomega::nhflow_bc_ikomega(lexer *p) : roughness(p)
{
    kappa=0.4;
}

nhflow_bc_ikomega::~nhflow_bc_ikomega()
{
}

void nhflow_bc_ikomega::bckomega_start(lexer *p, fdm_nhf *d, double *KIN, double *EPS, int gcval)
{
    /*
	int q;

	if(gcval==20)
	{
		QGC4LOOP
		if(p->gcb4[q][4]==5 || p->gcb4[q][4]==21 || p->gcb4[q][4]==22 || p->gcb4[q][4]==41 || p->gcb4[q][4]==42 || p->gcb4[q][4]==43)
		wall_law_kin(a,p,kin,eps,p->gcb4[q][0], p->gcb4[q][1], p->gcb4[q][2], p->gcb4[q][3], p->gcb4[q][4], p->gcb4[q][5],  p->gcd4[q]);
        
        QGCDF4LOOP
		wall_law_kin(a,p,kin,eps,p->gcdf4[q][0], p->gcdf4[q][1], p->gcdf4[q][2], p->gcdf4[q][3], p->gcdf4[q][4], p->gcdf4[q][5],  0.5*p->DXM);
	}

// ----------------- 
	if(gcval==30)
	{
		QGC4LOOP
		if(p->gcb4[q][4]==5 || p->gcb4[q][4]==21 || p->gcb4[q][4]==22 || p->gcb4[q][4]==41 || p->gcb4[q][4]==42 || p->gcb4[q][4]==43) 
		wall_law_omega(a,p,kin,eps,p->gcb4[q][0], p->gcb4[q][1], p->gcb4[q][2], p->gcb4[q][3], p->gcb4[q][4], p->gcb4[q][5],  p->gcd4[q]);
        
        QGCDF4LOOP
		wall_law_omega(a,p,kin,eps,p->gcdf4[q][0], p->gcdf4[q][1], p->gcdf4[q][2], p->gcdf4[q][3], p->gcdf4[q][4], p->gcdf4[q][5],  0.5*p->DXM);
	}*/
}

void nhflow_bc_ikomega::wall_law_kin(lexer *p, fdm_nhf *d, double *KIN, double *EPS,
                                     int ii,int jj,int kk,int cs,int bc, int id, double dist)
{
    /*
    double uvel,vvel,wvel;
    double zval;
	dist=0.5*p->DXM;
    
    if(p->j_dir==0)        
    dist = 0.5*(1.0/2.0)*(p->DXN[IP]+p->DZN[KP]);
        
    if(p->j_dir==1)
    dist = 0.5*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

	i=ii;
	j=jj;
	k=kk;
	
	ks=ks_val(p,a,ii,jj,kk,cs,bc);


        uvel=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
        vvel=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
        wvel=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
        
        if((bc==5 || (a->topo(i-1,j,k)<0.0 || a->topo(i+1,j,k-1)<0.0 || a->topo(i,j-1,k)<0.0 || a->topo(i,j+1,k)<0.0 || a->topo(i,j,k-1)<0.0)) && p->S10>0 && p->S16==4)
        {
        zval = a->bed(i,j) + p->T43*p->DZN[KP];
        
            uvel=p->ccipol1(a->u,p->XP[IP],p->YP[JP],zval);
            vvel=p->ccipol2(a->v,p->XP[IP],p->YP[JP],zval);
            wvel=p->ccipol3(a->w,p->XP[IP],p->YP[JP],zval);    
        }
        
        u_abs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);

		if(30.0*dist<ks)
		dist=ks/30.0;
		
        uplus = (1.0/kappa)*MAX(0.01,log(30.0*(dist/ks)));

	tau=(u_abs*u_abs)/pow((uplus>0.0?uplus:(1.0e20)),2.0);
	
	a->M.p[id] += (pow(p->cmu,0.75)*pow(fabs(kin(i,j,k)),0.5)*uplus)/dist;
	a->rhsvec.V[id] += (tau*u_abs)/dist;*/
}

void nhflow_bc_ikomega::wall_law_omega(lexer *p, fdm_nhf *d, double *KIN, double *EPS,
                                        int ii,int jj,int kk,int cs, int bc, int id, double dist)
{
    /*
	i=ii;
	j=jj;
	k=kk;
    
	dist=0.5*p->DXM;
    
    if(p->j_dir==0)        
    dist = 0.5*(1.0/2.0)*(p->DXN[IP]+p->DZN[KP]);
        
    if(p->j_dir==1)
    dist = 0.5*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

	eps_star = pow((kin(i,j,k)>(0.0)?(kin(i,j,k)):(0.0)),0.5) / (0.4*dist*pow(p->cmu, 0.25));

	eps(i,j,k) = eps_star;*/
}

void nhflow_bc_ikomega::bckin_matrix(lexer *p, fdm_nhf *d, double *KIN, double *EPS)
{
    /*
	int q;
    
    // set to zero inside direct forcing body
    if(p->X10==1)
    LOOP
    if(a->fb(i,j,k)<0.0)
    kin(i,j,k)=0.0;
    
    // set to zero inside solid forcing body
    if(p->G3==1)
    LOOP
    if(a->solid(i,j,k)<0.0 || a->topo(i,j,k)<0.0)
    kin(i,j,k)=0.0;

        n=0;
        LOOP
        {
            if(p->flag4[Im1JK]<0 || (p->X10==1 && a->fb(i-1,j,k)<0.0)
            || (p->G3==1 && (a->solid(i-1,j,k)<0.0 || a->topo(i-1,j,k)<0.0)))
            {
            a->rhsvec.V[n] -= a->M.s[n]*kin(i-1,j,k);
            a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0 || (p->X10==1 && a->fb(i+1,j,k)<0.0)
            || (p->G3==1 && (a->solid(i+1,j,k)<0.0 || a->topo(i+1,j,k)<0.0)))
            {
            a->rhsvec.V[n] -= a->M.n[n]*kin(i+1,j,k);
            a->M.n[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJm1K]<0 || (p->X10==1 && a->fb(i,j-1,k)<0.0)
            || (p->G3==1 && (a->solid(i,j-1,k)<0.0 || a->topo(i,j-1,k)<0.0)))
            {
            a->rhsvec.V[n] -= a->M.e[n]*kin(i,j-1,k);
            a->M.e[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJp1K]<0 || (p->X10==1 && a->fb(i,j+1,k)<0.0)
            || (p->G3==1 && (a->solid(i,j+1,k)<0.0 || a->topo(i,j+1,k)<0.0)))
            {
            a->rhsvec.V[n] -= a->M.w[n]*kin(i,j+1,k);
            a->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0 || (p->X10==1 && a->fb(i,j,k-1)<0.0)
            || (p->G3==1 && (a->solid(i,j,k-1)<0.0 || a->topo(i,j,k-1)<0.0)))
            {
            a->rhsvec.V[n] -= a->M.b[n]*kin(i,j,k-1);
            a->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0  || (p->X10==1 && a->fb(i,j,k+1)<0.0)
            || (p->G3==1 && (a->solid(i,j,k+1)<0.0 || a->topo(i,j,k+1)<0.0)))
            {
            a->rhsvec.V[n] -= a->M.t[n]*kin(i,j,k+1);
            a->M.t[n] = 0.0;
            }

        ++n;
        }
    
    // turn off inside direct forcing body
    if((p->X10==1) || (p->G3==1 && (a->solid(i,j,k)<0.0 || a->topo(i,j,k)<0.0)))
    {
    
        n=0;
        LOOP
        {
            if(a->fb(i,j,k)<0.0 || a->solid(i,j,k)<0.0 || a->topo(i,j,k)<0.0)
            {
            a->M.p[n]  =   1.0;

            a->M.n[n] = 0.0;
            a->M.s[n] = 0.0;

            a->M.w[n] = 0.0;
            a->M.e[n] = 0.0;

            a->M.t[n] = 0.0;
            a->M.b[n] = 0.0;
            
            a->rhsvec.V[n] = 0.0;
            }
            ++n;
        }
    }*/
}


void nhflow_bc_ikomega::bcomega_matrix(lexer *p, fdm_nhf *d, double *KIN, double *EPS)
{
    /*
	int q;
    
    // set to zero inside direct forcing body
    if(p->X10==1)
    LOOP
    eps(i,j,k)=0.0;
    
    // set to zero inside solid forcing body
    if(p->G3==1)
    LOOP
    if(a->solid(i,j,k)<0.0 || a->topo(i,j,k)<0.0)
    eps(i,j,k)=0.0;
        
        n=0;
        LOOP
        {
            if(p->flag4[Im1JK]<0 || (p->X10==1 && a->fb(i-1,j,k)<0.0)
            || (p->G3==1 && (a->solid(i-1,j,k)<0.0 || a->topo(i-1,j,k)<0.0)))
            {
            a->rhsvec.V[n] -= a->M.s[n]*eps(i-1,j,k);
            a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0 || (p->X10==1 && a->fb(i+1,j,k)<0.0)
            || (p->G3==1 && (a->solid(i+1,j,k)<0.0 || a->topo(i+1,j,k)<0.0)))
            {
            a->rhsvec.V[n] -= a->M.n[n]*eps(i+1,j,k);
            a->M.n[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJm1K]<0 || (p->X10==1 && a->fb(i,j-1,k)<0.0)
            || (p->G3==1 && (a->solid(i,j-1,k)<0.0 || a->topo(i,j-1,k)<0.0)))
            {
            a->rhsvec.V[n] -= a->M.e[n]*eps(i,j-1,k);
            a->M.e[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJp1K]<0 || (p->X10==1 && a->fb(i,j+1,k)<0.0)
            || (p->G3==1 && (a->solid(i,j+1,k)<0.0 || a->topo(i,j+1,k)<0.0)))
            {
            a->rhsvec.V[n] -= a->M.w[n]*eps(i,j+1,k);
            a->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0 || (p->X10==1 && a->fb(i,j,k-1)<0.0)
            || (p->G3==1 && (a->solid(i,j,k-1)<0.0 || a->topo(i,j,k-1)<0.0)))
            {
            a->rhsvec.V[n] -= a->M.b[n]*eps(i,j,k-1);
            a->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0  || (p->X10==1 && a->fb(i,j,k+1)<0.0)
            || (p->G3==1 && (a->solid(i,j,k+1)<0.0 || a->topo(i,j,k+1)<0.0)))
            {
            a->rhsvec.V[n] -= a->M.t[n]*eps(i,j,k+1);
            a->M.t[n] = 0.0;
            }

        ++n;
        }
    
    
    // turn off inside direct forcing body
    if((p->X10==1) || (p->G3==1 && (a->solid(i,j,k)<0.0 || a->topo(i,j,k)<0.0)))
    {
        n=0;
        LOOP
        {
            if(a->fb(i,j,k)<0.0 || a->solid(i,j,k)<0.0 || a->topo(i,j,k)<0.0)
            {
            a->M.p[n]  =   1.0;

            a->M.n[n] = 0.0;
            a->M.s[n] = 0.0;

            a->M.w[n] = 0.0;
            a->M.e[n] = 0.0;

            a->M.t[n] = 0.0;
            a->M.b[n] = 0.0;
            
            a->rhsvec.V[n] = 0.0;
            }
            ++n;
        }
    }*/
}
