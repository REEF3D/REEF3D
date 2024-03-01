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

#include"bc_ikomega.h"
#include"fdm.h"
#include"lexer.h"

bc_ikomega::bc_ikomega(lexer* p):roughness(p)
{
    kappa=0.4;
}

bc_ikomega::~bc_ikomega()
{
}

void bc_ikomega::bckomega_start(fdm* a,lexer* p,field& kin,field& eps,int gcval)
{
	int q;

	if(gcval==20)
	{
		QGC4LOOP
		if(p->gcb4[q][4]==5 || p->gcb4[q][4]==21 || p->gcb4[q][4]==22 || p->gcb4[q][4]==41 || p->gcb4[q][4]==42 || p->gcb4[q][4]==43)
		wall_law_kin(a,p,kin,eps,p->gcb4[q][0], p->gcb4[q][1], p->gcb4[q][2], p->gcb4[q][3], p->gcb4[q][4], p->gcb4[q][5],  p->gcd4[q]);
        
        QGCDF4LOOP
		wall_law_kin(a,p,kin,eps,p->gcdf4[q][0], p->gcdf4[q][1], p->gcdf4[q][2], p->gcdf4[q][3], p->gcdf4[q][4], p->gcdf4[q][5],  0.5*p->DXM);
        
        //wall_law_kin_df(a,p,kin,eps);
	}

// ----------------- 
	if(gcval==30)
	{
		QGC4LOOP
		if(p->gcb4[q][4]==5 || p->gcb4[q][4]==21 || p->gcb4[q][4]==22 || p->gcb4[q][4]==41 || p->gcb4[q][4]==42 || p->gcb4[q][4]==43  || (p->gcb4[q][4]==3 && p->gcb4[q][3]==6))
		wall_law_omega(a,p,kin,eps,p->gcb4[q][0], p->gcb4[q][1], p->gcb4[q][2], p->gcb4[q][3], p->gcb4[q][4], p->gcb4[q][5],  p->gcd4[q]);
        
        QGCDF4LOOP
		wall_law_omega(a,p,kin,eps,p->gcdf4[q][0], p->gcdf4[q][1], p->gcdf4[q][2], p->gcdf4[q][3], p->gcdf4[q][4], p->gcdf4[q][5],  0.5*p->DXM);
        
        //wall_law_omega_df(a,p,kin,eps);
	}
}

void bc_ikomega::wall_law_kin(fdm* a,lexer* p,field& kin,field& eps,int ii,int jj,int kk,int cs,int bc, int id, double dist)
{
    double uvel,vvel,wvel;
    double zval;
    
	i=ii;
	j=jj;
	k=kk;
    
    if(p->j_dir==0)
	dist = (1.0/4.0)*(p->DXN[IP] + p->DZN[KP]);
    
    if(p->j_dir==1)
    dist = (1.0/6.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
    
	ks=ks_val(p,a,ii,jj,kk,cs,bc);

        uvel=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
        vvel=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
        wvel=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
        
        if((bc==5 || (a->topo(i-1,j,k)<0.0 || a->topo(i+1,j,k-1)<0.0 || a->topo(i,j-1,k)<0.0 || a->topo(i,j+1,k)<0.0 || a->topo(i,j,k-1)<0.0)) && p->S10>0)
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
	a->rhsvec.V[id] += (tau*u_abs)/dist;
    
    
}

void bc_ikomega::wall_law_omega(fdm* a,lexer* p,field& kin,field& eps,int ii,int jj,int kk,int cs, int bc, int id, double dist)
{    
    i=ii;
	j=jj;
	k=kk;
	
    if(p->j_dir==0)
	dist=(1.0/4.0)*(p->DXN[IP] + p->DZN[KP]);
    
    if(p->j_dir==1)
    dist=(1.0/6.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);

	eps_star = pow((kin(i,j,k)>(0.0)?(kin(i,j,k)):(0.0)),0.5) / (0.4*dist*pow(p->cmu, 0.25));
    
    a->M.p[id] += 1.0e20;
	a->rhsvec.V[id] += eps_star*1.0e20;
    
    //eps(i,j,k) = eps_star;
}

void bc_ikomega::wall_law_kin_df(fdm* a,lexer* p,field& kin,field& eps)
{
    double uvel,vvel,wvel;
    double zval;
    double dirac, psi;
    
    n=0;
    LOOP
    {
    if (p->j_dir==0)
    psi = 1.1*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 
    
    if(p->j_dir==1)
    psi = 1.1*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);
    
    if(p->j_dir==0)
	dist = p->X41*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]);
    
    if(p->j_dir==1)
    dist = p->X41*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
    
    //dist = MAX(0.1*dist, fabs(MIN(a->solid(i,j,k),a->topo(i,j,k))));
    
	//ks=ks_val(p,a,i,j,k,cs,bc);
    ks=p->B50;

        uvel=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
        vvel=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
        wvel=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
        /*
        if((a->topo(i-1,j,k)<0.0 || a->topo(i+1,j,k-1)<0.0 || a->topo(i,j-1,k)<0.0 || a->topo(i,j+1,k)<0.0 || a->topo(i,j,k-1)<0.0) && p->S10>0 && p->S16==4)
        {
        zval = a->bed(i,j) + p->T43*p->DZN[KP];
        
            uvel=p->ccipol1(a->u,p->XP[IP],p->YP[JP],zval);
            vvel=p->ccipol2(a->v,p->XP[IP],p->YP[JP],zval);
            wvel=p->ccipol3(a->w,p->XP[IP],p->YP[JP],zval);   

        ks=p->S21*p->S20;
        }*/
        
        if((a->topo(i-1,j,k)<0.0 || a->topo(i+1,j,k-1)<0.0 || a->topo(i,j-1,k)<0.0 || a->topo(i,j+1,k)<0.0 || a->topo(i,j,k-1)<0.0) && p->S10>0 && p->S16==4)
        ks=p->S21*p->S20;
        
        u_abs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);

		if(30.0*dist<ks)
		dist=ks/30.0;
		
        uplus = (1.0/kappa)*MAX(0.01,log(30.0*(dist/ks)));

	tau=(u_abs*u_abs)/pow((uplus>0.0?uplus:(1.0e20)),2.0);
    
    dirac = 0.0;
    if(fabs(MIN(a->solid(i,j,k),a->topo(i,j,k)))<psi)
    dirac = (0.5/psi)*(1.0 + cos((PI*(MIN(a->solid(i,j,k),a->topo(i,j,k))))/psi));
        
    a->M.p[n]      += dirac*(pow(p->cmu,0.75)*pow(fabs(kin(i,j,k)),0.5)*uplus)/dist;
    a->rhsvec.V[n] += dirac*(tau*u_abs)/dist;
    a->test(i,j,k) = dirac;
    
    //kin(i,j,k) = dirac*tau*tau/pow(p->cmu,2.0);
    ++n;
    }
}

void bc_ikomega::wall_law_omega_df(fdm* a,lexer* p,field& kin,field& eps)
{
    double dirac, psi;
    
    n=0;
    LOOP
    {
    if (p->j_dir==0)
    psi = 1.1*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 
    
    if(p->j_dir==1)
    psi = 1.1*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);
    
    if(p->j_dir==0)
	dist = p->X41*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]);
    
    if(p->j_dir==1)
    dist = p->X41*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
    
    //dist = MAX(0.1*dist, fabs(MIN(a->solid(i,j,k),a->topo(i,j,k))));
    
    dirac = 0.0;
    if(fabs(MIN(a->solid(i,j,k),a->topo(i,j,k)))<psi)
    dirac = (0.5/psi)*(1.0 + cos((PI*(MIN(a->solid(i,j,k),a->topo(i,j,k))))/psi));
    
    
    eps_star = pow((kin(i,j,k)>(0.0)?(kin(i,j,k)):(0.0)),0.5)/(0.4*dist*pow(p->cmu, 0.25));
    
    a->M.p[n]      += dirac*1.0e20;
	a->rhsvec.V[n] += dirac*eps_star*1.0e20;
    
    //eps(i,j,k) = dirac*eps_star;
    ++n;
    }
}

void bc_ikomega::bckin_matrix(fdm* a,lexer* p,field& kin,field& eps)
{
        n=0;
        LOOP
        {
            if(p->flag4[Im1JK]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[Im1JK]<0))
            {
            if(p->BC[Im1JK]!=1)
            a->rhsvec.V[n] -= a->M.s[n]*kin(i,j,k);
            a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[Ip1JK]<0))
            {
            if(p->BC[Ip1JK]!=1)
            a->rhsvec.V[n] -= a->M.n[n]*kin(i,j,k);
            a->M.n[n] = 0.0;
            }
            
            if((p->flag4[IJm1K]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJm1K]<0)) && p->j_dir==1 && p->BC[IJm1K]==0)
            {
            a->rhsvec.V[n] -= a->M.e[n]*kin(i,j,k);
            a->M.e[n] = 0.0;
            }
            
            if((p->flag4[IJp1K]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJp1K]<0)) && p->j_dir==1 && p->BC[IJp1K]==0)
            {
            a->rhsvec.V[n] -= a->M.w[n]*kin(i,j,k);
            a->M.w[n] = 0.0;
            }
            
            if((p->flag4[IJKm1]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJKm1]<0)) && p->BC[IJKm1]==0)
            {
            a->rhsvec.V[n] -= a->M.b[n]*kin(i,j,k);
            a->M.b[n] = 0.0;
            }
            
            if((p->flag4[IJKp1]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJKp1]<0)) && p->BC[IJKp1]==0)
            {
            a->rhsvec.V[n] -= a->M.t[n]*kin(i,j,k);
            a->M.t[n] = 0.0;
            }
            
        ++n;
        }
        
        
    // turn off inside direct forcing body
    if(p->X10==1 || p->G3==1)
    {
        n=0;
        LOOP
        {
            if(p->flagsf4[IJK]<0)
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
    }
}

void bc_ikomega::bcomega_matrix(fdm* a,lexer* p,field& kin,field& eps)
{
    
    // bc
        n=0;
        LOOP
        {

            if(p->flag4[Im1JK]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[Im1JK]<0))
            {
            if(p->BC[Im1JK]!=1)
            a->rhsvec.V[n] -= a->M.s[n]*eps(i,j,k);
            a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[Ip1JK]<0))
            {
            if(p->BC[Ip1JK]!=1)
            a->rhsvec.V[n] -= a->M.n[n]*eps(i,j,k);
            a->M.n[n] = 0.0;
            }
            
            if((p->flag4[IJm1K]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJm1K]<0)) && p->j_dir==1 && p->BC[IJm1K]==0)
            {
            a->rhsvec.V[n] -= a->M.e[n]*eps(i,j,k);
            a->M.e[n] = 0.0;
            }
            
            if((p->flag4[IJp1K]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJp1K]<0)) && p->j_dir==1 && p->BC[IJp1K]==0)
            {
            a->rhsvec.V[n] -= a->M.w[n]*eps(i,j,k);
            a->M.w[n] = 0.0;
            }
            
            if((p->flag4[IJKm1]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJKm1]<0)) && p->BC[IJKm1]==0)
            {
            a->rhsvec.V[n] -= a->M.b[n]*eps(i,j,k);
            a->M.b[n] = 0.0;
            }
            
            if((p->flag4[IJKp1]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJKp1]<0)) && p->BC[IJKp1]==0)
            {
            a->rhsvec.V[n] -= a->M.t[n]*eps(i,j,k);
            a->M.t[n] = 0.0;
            }
            ++n;
            
        }
        
    // turn off inside direct forcing body
    if(p->X10==1 || p->G3==1)
    {
        n=0;
        LOOP
        {
            if(p->flagsf4[IJK]<0)
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
    }
}

