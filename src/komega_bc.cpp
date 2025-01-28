/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"komega_bc.h"
#include"fdm.h"
#include"lexer.h"

komega_bc::komega_bc(lexer* p):roughness(p)
{
    kappa=0.4;
}

komega_bc::~komega_bc()
{
}

void komega_bc::bckomega_start(fdm* a,lexer* p,field& kin,field& eps,int gcval)
{
	int q;

	if(gcval==20)
	{
		QGC4LOOP
		if(p->gcb4[q][4]==5 || p->gcb4[q][4]==21 || p->gcb4[q][4]==22 || p->gcb4[q][4]==41 || p->gcb4[q][4]==42 || p->gcb4[q][4]==43)
		wall_law_kin(a,p,kin,eps,p->gcb4[q][0], p->gcb4[q][1], p->gcb4[q][2], p->gcb4[q][3], p->gcb4[q][4], p->gcb4[q][5],  p->gcd4[q]);
        
        QGCDF4LOOP
		wall_law_kin(a,p,kin,eps,p->gcdf4[q][0], p->gcdf4[q][1], p->gcdf4[q][2], p->gcdf4[q][3], p->gcdf4[q][4], p->gcdf4[q][5],  0.5*p->DXM);
        
        if(p->S10==2 || p->B269==1 || p->B269==2)
        vrans_wall_law_kin(p,a,kin,eps);
	}

// ----------------- 
	if(gcval==30)
	{
		QGC4LOOP
		if(p->gcb4[q][4]==5 || p->gcb4[q][4]==21 || p->gcb4[q][4]==22 || p->gcb4[q][4]==41 || p->gcb4[q][4]==42 || p->gcb4[q][4]==43  || (p->gcb4[q][4]==3 && p->gcb4[q][3]==6))
		wall_law_omega(a,p,kin,eps,p->gcb4[q][0], p->gcb4[q][1], p->gcb4[q][2], p->gcb4[q][3], p->gcb4[q][4], p->gcb4[q][5],  p->gcd4[q]);
        
        QGCDF4LOOP
		wall_law_omega(a,p,kin,eps,p->gcdf4[q][0], p->gcdf4[q][1], p->gcdf4[q][2], p->gcdf4[q][3], p->gcdf4[q][4], p->gcdf4[q][5],  0.5*p->DXM);
        
        if(p->S10==2 || p->B269==1 || p->B269==2)
        vrans_wall_law_omega(p,a,kin,eps);
	}
}

void komega_bc::wall_law_kin(fdm* a,lexer* p,field& kin,field& eps,int ii,int jj,int kk,int cs,int bc, int id, double dist)
{
    double uvel,vvel,wvel;
    double zval;
    
	i=ii;
	j=jj;
	k=kk;
    
    if(cs==1 || cs==4)
    dist = 0.5*p->DXN[IP];
    
    if(cs==2 || cs==3)
    dist = 0.5*p->DYN[JP];
    
    if(cs==5 || cs==6)
    dist = 0.5*p->DZN[KP];
    
	ks=ks_val(p,a,ii,jj,kk,cs,bc);

        uvel=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
        vvel=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
        wvel=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
        
        if((bc==5 || (a->topo(i-1,j,k)<0.0 || a->topo(i+1,j,k)<0.0 || a->topo(i,j-1,k)<0.0 || a->topo(i,j+1,k)<0.0 || a->topo(i,j,k-1)<0.0)) && p->S10>0)
        {
        zval = a->bed(i,j) + p->T44*p->DZN[KP];
        
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

void komega_bc::wall_law_omega(fdm* a,lexer* p,field& kin,field& eps,int ii,int jj,int kk,int cs, int bc, int id, double dist)
{    
    i=ii;
	j=jj;
	k=kk;
	
    if(cs==1 || cs==4)
    dist = 0.5*p->DXN[IP];
    
    if(cs==2 || cs==3)
    dist = 0.5*p->DYN[JP];
    
    if(cs==5 || cs==6)
    dist = 0.5*p->DZN[KP];

	eps_star = pow((kin(i,j,k)>(0.0)?(kin(i,j,k)):(0.0)),0.5) / (0.4*dist*pow(p->cmu, 0.25));
    
    a->M.p[id] += 1.0e20;
    a->rhsvec.V[id] += eps_star*1.0e20;
}

void komega_bc::bckin_matrix(fdm* a,lexer* p,field& kin,field& eps)
{
        n=0;
        LOOP
        {
            if(p->flag4[Im1JK]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[Im1JK]<0))
            {
            if(p->IO[Im1JK]!=1)
            a->rhsvec.V[n] -= a->M.s[n]*kin(i,j,k);
            a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[Ip1JK]<0))
            {
            if(p->IO[Ip1JK]!=1)
            a->rhsvec.V[n] -= a->M.n[n]*kin(i,j,k);
            a->M.n[n] = 0.0;
            }
            
            if((p->flag4[IJm1K]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJm1K]<0)) && p->j_dir==1 && p->IO[IJm1K]==0)
            {
            a->rhsvec.V[n] -= a->M.e[n]*kin(i,j,k);
            a->M.e[n] = 0.0;
            }
            
            if((p->flag4[IJp1K]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJp1K]<0)) && p->j_dir==1 && p->IO[IJp1K]==0)
            {
            a->rhsvec.V[n] -= a->M.w[n]*kin(i,j,k);
            a->M.w[n] = 0.0;
            }
            
            if((p->flag4[IJKm1]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJKm1]<0)) && p->IO[IJKm1]==0)
            {
            a->rhsvec.V[n] -= a->M.b[n]*kin(i,j,k);
            a->M.b[n] = 0.0;
            }
            
            if((p->flag4[IJKp1]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJKp1]<0)) && p->IO[IJKp1]==0)
            {
            a->rhsvec.V[n] -= a->M.t[n]*kin(i,j,k);
            a->M.t[n] = 0.0;
            }
            
        ++n;
        }
        
        
    // turn off inside direct forcing body
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

void komega_bc::bcomega_matrix(fdm* a,lexer* p,field& kin,field& eps)
{
    
    // bc
        n=0;
        LOOP
        {

            if(p->flag4[Im1JK]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[Im1JK]<0))
            {
            if(p->IO[Im1JK]!=1)
            a->rhsvec.V[n] -= a->M.s[n]*eps(i,j,k);
            a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[Ip1JK]<0))
            {
            if(p->IO[Ip1JK]!=1)
            a->rhsvec.V[n] -= a->M.n[n]*eps(i,j,k);
            a->M.n[n] = 0.0;
            }
            
            if((p->flag4[IJm1K]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJm1K]<0)) && p->j_dir==1 && p->IO[IJm1K]==0)
            {
            a->rhsvec.V[n] -= a->M.e[n]*eps(i,j,k);
            a->M.e[n] = 0.0;
            }
            
            if((p->flag4[IJp1K]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJp1K]<0)) && p->j_dir==1 && p->IO[IJp1K]==0)
            {
            a->rhsvec.V[n] -= a->M.w[n]*eps(i,j,k);
            a->M.w[n] = 0.0;
            }
            
            if((p->flag4[IJKm1]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJKm1]<0)) && p->IO[IJKm1]==0)
            {
            a->rhsvec.V[n] -= a->M.b[n]*eps(i,j,k);
            a->M.b[n] = 0.0;
            }
            
            if((p->flag4[IJKp1]<0 || (p->flagsf4[IJK]>0 && p->flagsf4[IJKp1]<0)) && p->IO[IJKp1]==0)
            {
            a->rhsvec.V[n] -= a->M.t[n]*eps(i,j,k);
            a->M.t[n] = 0.0;
            }
            ++n;
            
        }
        
    // turn off inside direct forcing body
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

void komega_bc::vrans_wall_law_kin(lexer *p,fdm *a,field &kin,field &eps)
{
    double uvel,vvel,wvel;
    double fac = 3.0;
    
    n=0;
    LOOP
    {
        if(a->porosity(i,j,k)>=0.99 &&  (a->porosity(i-1,j,k)<0.99 || a->porosity(i+1,j,k)<0.99 
                                        || a->porosity(i,j-1,k)<0.99 || a->porosity(i,j+1,k)<0.99)
                                        || a->porosity(i,j,k-1)<0.99 || a->porosity(i,j,k+1)<0.99)
        {
        if(p->j_dir==0)
        dist = (1.0/4.0)*(p->DXN[IP] + p->DZN[KP]);
        
        if(p->j_dir==1)
        dist = (1.0/6.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
        
        if(a->porosity(i-1,j,k)<0.99)
        ks=fac*a->porpart(i-1,j,k);
        
        if(a->porosity(i+1,j,k)<0.99)
        ks=fac*a->porpart(i+1,j,k);
        
        if(fac*a->porosity(i,j-1,k)<0.99)
        ks=a->porpart(i,j-1,k);
        
        if(a->porosity(i,j+1,k)<0.99)
        ks=fac*a->porpart(i,j+1,k);
        
        if(a->porosity(i,j,k-1)<0.99)
        ks=fac*a->porpart(i,j,k);
        
        if(a->porosity(i,j,k+1)<0.99)
        ks=fac*a->porpart(i,j,k);

            uvel=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
            vvel=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
            wvel=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
            
            u_abs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);

            if(30.0*dist<ks)
            dist=ks/30.0;
            
            uplus = (1.0/kappa)*MAX(0.01,log(30.0*(dist/ks)));

        tau=(u_abs*u_abs)/pow((uplus>0.0?uplus:(1.0e20)),2.0);
        
        a->M.p[n] += (pow(p->cmu,0.75)*pow(fabs(kin(i,j,k)),0.5)*uplus)/dist;
        a->rhsvec.V[n] += (tau*u_abs)/dist;
        }
    ++n;
    }
    
}

void komega_bc::vrans_wall_law_omega(lexer *p,fdm *a,field &kin,field &eps)
{
    double uvel,vvel,wvel;
    double fac = 3.0;
    
    n=0;
    LOOP
    {
        if(a->porosity(i,j,k)>=0.99 &&  (a->porosity(i-1,j,k)<0.99 || a->porosity(i+1,j,k)<0.99 
                                        || a->porosity(i,j-1,k)<0.99 || a->porosity(i,j+1,k)<0.99)
                                        || a->porosity(i,j,k-1)<0.99 || a->porosity(i,j,k+1)<0.99)
        {
        if(p->j_dir==0)
        dist=(1.0/4.0)*(p->DXN[IP] + p->DZN[KP]);
        
        if(p->j_dir==1)
        dist=(1.0/6.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);

        eps_star = pow((kin(i,j,k)>(0.0)?(kin(i,j,k)):(0.0)),0.5) / (0.4*dist*pow(p->cmu, 0.25));
        
        a->M.p[n] += 1.0e20;
        a->rhsvec.V[n] += eps_star*1.0e20;
        }
    ++n;
    }
    
}