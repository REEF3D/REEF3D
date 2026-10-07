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

#include"nhflow_kepsilon_bc.h"
#include"nhflow_wall.h"
#include"fdm_nhf.h"
#include"lexer.h"

// inflow ghost cell with prescribed turbulence profile (B 60 >= 1, see nhflow_rans_io::inflow, which writes the i-1 ghosts): use the ghost value
#define TURBIN(X) (p->B60>=1 && p->IO[X]==1 && p->DF[X]>0)
 
nhflow_kepsilon_bc::nhflow_kepsilon_bc(lexer *p) : roughness(p)
{
    kappa=0.4;
}

nhflow_kepsilon_bc::~nhflow_kepsilon_bc()
{
}

void nhflow_kepsilon_bc::bckepsilon_start(lexer *p, fdm_nhf *d, double *KIN, double *EPS, int gcval)
{
	if(gcval==20)
    wall_law_kin(p,d,KIN,EPS);
        
	if(gcval==30)
	wall_law_epsilon(p,d,KIN,EPS);
}

void nhflow_kepsilon_bc::wall_law_kin(lexer *p, fdm_nhf *d, double *KIN, double *EPS)
{
    // wall function for k at the nearest wall that has one (nhflow_wall.h): production tau*u/y and
    // dissipation cmu^3/4 k^1/2 u+/y (implicit); y >= ks/30
    double ut;
    
    count=0;
    LOOP
    {
        if(nhflow_turb_wall(p,d,i,j,k,dist,ks,ut)==1)
        {
        if(30.0*dist<ks)
        dist=ks/30.0;
        
        uplus = (1.0/kappa)*MAX(1.0,log(30.0*(dist/ks)));
        
        tau = (ut*ut)/(uplus*uplus);
        
        d->M.p[count] += (pow(p->cmu,0.75)*pow(fabs(KIN[IJK]),0.5)*uplus)/dist;
        d->rhsvec.V[count] += (tau*ut)/dist;
        }
    ++count;
    }
}

void nhflow_kepsilon_bc::wall_law_epsilon(lexer *p, fdm_nhf *d, double *KIN, double *EPS)
{
    // wall value at the nearest wall that has a wall function (nhflow_wall.h)
    double ut;
    
    count=0;
    LOOP
    {
        if(nhflow_turb_wall(p,d,i,j,k,dist,ks,ut)==1)
        {
        eps_star = (pow(p->cmu, 0.75)*pow(MAX(KIN[IJK],0.0),1.5)) / (kappa*dist);
        
        d->M.p[count] += 1.0e20;
        d->rhsvec.V[count] += eps_star*1.0e20;
        }
    ++count;
    }
}

void nhflow_kepsilon_bc::bckin_matrix(lexer *p, fdm_nhf *d, double *KIN, double *EPS)
{
    // sharp thin bodies (X 330 'mobility sharp'): zero gradient across the body (the wall function, nhflow_wall.h)
    if(d->thinbody!=nullptr)
    d->thinbody->matrix_walls(p,d,KIN);

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
            if(p->flag4[IJK]>0 && p->DF[IJK]>0)
            {
            if((p->flag4[Im1JK]<0 || p->DF[Im1JK]<0))// && inflow==0)
            {
            if(TURBIN(Im1JK)) d->rhsvec.V[n] -= d->M.s[n]*KIN[Im1JK];   // discharge inflow profile (Dirichlet)
            else d->M.p[n] += d->M.s[n];   // zero gradient, implicit (was lagged)
            d->M.s[n] = 0.0;
            }
            
            if((p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0))// && outflow==0)
            {
            d->M.p[n] += d->M.n[n];
            d->M.n[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJm1K]<0 || p->DF[IJm1K]<0)
            {
            d->M.p[n] += d->M.e[n];
            d->M.e[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJp1K]<0 || p->DF[IJp1K]<0)
            {
            d->M.p[n] += d->M.w[n];
            d->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0 || p->DF[IJKm1]<0)
            {
            d->M.p[n] += d->M.b[n];
            d->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0 || p->DF[IJKp1]<0)
            {
            d->M.p[n] += d->M.t[n];
            d->M.t[n] = 0.0;
            }
            }

        ++n;
        }
        
        // wet/dry front: a dry neighbour column is not a wall, zero gradient (implicit) instead of k = eps = 0
        n=0;
        LOOP
        {
            if(p->flag4[IJK]>0 && p->DF[IJK]>0 && p->wet[IJ]==1)
            {
            if(p->flag4[Im1JK]>0 && p->DF[Im1JK]>0 && p->wet[Im1J]==0)
            {
            d->M.p[n] += d->M.s[n];
            d->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]>0 && p->DF[Ip1JK]>0 && p->wet[Ip1J]==0)
            {
            d->M.p[n] += d->M.n[n];
            d->M.n[n] = 0.0;
            }
            
            if(p->j_dir==1 && p->flag4[IJm1K]>0 && p->DF[IJm1K]>0 && p->wet[IJm1]==0)
            {
            d->M.p[n] += d->M.e[n];
            d->M.e[n] = 0.0;
            }
            
            if(p->j_dir==1 && p->flag4[IJp1K]>0 && p->DF[IJp1K]>0 && p->wet[IJp1]==0)
            {
            d->M.p[n] += d->M.w[n];
            d->M.w[n] = 0.0;
            }
            }
        ++n;
        }
        
        n=0;
        LOOP
        {
            if(p->DF[IJK]<0 || p->wet[IJ]==0)
            {   
            KIN[IJK] = 0.0;
            
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

void nhflow_kepsilon_bc::bcepsilon_matrix(lexer *p, fdm_nhf *d, double *KIN, double *EPS)
{
    // sharp thin bodies (X 330 'mobility sharp'): zero gradient across the body (the wall function, nhflow_wall.h)
    if(d->thinbody!=nullptr)
    d->thinbody->matrix_walls(p,d,EPS);

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
            if(p->flag4[IJK]>0 && p->DF[IJK]>0)
            {
            // s
            if(p->flag4[Im1JK]<0)// && inflow==0)
            {
            if(TURBIN(Im1JK)) d->rhsvec.V[n] -= d->M.s[n]*EPS[Im1JK];   // discharge inflow profile (Dirichlet)
            else d->M.p[n] += d->M.s[n];   // zero gradient, implicit (was lagged)
            d->M.s[n] = 0.0;
            }
            
            if(p->DF[Im1JK]<0)
            {
            if(TURBIN(Im1JK)) d->rhsvec.V[n] -= d->M.s[n]*EPS[Im1JK];   // discharge inflow profile (Dirichlet)
            else d->M.p[n] += d->M.s[n];   // zero gradient, implicit (was lagged)
            d->M.s[n] = 0.0;
            }
            
            // n
            if(p->flag4[Ip1JK]<0)// && outflow==0)
            {
            d->M.p[n] += d->M.n[n];
            d->M.n[n] = 0.0;
            }
            
            if(p->DF[Ip1JK]<0)
            {
            d->M.p[n] += d->M.n[n];
            d->M.n[n] = 0.0;
            }
            
            // e
            if(p->j_dir==1)
            if(p->flag4[IJm1K]<0)
            {
            d->M.p[n] += d->M.e[n];
            d->M.e[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->DF[IJm1K]<0)
            {
            d->M.p[n] += d->M.e[n];
            d->M.e[n] = 0.0;
            }
            
            // w
            if(p->j_dir==1)
            if(p->flag4[IJp1K]<0)
            {
            d->M.p[n] += d->M.w[n];
            d->M.w[n] = 0.0;
            }
            
            if(p->j_dir==1)
            if(p->DF[IJp1K]<0)
            {
            d->M.p[n] += d->M.w[n];
            d->M.w[n] = 0.0;
            }
            
            // b
            if(p->flag4[IJKm1]<0)
            {
            d->M.p[n] += d->M.b[n];
            d->M.b[n] = 0.0;
            }
            
            if(p->DF[IJKm1]<0)
            {
            d->M.p[n] += d->M.b[n];
            d->M.b[n] = 0.0;
            }
            
            // t
            if(p->flag4[IJKp1]<0)
            {
            d->M.p[n] += d->M.t[n];
            d->M.t[n] = 0.0;
            }
            
            if(p->DF[IJKp1]<0)
            {
            d->M.p[n] += d->M.t[n];
            d->M.t[n] = 0.0;
            }
            }

        ++n;
        }
        
        // turn off inside direct forcing body
        // wet/dry front: a dry neighbour column is not a wall, zero gradient (implicit) instead of k = eps = 0
        n=0;
        LOOP
        {
            if(p->flag4[IJK]>0 && p->DF[IJK]>0 && p->wet[IJ]==1)
            {
            if(p->flag4[Im1JK]>0 && p->DF[Im1JK]>0 && p->wet[Im1J]==0)
            {
            d->M.p[n] += d->M.s[n];
            d->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]>0 && p->DF[Ip1JK]>0 && p->wet[Ip1J]==0)
            {
            d->M.p[n] += d->M.n[n];
            d->M.n[n] = 0.0;
            }
            
            if(p->j_dir==1 && p->flag4[IJm1K]>0 && p->DF[IJm1K]>0 && p->wet[IJm1]==0)
            {
            d->M.p[n] += d->M.e[n];
            d->M.e[n] = 0.0;
            }
            
            if(p->j_dir==1 && p->flag4[IJp1K]>0 && p->DF[IJp1K]>0 && p->wet[IJp1]==0)
            {
            d->M.p[n] += d->M.w[n];
            d->M.w[n] = 0.0;
            }
            }
        ++n;
        }
        
        n=0;
        LOOP
        {
            if(p->DF[IJK]<0 || p->wet[IJ]==0)
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
