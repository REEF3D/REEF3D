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
#include"nhflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

nhflow_f::nhflow_f(lexer *p, fdm_nhf *d, ghostcell *pgc) 
{
    margin=3;
}

nhflow_f::~nhflow_f()
{
}

void nhflow_f::ini(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow)
{
}

void nhflow_f::kinematic_fsf(lexer *p, fdm_nhf *d, double *U, double *V, double *W, slice &eta1, slice &eta2, double alpha)
{
    double wval,w_n,udetax;
    double Pval,Qval;
    double detax;
    double uvel1,uvel2;
    double zloc1,zloc2;

    
    // Kinematic Free Surface BC
    GC4LOOP
    if(p->gcb4[n][3]==6 && p->gcb4[n][4]==3)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
    
    if(p->A515==1 && p->A540==1)
    {
    Pval = 0.5*(U[IJK]+U[Im1JK]);
    Qval = 0.5*(V[IJK]+V[IJm1K]);
    
    wval = (d->eta(i,j) - d->eta_n(i,j))/(p->dt)
    
         + MAX(0.0,Pval)*((eta1(i,j)-eta1(i-1,j))/(p->DXP[IP]))
         + MIN(0.0,Pval)*((eta1(i+1,j)-eta1(i,j))/(p->DXP[IP1]))
    
         + MAX(0.0,Qval)*((eta1(i,j)-eta1(i,j-1))/(p->DYP[JP]))
         + MIN(0.0,Qval)*((eta1(i,j+1)-eta1(i,j))/(p->DYP[JP1]));
    }
         
    if(p->A515==1 && p->A540==2)
    {
    Pval = 0.5*(d->U[IJK]+d->U[Im1JK]);
    Qval = 0.5*(d->V[IJK]+d->V[IJm1K]);
        
    wval = (d->eta(i,j) - d->eta_n(i,j))/p->dt
        
         + MAX(0.0,Pval)*((d->eta(i,j)-d->eta(i-1,j))/(p->DXP[IP]))
         + MIN(0.0,Pval)*((d->eta(i+1,j)-d->eta(i,j))/(p->DXP[IP1]))
         
         + MAX(0.0,Qval)*((d->eta(i,j)-d->eta(i,j-1))/(p->DYP[JP]))
         + MIN(0.0,Qval)*((d->eta(i,j+1)-d->eta(i,j))/(p->DYP[JP1]));
    }
    
    if(p->A515==2 && p->A540==1)
    {
    wval = (d->eta(i,j) - d->eta_n(i,j))/p->dt
    
         + (eta1(i+1,j)-eta1(i-1,j))/(p->DXP[IP]+p->DXP[IP1])

         + (eta1(i,j+1)-eta1(i,j-1))/(p->DYP[JP]+p->DYP[JP1]);
    }
         
    if(p->A515==2 && p->A540==2)
    {
    wval = (d->eta(i,j) - d->eta_n(i,j))/p->dt
    
         + (d->eta(i+1,j)-d->eta(i-1,j))/(p->DXP[IP]+p->DXP[IP1])

         + (d->eta(i,j+1)-d->eta(i,j-1))/(p->DYP[JP]+p->DYP[JP1]);
    }
    
    if(p->A515==3 && p->A540==1)
    {
    Pval = 0.5*(U[IJK]+U[Im1JK]);
    Qval = 0.5*(V[IJK]+V[IJm1K]);
    
    wval = (d->eta(i,j) - d->eta_n(i,j))/(p->dt)
    
         + MAX(0.0,Pval)*(((p->A223*eta1(i,j)+(1.0-p->A223)*eta2(i,j))-(p->A223*eta1(i-1,j)+(1.0-p->A223)*eta2(i-1,j)))/(p->DXP[IP]))
         + MIN(0.0,Pval)*(((p->A223*eta1(i+1,j)+(1.0-p->A223)*eta2(i+1,j))-(p->A223*eta1(i+1,j)+(1.0-p->A223)*eta2(i+1,j)))/(p->DXP[IP1]))
    
         + MAX(0.0,Qval)*(((p->A223*eta1(i,j)+(1.0-p->A223)*eta2(i,j))-(p->A223*eta1(i,j-1)+(1.0-p->A223)*eta2(i,j-1)))/(p->DYP[JP]))
         + MIN(0.0,Qval)*(((p->A223*eta1(i,j+1)+(1.0-p->A223)*eta2(i,j+1))-(p->A223*eta1(i,j)+(1.0-p->A223)*eta2(i,j)))/(p->DYP[JP1]));
    }
         

    if(p->A515==11 && p->A540==1)
    {
    zloc1 = 0.5*p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[Ip1J]);
    zloc2 = 0.5*p->DZN[KP]*0.5*(p->sigz[Im1J]+p->sigz[IJ]);
    
    uvel1 = U[Im1JK] + (U[Im1JK]-U[Im1JKm1])/(p->DZN[KP]*0.5*(p->sigz[Im1J]+p->sigz[IJ]))*zloc1;
    uvel2 = U[IJK] + (U[IJK]-U[IJKm1])/(p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[Ip1J]))*zloc2;
    
    Pval = 0.5*(uvel1 + uvel2);
    Qval = 0.5*(d->V[IJK]+d->V[IJm1K]);
    
    wval = (d->eta(i,j) - d->eta_n(i,j))/(p->dt)
    
         + MAX(0.0,Pval)*((eta1(i,j)-eta1(i-1,j))/(p->DXP[IP]))
         + MIN(0.0,Pval)*((eta1(i+1,j)-eta1(i,j))/(p->DXP[IP1]))
         
         + MAX(0.0,Qval)*((eta1(i,j)-eta1(i,j-1))/(p->DYP[JP]))
         + MIN(0.0,Qval)*((eta1(i,j+1)-eta1(i,j))/(p->DYP[JP1]));
    }
    
    if(p->A515==11 && p->A540==2)
    {
    zloc1 = 0.5*p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[Ip1J]);
    zloc2 = 0.5*p->DZN[KP]*0.5*(p->sigz[Im1J]+p->sigz[IJ]);
    
    uvel1 = U[Im1JK] + (U[Im1JK]-U[Im1JKm1])/(p->DZN[KP]*0.5*(p->sigz[Im1J]+p->sigz[IJ]))*zloc1;
    uvel2 = U[IJK] + (U[IJK]-U[IJKm1])/(p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[Ip1J]))*zloc2;
    
    Pval = 0.5*(uvel1 + uvel2);
    Qval = 0.5*(d->V[IJK]+d->V[IJm1K]);
    
    wval = (d->eta(i,j) - d->eta_n(i,j))/(p->dt)
    
         + MAX(0.0,Pval)*((d->eta(i,j)-d->eta(i-1,j))/(p->DXP[IP]))
         + MIN(0.0,Pval)*((d->eta(i+1,j)-d->eta(i,j))/(p->DXP[IP1]))
         
         + MAX(0.0,Qval)*((d->eta(i,j)-d->eta(i,j-1))/(p->DYP[JP]))
         + MIN(0.0,Qval)*((d->eta(i,j+1)-d->eta(i,j))/(p->DYP[JP1]));
    }
         
        W[IJKp1] = wval;
        W[IJKp2] = wval;
        W[IJKp3] = wval;
    }
    
    // Kinematic Bed BC
    GC4LOOP
    if(p->gcb4[n][3]==5 && p->gcb4[n][4]==21)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
    
    Pval = 0.5*(U[IJK]+U[Im1JK]);
    Qval = 0.5*(V[IJK]+V[IJm1K]);
    

    wval = - MAX(0.0,Pval)*((d->depth(i,j)-d->depth(i-1,j))/(p->DXP[IP]))
           - MIN(0.0,Pval)*((d->depth(i+1,j)-d->depth(i,j))/(p->DXP[IP1]))
           
           - MAX(0.0,Qval)*((d->depth(i,j)-d->depth(i,j-1))/(p->DYP[JP]))
           - MIN(0.0,Qval)*((d->depth(i,j+1)-d->depth(i,j))/(p->DYP[JP1]));
    
       
    /*wval = - Pval*((d->depth(i+1,j)-d->depth(i-1,j))/(p->DXP[IP]+p->DXP[IP1]))
           - Qval*((d->depth(i,j+1)-d->depth(i,j-1))/(p->DYP[JP]+p->DYP[JP1]));*/
    

    //wval =0.0;
        
        W[IJKm1] = wval;
        W[IJKm2] = wval;
        W[IJKm3] = wval;
        
        w_n = d->wbed(i,j);
        d->wbed(i,j) = wval;
        
        d->dwdt(i,j) = (wval - w_n)/(alpha*p->dt);
    }
    
}


