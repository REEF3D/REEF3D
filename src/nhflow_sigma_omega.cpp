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

#include"nhflow_sigma.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

#define WLVL (fabs(d->WL(i,j))>0.00005?d->WL(i,j):1.0e20)
#define HX (fabs(d->hx(i,j))>1.0e-20?d->hx(i,j):1.0e20)
#define HXP (fabs(0.5*(d->WL(i,j)+d->WL(i+1,j)))>1.0e-20?0.5*(d->WL(i,j)+d->WL(i+1,j)):1.0e20)
#define HY (fabs(d->hy(i,j))>1.0e-20?d->hy(i,j):1.0e20)

void nhflow_sigma::omega_update(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, double *U, double *V, double *W)
{ 
    double wval,Pval,Qval,Rval,fac;
    
    if(p->A517==1)
    FLOOP
    {
        if(U[IJK]>=0.0)
        Pval=0.5*(U[Im1JK] + U[IJK]);
            
        if(U[IJK]<0.0)
        Pval=0.5*(U[IJK] + U[Ip1JK]);
        
        
        if(V[IJK]>=0.0)
        Qval=0.5*(V[IJm1K] + V[IJK]);
            
        if(V[IJK]<0.0)
        Qval=0.5*(V[IJK] + V[IJp1K]);
        
        
        if(W[IJK]>=0.0)
        Rval=0.5*(W[IJKm1] + W[IJK]);
            
        if(W[IJK]<0.0)
        Rval=0.5*(W[IJK] + W[IJKp1]);

    // omega
        d->omegaF[FIJK] =  d->WL(i,j)*(p->sigt[FIJK]
                        
                        +  Pval*p->sigx[FIJK]
                        
                        +  Qval*p->sigy[FIJK]
                        
                        +  Rval*p->sigz[IJ]);
    }
    
    if(p->A517==2)
    {
        FLOOP
        {
            
            Pval = 0.5*(U[IJK]+U[IJKm1]);
            Qval = 0.5*(V[IJK]+V[IJKm1]);
            Rval = 0.5*(W[IJK]+W[IJKm1]);
        
            // omega
            d->omegaF[FIJK] =  WL(i,j)*(p->sigt[FIJK]
                            
                            +  Pval*p->sigx[FIJK]
                            
                            +  Qval*p->sigy[FIJK]
                            
                            +  Rval*p->sigz[IJ]);
        }
        
        
        ILOOP 
        JLOOP
        if(p->wet[IJ]==0 || p->wet[Im1J]==0 || p->wet[Ip1J]==0 || p->wet[IJm1]==0 || p->wet[IJp1]==0
        || p->wet[Im2J]==0 || p->wet[Ip2J]==0 || p->wet[IJm2]==0 || p->wet[IJp2]==0)
        {
        FKLOOP 
        FPCHECK 
        d->omegaF[FIJK] = 0.0;  
        
        KLOOP 
        PCHECK 
        d->omegaF[FIJKp1] =   d->omegaF[FIJK]
                        
                        - p->DZN[KP]*(d->detadt(i,j) 
                        
                        + (d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP]  + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir);
        }
        
    }
    
    if(p->A517==3)
    {
        FLOOP
        d->omegaF[FIJK] = 0.0;
        
        
        LOOP
        {
        d->omegaF[FIJKp1] =   d->omegaF[FIJK]
                            
                            - p->DZN[KP]*(d->detadt(i,j) 
                            
                            + (d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP]  + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir);
        }
    }
    
    if(p->A517==4)
    {
        FLOOP
        d->omegaF[FIJK] = 0.0;
        
        LOOP
        {
        d->omegaF[FIJKp1] =   d->omegaF[FIJK]
                            
                            - p->DZN[KP]*(d->detadt_n(i,j) 
                            
                            + (d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP]  + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir);
        }
    }

      
    GC4LOOP
    if(p->gcb4[n][3]==6 && p->gcb4[n][4]==3)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
    
        k+=1;
        d->omegaF[FIJK] =  0.0;
        d->omegaF[FIJKp1] =  0.0;
        d->omegaF[FIJKp2] =  0.0;
        d->omegaF[FIJKp3] =  0.0;
        
    }
    
    GC4LOOP
    if(p->gcb4[n][3]==5 && p->gcb4[n][4]==21)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
        
        d->omegaF[FIJK] =  0.0;
        d->omegaF[FIJKm1] =  0.0;
        d->omegaF[FIJKm2] =  0.0;
        d->omegaF[FIJKm3] =  0.0;
        
    }
    

    FLOOP
    if(p->wet[IJ]==0)
    d->omegaF[FIJK] = 0.0;
    
    pgc->start7S(p,d->omegaF,17);
    
}




