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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_reinidisc_fsf.h"

nhflow_forcing::nhflow_forcing(lexer *p) : epsi(1.6)
{
    if(p->A561>0 || p->A564>0)
    {
    p->Iarray(IO,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(CL,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(CR,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(FX,p->imax*p->jmax*(p->kmax+2));
    p->Darray(FY,p->imax*p->jmax*(p->kmax+2));
    p->Darray(FZ,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(FRK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(dt,p->imax*p->jmax*(p->kmax+2));
    p->Darray(L,p->imax*p->jmax*(p->kmax+2));
    
    prdisc = new nhflow_reinidisc_fsf(p);
    }
}

nhflow_forcing::~nhflow_forcing()
{
}

void nhflow_forcing::forcing(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, double *UH, double *VH, double *WH, slice &WL)
{
    if(p->A561>0 || p->A564>0)
    {
    // update direct forcing function
    ray_cast(p, d, pgc);
    reini_RK2(p, d, pgc, d->SOLID);
    
    // update Heaviside
    LOOP
    d->FHB[IJK] = 0.0;

    pgc->start5V(p,d->FHB,1);
    
    LOOP
    {
        H = Hsolidface(p,d,0,0,0);
        d->FHB[IJK] = min(d->FHB[IJK] + H, 1.0); 
        
        uf = 0.0;
        
        d->FX[IJK] += d->FHB[IJK]*(uf*WL(i,j) - UH[IJK])/(alpha*p->dt);  
    }
    
    LOOP
    {
        vf = 0.0;

        d->FY[IJK] += d->FHB[IJK]*(vf*WL(i,j) - VH[IJK])/(alpha*p->dt);  
    }
    
    LOOP
    {
        wf = 0.0;

        d->FZ[IJK] += d->FHB[IJK]*(wf*WL(i,j) - WH[IJK])/(alpha*p->dt);  
    }
    	
    pgc->start5V(p,d->FX,1);
    pgc->start5V(p,d->FY,1);
    pgc->start5V(p,d->FZ,1);
    pgc->start5V(p,d->FHB,1);
    
// Calculate forcing fields
    
    // add forcing term to RHS
    
    LOOP
    {
        UH[IJK] += alpha*p->dt*CPORNH*d->FX[IJK];
        
        if(p->count<10)
        d->maxF = MAX(fabs(alpha*CPORNH*d->FX[IJK]), d->maxF);
        
        p->fbmax = MAX(fabs(alpha*CPORNH*d->FX[IJK]), p->fbmax);
    }
    
    LOOP
    {
        VH[IJK] += alpha*p->dt*CPORNH*d->FY[IJK];
        
        if(p->count<10)
        d->maxG = MAX(fabs(alpha*CPORNH*d->FY[IJK]), d->maxG);
        
        p->fbmax = MAX(fabs(alpha*CPORNH*d->FY[IJK]), p->fbmax);
    }
    
    LOOP
    {
        WH[IJK] += alpha*p->dt*CPORNH*d->FZ[IJK];
        
        if(p->count<10)
        d->maxH = MAX(fabs(alpha*CPORNH*d->FZ[IJK]), d->maxH);
        
        p->fbmax = MAX(fabs(alpha*CPORNH*d->FZ[IJK]), p->fbmax);
    }
    
    }
}

void nhflow_forcing::forcing_ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->A561>0 || p->A564>0)
    {
    LOOP
    p->ZSP[IJK]  = p->ZP[KP]*d->WL(i,j) + d->bed(i,j);
    
    pgc->start5V(p,p->ZSP,1);
    
    objects_create(p, pgc);
    
    geometry_refinement(p,pgc);
    
    ray_cast(p, d, pgc);
    
    reini_RK2(p, d, pgc, d->SOLID);
    }
}