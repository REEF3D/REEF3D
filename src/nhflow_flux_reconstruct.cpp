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

#include"nhflow_flux_reconstruct.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"

nhflow_flux_reconstruct::nhflow_flux_reconstruct(lexer* p, patchBC_interface *ppBC) : dfdx(p), dfdy(p)
{
    pBC = ppBC;
    
    p->Darray(DFDX,p->imax*p->jmax*(p->kmax+2));
    p->Darray(DFDY,p->imax*p->jmax*(p->kmax+2));

}

nhflow_flux_reconstruct::~nhflow_flux_reconstruct()
{
}

void nhflow_flux_reconstruct::reconstruct_2D(lexer* p, ghostcell *pgc, fdm_nhf*, slice& f, slice &fs, slice &fn, slice &fe, slice &fw)
{
    // gradient
    SLICELOOP4
    {
    dfdx_plus = (f(i+1,j) - f(i,j))/p->DXP[IP];
    dfdx_min  = (f(i,j) - f(i-1,j))/p->DXP[IM1];
    
    dfdy_plus = (f(i,j+1) - f(i,j))/p->DYP[JP];
    dfdy_min  = (f(i,j) - f(i,j-1))/p->DYP[JM1];
    
    dfdx(i,j) = limiter(dfdx_plus,dfdx_min);
    dfdy(i,j) = limiter(dfdy_plus,dfdy_min);
    }
    // reconstruct
    SLICELOOP1  
    {
    fs(i,j) = f(i,j) + 0.5*p->DXP[IP]*dfdx(i,j); 
    fn(i,j) = f(i+1,j) - 0.5*p->DXP[IP1]*dfdx(i+1,j);
    }

    SLICELOOP2 
    {
    fe(i,j) = f(i,j) + 0.5*p->DYP[JP]*dfdy(i,j); 
    fw(i,j) = f(i,j+1) - 0.5*p->DYP[JP1]*dfdy(i,j+1); 
    }
}

void nhflow_flux_reconstruct::reconstruct_3D(lexer* p, ghostcell *pgc, fdm_nhf*, double *Fx, double *Fy, double *Fs, double *Fn, double *Fe, double *Fw)
{
    // gradient
    ULOOP
    {
    dfdx_plus = (Fx[Ip1JK] - Fx[IJK])/p->DXP[IP];
    dfdx_min  = (Fx[IJK] - Fx[Im1JK])/p->DXP[IM1];
    
    DFDX[IJK] = limiter(dfdx_plus,dfdx_min);
    }
    
    VLOOP
    {
    dfdy_plus = (Fy[IJp1K] - Fy[IJK])/p->DYP[JP];
    dfdy_min  = (Fy[IJK] - Fy[IJm1K])/p->DYP[JM1];

    DFDY[IJK] = limiter(dfdy_plus,dfdy_min);
    }
    
    pgc->start1V(p,DFDX,10);
    pgc->start2V(p,DFDY,11);
    
    // reconstruct
    ULOOP 
    {
    Fs[IJK] = Fx[IJK]     + 0.5*p->DXP[IP]*DFDX[IJK]; 
    Fn[IJK] = Fx[Ip1JK]   - 0.5*p->DXP[IP1]*DFDX[Ip1JK];
    }
    
    VLOOP
    {
    Fe[IJK] = Fy[IJK]     + 0.5*p->DYP[IP]*DFDY[IJK]; 
    Fw[IJK] = Fy[IJp1K]   - 0.5*p->DYP[JP1]*DFDY[IJp1K];
    }
    
    
	pgc->start1V(p,Fs,10);
    pgc->start1V(p,Fs,10);
    pgc->start2V(p,Fe,11);
    pgc->start2V(p,Fw,11);
}

double nhflow_flux_reconstruct::limiter(double v1, double v2)
{
    denom = fabs(v1) + fabs(v2);
    
    denom = fabs(denom)>1.0e-10?denom:1.0e10;
    
    val =  (v1*fabs(v2) + fabs(v1)*v2)/denom;

    return val;	
}


