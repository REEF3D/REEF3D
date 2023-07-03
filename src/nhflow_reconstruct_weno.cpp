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

#include"nhflow_reconstruct_weno.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"

nhflow_reconstruct_weno::nhflow_reconstruct_weno(lexer* p, patchBC_interface *ppBC) : nhflow_gradient(p),dfdx(p), dfdy(p)
{
    pBC = ppBC;
    
    p->Darray(DFDXs,p->imax*p->jmax*(p->kmax+2));
    p->Darray(DFDXn,p->imax*p->jmax*(p->kmax+2));
}

nhflow_reconstruct_weno::~nhflow_reconstruct_weno()
{
}

void nhflow_reconstruct_weno::reconstruct_2D(lexer* p, ghostcell *pgc, fdm_nhf*, slice& f, slice &fs, slice &fn, slice &fe, slice &fw)
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
    
    pgc->gcsl_start1(p,dfdx,10);
    pgc->gcsl_start2(p,dfdy,11);
    
    // reconstruct
    SLICELOOP1  
    {
    fs(i,j) = f(i,j)   + 0.5*p->DXP[IP]*dfdx(i,j); 
    fn(i,j) = f(i+1,j) - 0.5*p->DXP[IP1]*dfdx(i+1,j);
    }

    SLICELOOP2 
    {
    fe(i,j) = f(i,j)   + 0.5*p->DYP[JP]*dfdy(i,j); 
    fw(i,j) = f(i,j+1) - 0.5*p->DYP[JP1]*dfdy(i,j+1); 
    }
    
}

void nhflow_reconstruct_weno::reconstruct_2D_x(lexer* p, ghostcell *pgc, fdm_nhf*, slice& f, slice &fs, slice &fn)
{
    // gradient
    SLICELOOP4
    {
    dfdx_plus = (f(i+1,j) - f(i,j))/p->DXP[IP];
    dfdx_min  = (f(i,j) - f(i-1,j))/p->DXP[IM1];
    
    dfdx(i,j) = limiter(dfdx_plus,dfdx_min);
    }
    
    pgc->gcsl_start1(p,dfdx,10);
    
    // reconstruct
    SLICELOOP1  
    {
    fs(i,j) = f(i,j)   + 0.5*p->DXP[IP]*dfdx(i,j); 
    fn(i,j) = f(i+1,j) - 0.5*p->DXP[IP1]*dfdx(i+1,j);
    }
}

void nhflow_reconstruct_weno::reconstruct_2D_y(lexer* p, ghostcell *pgc, fdm_nhf*, slice& f, slice &fe, slice &fw)
{
    // gradient
    SLICELOOP4
    {
    dfdy_plus = (f(i,j+1) - f(i,j))/p->DYP[JP];
    dfdy_min  = (f(i,j) - f(i,j-1))/p->DYP[JM1];
    
    dfdy(i,j) = limiter(dfdy_plus,dfdy_min);
    }
    
    pgc->gcsl_start2(p,dfdy,11);
    
    // reconstruct
    SLICELOOP2 
    {
    fe(i,j) = f(i,j)   + 0.5*p->DYP[JP]*dfdy(i,j); 
    fw(i,j) = f(i,j+1) - 0.5*p->DYP[JP1]*dfdy(i,j+1); 
    }
}

void nhflow_reconstruct_weno::reconstruct_3D_x(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fx, double *Fs, double *Fn)
{
    // gradient
    LOOP
    {
    DFDXs[IJK] = dwenox(Fx,1.0);
    DFDXn[IJK] = dwenox(Fx,-1.0);
    }

    pgc->start1V(p,DFDXs,10);
    pgc->start1V(p,DFDXn,10);
    
    // reconstruct
    ULOOP 
    {
    Fs[IJK] = (Fx[IJK]    + 0.5*p->DXP[IP]*DFDXs[IJK]); 
    Fn[IJK] = (Fx[Ip1JK]  - 0.5*p->DXP[IP1]*DFDXn[Ip1JK]);
    }
}

void nhflow_reconstruct_weno::reconstruct_3D_y(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fy, double *Fe, double *Fw)
{
    // gradient
    if(p->j_dir==1)
    {
        
    // gradient
    LOOP
    {
    DFDXs[IJK] = dwenoy(Fy,1.0);
    DFDXn[IJK] = dwenoy(Fy,-1.0);
    }

    pgc->start2V(p,DFDXs,11);
    pgc->start2V(p,DFDXn,11);
    
    // reconstruct
    VLOOP
    {
    Fe[IJK] = (Fy[IJK]    + 0.5*p->DYP[JP]*DFDXs[IJK]); 
    Fw[IJK] = (Fy[IJp1K]  - 0.5*p->DYP[JP1]*DFDXn[IJp1K]);
    }
    
    }
}

void nhflow_reconstruct_weno::reconstruct_3D_z(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fz, double *Fb, double *Ft)
{
    // gradient
    LOOP
    {
    dfdz_plus = (Fz[IJKp1] - Fz[IJK])/p->DZP[KP];
    dfdz_min  = (Fz[IJK] - Fz[IJKm1])/p->DZP[KM1];
    
    DFDXs[IJK] = limiter(dfdz_plus,dfdz_min);
    }
    
    pgc->start3V(p,DFDXs,12);
    
    // reconstruct
    LOOP 
    {
    Fb[IJK] = (Fz[IJK]    + 0.5*p->DZP[KP]*DFDXs[IJK]); 
    Ft[IJK] = (Fz[IJKp1]  - 0.5*p->DZP[KP1]*DFDXs[IJKp1]);
    }
}

double nhflow_reconstruct_weno::limiter(double v1, double v2)
{
    val=0.0;
    
    if(p->A514==1)
    {
    denom = fabs(v1) + fabs(v2);
    
    denom = fabs(denom)>1.0e-10?denom:1.0e10;
    
    val =  (v1*fabs(v2) + fabs(v1)*v2)/denom;
    }

    if(p->A514==2)
    {
    
    r=v2/(fabs(v1)>1.0e-10?v1:1.0e20);

    if(r<0.0)
    phi = 0.0;
    
    if(r>=0.0 && r<0.5)
    phi = 2.0*r;
    
    if(r>=0.5 && r<1.0)
    phi = 1.0;
    
    if(r>=1.0)
    phi = MIN(MIN(r,2.0), 2.0/(1.0+r));
    
    val = 0.5*phi*(v1+v2);
    }
    
    if(p->A514==3)
    {
    r=v2/(fabs(v1)>1.0e-10?v1:1.0e20);
    
    phi = (r*r + r)/(r*r+1.0);
    
    val = 0.5*phi*(v1+v2);
    }
    
    return val;
}


