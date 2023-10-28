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

#include"nhflow_reconstruct_wenograd.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"

nhflow_reconstruct_wenograd::nhflow_reconstruct_wenograd(lexer* p, patchBC_interface *ppBC) : nhflow_gradient(p),dfdxs(p), dfdxn(p), dfdye(p), dfdyw(p)
{
    pBC = ppBC;
    
    p->Darray(DFDXs,p->imax*p->jmax*(p->kmax+2));
    p->Darray(DFDXn,p->imax*p->jmax*(p->kmax+2));
}

nhflow_reconstruct_wenograd::~nhflow_reconstruct_wenograd()
{
}
  

void nhflow_reconstruct_wenograd::reconstruct_2D_x(lexer* p, ghostcell *pgc, fdm_nhf*, slice& f, slice &fs, slice &fn)
{
    // gradient
    SLICELOOP4
    {
    dfdx_plus = (f(i+1,j) - f(i,j))/p->DXP[IP];
    dfdx_min  = (f(i,j) - f(i-1,j))/p->DXP[IM1];
    
    dfdxs(i,j) = dslwenox(f,1.0);
    dfdxn(i,j) = dslwenox(f,-1.0);
    }
    
    pgc->gcsl_start1(p,dfdxs,10);
    pgc->gcsl_start1(p,dfdxn,10);
    
    // reconstruct
    SLICELOOP1  
    {
    fs(i,j) = f(i,j)   + 0.5*p->DXP[IP]*dfdxs(i,j); 
    fn(i,j) = f(i+1,j) - 0.5*p->DXP[IP1]*dfdxn(i+1,j);
    }
    
    pgc->gcsl_start1(p,fs,1);
    pgc->gcsl_start1(p,fn,1);
}

void nhflow_reconstruct_wenograd::reconstruct_2D_y(lexer* p, ghostcell *pgc, fdm_nhf*, slice& f, slice &fe, slice &fw)
{
    if(p->j_dir==1)
    {
    // gradient
    SLICELOOP4
    {   
    dfdy_plus = (f(i,j+1) - f(i,j))/p->DYP[JP];
    dfdy_min  = (f(i,j) - f(i,j-1))/p->DYP[JM1];
    
    dfdye(i,j) = dslwenoy(f,1.0);
    dfdyw(i,j) = dslwenoy(f,-1.0);
    }
    
    pgc->gcsl_start2(p,dfdye,11);
    pgc->gcsl_start2(p,dfdyw,11);
    
    // reconstruct
    SLICELOOP2 
    {
    fe(i,j) = f(i,j)   + 0.5*p->DYP[JP]*dfdye(i,j); 
    fw(i,j) = f(i,j+1) - 0.5*p->DYP[JP1]*dfdyw(i,j+1); 
    }
    
    pgc->gcsl_start2(p,fe,1);
    pgc->gcsl_start2(p,fw,1);
    }
}

void nhflow_reconstruct_wenograd::reconstruct_2D_WL(lexer* p, ghostcell *pgc, fdm_nhf *d)
{
    // water level  
    SLICELOOP1
    d->dfx(i,j) = 0.5*(d->depth(i+1,j)+d->depth(i,j));
    
    SLICELOOP2
    d->dfy(i,j) = 0.5*(d->depth(i,j+1)+d->depth(i,j));
    
    pgc->gcsl_start1(p,d->dfx,1);
    pgc->gcsl_start2(p,d->dfy,1);

    SLICELOOP1
    {
    d->Ds(i,j) = MAX(d->ETAs(i,j) + 0.5*(d->depth(i+1,j)+d->depth(i,j)), p->A544);
    d->Dn(i,j) = MAX(d->ETAn(i,j) + 0.5*(d->depth(i+1,j)+d->depth(i,j)), p->A544);
    }
    
    SLICELOOP2
    {
    d->De(i,j) = MAX(d->ETAe(i,j)  + 0.5*(d->depth(i,j+1)+d->depth(i,j)), p->A544);
    d->Dw(i,j) = MAX(d->ETAw(i,j)  + 0.5*(d->depth(i,j+1)+d->depth(i,j)), p->A544);
    }
    
    pgc->gcsl_start1(p,d->Ds,1);
    pgc->gcsl_start1(p,d->Dn,1);
    pgc->gcsl_start2(p,d->De,1);
    pgc->gcsl_start2(p,d->Dw,1);
}

void nhflow_reconstruct_wenograd::reconstruct_3D_x(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fx, double *Fs, double *Fn)
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

void nhflow_reconstruct_wenograd::reconstruct_3D_y(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fy, double *Fe, double *Fw)
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

void nhflow_reconstruct_wenograd::reconstruct_3D_z(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fz, double *Fb, double *Ft)
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
    WLOOP 
    {
    Fb[IJK] = (Fz[IJK]    + 0.5*p->DZN[KP]*DFDXs[IJK]); 
    Ft[IJK] = (Fz[IJKp1]  - 0.5*p->DZN[KP1]*DFDXs[IJKp1]);
    }
}

double nhflow_reconstruct_wenograd::limiter(double v1, double v2)
{
    val=0.0;
    
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

    
    return val;
}


