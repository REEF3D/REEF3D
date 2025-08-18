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

#include"nhflow_reconstruct_void.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"

nhflow_reconstruct_void::nhflow_reconstruct_void(lexer* p, patchBC_interface *ppBC) : nhflow_gradient(p),dfdx(p), dfdy(p)
{
    pBC = ppBC;
    
    p->Darray(DFDX,p->imax*p->jmax*(p->kmax+2));
}

nhflow_reconstruct_void::~nhflow_reconstruct_void()
{
}

void nhflow_reconstruct_void::reconstruct_2D_x(lexer* p, ghostcell *pgc, fdm_nhf*, slice& f, slice &fs, slice &fn)
{

}

void nhflow_reconstruct_void::reconstruct_2D_y(lexer* p, ghostcell *pgc, fdm_nhf*, slice& f, slice &fe, slice &fw)
{

}

void nhflow_reconstruct_void::reconstruct_2D_WL(lexer* p, ghostcell *pgc, fdm_nhf *d)
{
    
}

void nhflow_reconstruct_void::reconstruct_3D_x(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fx, double *Fs, double *Fn)
{

}

void nhflow_reconstruct_void::reconstruct_3D_y(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fy, double *Fe, double *Fw)
{
    
}

void nhflow_reconstruct_void::reconstruct_3D_z(lexer* p, ghostcell *pgc, fdm_nhf *d, double *Fz, double *Fb, double *Ft)
{

}

double nhflow_reconstruct_void::limiter(double v1, double v2)
{
    val=0.0;
    
    // van Leer
    if(p->A514==1)
    {
    denom = fabs(v1) + fabs(v2);
    
    denom = fabs(denom)>1.0e-10?denom:1.0e10;
    
    val =  (v1*fabs(v2) + fabs(v1)*v2)/denom;
    }
    
    // Superbee
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
    
    // van Albada
    if(p->A514==3)
    {
    r=v2/(fabs(v1)>1.0e-10?v1:1.0e20);
    
    phi = (r*r + r)/(r*r+1.0);
    
    val = 0.5*phi*(v1+v2);
    }
    
    if(p->wet[IJ]==0 || p->wet[Ip1J]==0 || p->wet[Im1J]==0 || p->wet[IJp1]==0 || p->wet[IJm1]==0)
    val=0.0;
    
    return val;
}


