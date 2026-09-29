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

#include"sflow_reconstruct_hires.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm2D.h"
#include"patchBC_interface.h"

sflow_reconstruct_hires::sflow_reconstruct_hires(lexer* p, patchBC_interface *ppBC) : dfdx(p), dfdy(p)
{
    pBC = ppBC;
}

sflow_reconstruct_hires::~sflow_reconstruct_hires()
{
}

void sflow_reconstruct_hires::reconstruct_x(lexer* p, ghostcell *pgc, fdm2D *b, slice& f, slice &fs, slice &fn)
{
    SLICELOOP4
    dfdx(i,j) = 0.0;
    
    // limited gradient
    SLICELOOP4
    WETDRY
    {
    dfdx_plus = (f(i+1,j) - f(i,j))/p->DXP[IP];
    dfdx_min  = (f(i,j) - f(i-1,j))/p->DXP[IM1];
    
    dfdx(i,j) = limiter(p,dfdx_plus,dfdx_min);
    }
    
    pgc->gcsl_start4(p,dfdx,1);
    
    // face states
    SLICELOOP1  
    {
    fs(i,j) = f(i,j)   + 0.5*p->DXN[IP]*dfdx(i,j); 
    fn(i,j) = f(i+1,j) - 0.5*p->DXN[IP1]*dfdx(i+1,j);
    }
}

void sflow_reconstruct_hires::reconstruct_y(lexer* p, ghostcell *pgc, fdm2D *b, slice& f, slice &fe, slice &fw)
{
    if(p->j_dir==1)
    {
    SLICELOOP4
    dfdy(i,j) = 0.0;
    
    // limited gradient
    SLICELOOP4
    WETDRY
    {
    dfdy_plus = (f(i,j+1) - f(i,j))/p->DYP[JP];
    dfdy_min  = (f(i,j) - f(i,j-1))/p->DYP[JM1];
    
    dfdy(i,j) = limiter(p,dfdy_plus,dfdy_min);
    }
    
    pgc->gcsl_start4(p,dfdy,1);
    
    // face states
    SLICELOOP2 
    {
    fe(i,j) = f(i,j)   + 0.5*p->DYN[JP]*dfdy(i,j); 
    fw(i,j) = f(i,j+1) - 0.5*p->DYN[JP1]*dfdy(i,j+1); 
    }
    }
}

void sflow_reconstruct_hires::reconstruct_WL(lexer* p, ghostcell *pgc, fdm2D *b)
{
    // still water depth at the faces
    SLICELOOP1
    b->dfx(i,j) = 0.5*(b->depth(i+1,j)+b->depth(i,j));

    SLICELOOP2
    b->dfy(i,j) = 0.5*(b->depth(i,j+1)+b->depth(i,j));
    
    pgc->gcsl_start1(p,b->dfx,1);
    pgc->gcsl_start2(p,b->dfy,1);
    
    // water depth at the faces
    SLICELOOP1
    {
    b->Ds(i,j) = MAX(b->ETAs(i,j) + b->dfx(i,j), p->A244);
    b->Dn(i,j) = MAX(b->ETAn(i,j) + b->dfx(i,j), p->A244);
    }
    
    SLICELOOP2
    {
    b->De(i,j) = MAX(b->ETAe(i,j) + b->dfy(i,j), p->A244);
    b->Dw(i,j) = MAX(b->ETAw(i,j) + b->dfy(i,j), p->A244);
    }
}

double sflow_reconstruct_hires::limiter(lexer *p, double v1, double v2)
{
    val=0.0;
    
    // van Leer
    if(p->A211==1)
    {
    denom = fabs(v1) + fabs(v2);
    denom = fabs(denom)>1.0e-10?denom:1.0e10;
    
    val =  (v1*fabs(v2) + fabs(v1)*v2)/denom;
    }
    
    // Superbee
    if(p->A211==2)
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
    if(p->A211==3)
    {
    r=v2/(fabs(v1)>1.0e-10?v1:1.0e20);
    
    phi = (r*r + r)/(r*r+1.0);
    
    val = 0.5*phi*(v1+v2);
    }
    
    // first-order at the wet-dry interface
    if(p->wet[IJ]==0 || p->wet[Ip1J]==0 || p->wet[Im1J]==0 || p->wet[IJp1]==0 || p->wet[IJm1]==0)
    val=0.0;
    
    return val;
}
