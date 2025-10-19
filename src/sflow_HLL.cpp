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

#include"nhflow_HLL.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm2D.h"
#include"slice.h"
#include"patchBC_interface.h"
#include"nhflow_reconstruct_hires.h"
#include"nhflow_signal_speed.h"
#include"nhflow_flux_build_f.h"

nhflow_HLL::nhflow_HLL (lexer *p, ghostcell *ppgc, patchBC_interface *ppBC) 
{
    pgc = ppgc;
    pBC = ppBC;
    
    pflux = new nhflow_flux_build_f(p,pgc,pBC);
}

nhflow_HLL::~nhflow_HLL()
{
}

void nhflow_HLL::precalc(lexer* p, fdm2D *b, int ipolL, slice &eta)
{
}

void nhflow_HLL::start(lexer *&p, fdm2D *&b, int ipol, slice &eta)
{
    if(ipol==1)
    aij_U(p,d,1);

    if(ipol==2 && p->j_dir==1)
    aij_V(p,d,2);

    if(ipol==3)
    aij_W(p,d,3);
    
    if(ipol==4)
    aij_E(p,d,4);
}

void nhflow_HLL::aij_U(lexer *&p,fdm2D *&b, int ipol)
{
    // HLL flux 
    pflux->start_U(p,d,pgc);
    HLL(p,d,b->UHs,b->UHn,b->UHe,b->UHw);
    
    pgc->start1V(p,b->Fx,10);
    pgc->start2V(p,b->Fy,10);
    
    LOOP
    WETDRY
    {
    b->F(i,j) -= ((b->Fx(i,j) - b->Fx[Im1JK])/p->DXN[IP] 
                + (b->Fy(i,j) - b->Fy[IJm1K])/p->DYN[JP]*p->y_dir;
    }    
}

void nhflow_HLL::aij_V(lexer *&p, fdm2D *&b, int ipol)
{
    // HLL flux 
    pflux->start_V(p,d,pgc);
    HLL(p,d,b->VHs,b->VHn,b->VHe,b->VHw);
    
    pgc->start1V(p,b->Fx,11);
    pgc->start2V(p,b->Fy,11);
    
    LOOP
    WETDRY
    {
    b->G(i,j) -= ((b->Fx(i,j) - b->Fx[Im1JK])/p->DXN[IP] 
                + (b->Fy(i,j) - b->Fy[IJm1K])/p->DYN[JP]*p->y_dir;
    }    
}

void nhflow_HLL::aij_W(lexer *&p,fdm2D *&b, int ipol)
{
    // HLL flux 
    pflux->start_W(p,d,pgc);
    HLL(p,d,b->WHs,b->WHn,b->WHe,b->WHw);
    
    pgc->start1V(p,b->Fx,12);
    pgc->start2V(p,b->Fy,12);
    
    LOOP
    WETDRY
    {
    b->H(i,j) -= ((b->Fx(i,j) - b->Fx[Im1JK])/p->DXN[IP] 
                + (b->Fy(i,j) - b->Fy[IJm1K])/p->DYN[JP]*p->y_dir;
    }    
}

void nhflow_HLL::aij_E(lexer *&p, fdm2D *&b, int ipol)
{
    // HLL flux 
    pflux->start_E(p,d,pgc);
    
    HLL_E(p,d);  // -----
    
    LOOP
    WETDRY
    {
    if(p->wet[Ip1J]==0)
    b->Fx(i,j) = 0.0;
    
    if(p->wet[Im1J]==0)
    b->Fx[Im1JK] = 0.0;
    
    if(p->wet[IJp1]==0)
    b->Fy(i,j) = 0.0;
    
    if(p->wet[IJm1]==0)
    b->Fy[IJm1K] = 0.0;
    }
    
    pgc->start1V(p,b->Fx,14);
    pgc->start2V(p,b->Fy,14); 
}

void nhflow_HLL::HLL(lexer *&p,fdm2D *&b, double *Us, double *Un, double *Ue, double *Uw)
{    
    // HLL flux
    ULOOP
    {
        if(b->Ss(i,j)>=0.0)
        b->Fx(i,j) = b->Fs(i,j);
        
        else
        if(b->Sn(i,j)<=0.0)
        b->Fx(i,j) = b->Fn(i,j);
        
        else
        {
        denom = b->Sn(i,j)-b->Ss(i,j);
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        b->Fx(i,j) = (b->Sn(i,j)*b->Fs(i,j) - b->Ss(i,j)*b->Fn(i,j) + b->Sn(i,j)*b->Ss(i,j)*(Un(i,j) - Us(i,j)))/denom;
        }
    }
    
    // HLL flux y-dir
    if(p->j_dir==1)
    {
    VLOOP
    {
        if(b->Se(i,j)>=0.0)
        b->Fy(i,j) = b->Fe(i,j);
        
        else
        if(b->Sw(i,j)<=0.0)
        b->Fy(i,j) = b->Fw(i,j);
        
        else
        {
        denom = b->Sw(i,j)-b->Se(i,j);
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        b->Fy(i,j) = (b->Sw(i,j)*b->Fe(i,j) - b->Se(i,j)*b->Fw(i,j) + b->Sw(i,j)*b->Se(i,j)*(Uw(i,j) - Ue(i,j)))/denom;
        }
    }
    }
}

void nhflow_HLL::HLL_E(lexer *&p, fdm2D *&b)
{
    // HLL flux
    ULOOP
    {
        if(b->Ss(i,j)>=0.0)
        b->Fx(i,j) = b->Fs(i,j);
        
        else
        if(b->Sn(i,j)<=0.0)
        b->Fx(i,j) = b->Fn(i,j);
        
        else
        {
        denom = b->Sn(i,j)-b->Ss(i,j);
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        b->Fx(i,j) = (b->Sn(i,j)*b->Fs(i,j) - b->Ss(i,j)*b->Fn(i,j) + b->Sn(i,j)*b->Ss(i,j)*(b->Dn(i,j) - b->Ds(i,j)))/denom;
        }
    }
    
    // HLL flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
        if(b->Se(i,j)>=0.0)
        b->Fy(i,j) = b->Fe(i,j);
        
        else
        if(b->Sw(i,j)<=0.0)
        b->Fy(i,j) = b->Fw(i,j);
        
        else
        {
        denom = b->Sw(i,j)-b->Se(i,j);
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        b->Fy(i,j) = (b->Sw(i,j)*b->Fe(i,j) - b->Se(i,j)*b->Fw(i,j) + b->Sw(i,j)*b->Se(i,j)*(b->Dw(i,j) - b->De(i,j)))/denom;
        }
    }
}
