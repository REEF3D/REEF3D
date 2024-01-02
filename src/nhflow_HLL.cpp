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

#include"nhflow_HLL.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
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

void nhflow_HLL::precalc(lexer* p, fdm_nhf* d, int ipolL, slice &eta)
{
}

void nhflow_HLL::start(lexer* p, fdm_nhf *&d, int ipol, slice &eta)
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

double nhflow_HLL::aij_U(lexer *p,fdm_nhf *&d, int ipol)
{
    // HLL flux 
    pflux->start_U(p,d,pgc);
    HLL(p,d,d->UHs,d->UHn,d->UHe,d->UHw);
    
    LOOP
    WETDRY
    {
    if(p->wet[IJp1]==0 && p->flag2[IJp1K]>0)
    d->Fy[IJK] = 0.0;
    
    if(p->wet[IJm1]==0 && p->flag2[IJm1K]>0)
    d->Fy[IJm1K] = 0.0;
    }
    
    pgc->start1V(p,d->Fx,10);
    pgc->start2V(p,d->Fy,10);
    pgc->start3V(p,d->Fz,10);
    
    LOOP
    WETDRY
    {
    d->F[IJK] -= ((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP] 
                + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir
                + (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    }    
}

double nhflow_HLL::aij_V(lexer* p, fdm_nhf *&d, int ipol)
{
    // HLL flux 
    pflux->start_V(p,d,pgc);
    HLL(p,d,d->VHs,d->VHn,d->VHe,d->VHw);
    
    LOOP
    WETDRY
    {
    if(p->wet[Ip1J]==0 && p->flag1[Ip1JK]>0)
    d->Fx[IJK] = 0.0;
    
    if(p->wet[Im1J]==0 && p->flag1[Im1JK]>0)
    d->Fx[Im1JK] = 0.0;
    }
    
    pgc->start1V(p,d->Fx,11);
    pgc->start2V(p,d->Fy,11);
    pgc->start3V(p,d->Fz,11);
    
    LOOP
    WETDRY
    {
    d->G[IJK] -= ((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP] 
                + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir
                + (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    }    
}

double nhflow_HLL::aij_W(lexer *p,fdm_nhf *&d, int ipol)
{
    // HLL flux 
    pflux->start_W(p,d,pgc);
    HLL(p,d,d->WHs,d->WHn,d->WHe,d->WHw);
    
    pgc->start1V(p,d->Fx,12);
    pgc->start2V(p,d->Fy,12);
    pgc->start3V(p,d->Fz,12);
    
    LOOP
    WETDRY
    {
    d->H[IJK] -= ((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP] 
                + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir
                + (d->Fz[IJK] - d->Fz[IJKm1])/p->DZN[KP]);
    }    
}

double nhflow_HLL::aij_E(lexer *p, fdm_nhf *&d, int ipol)
{
    // HLL flux 
    pflux->start_E(p,d,pgc);
    
    HLL_E(p,d);  // -----
    
    LOOP
    WETDRY
    {
    if(p->wet[Ip1J]==0)
    d->Fx[IJK] = 0.0;
    
    if(p->wet[Im1J]==0)
    d->Fx[Im1JK] = 0.0;
    
    if(p->wet[IJp1]==0)
    d->Fy[IJK] = 0.0;
    
    if(p->wet[IJm1]==0)
    d->Fy[IJm1K] = 0.0;
    }
    
    pgc->start1V(p,d->Fx,14);
    pgc->start2V(p,d->Fy,14); 
}

double nhflow_HLL::HLL(lexer *p,fdm_nhf *&d, double *Us, double *Un, double *Ue, double *Uw)
{
    // HLL flux
    ULOOP
    {
        if(d->Ss[IJK]>=0.0)
        d->Fx[IJK] = d->Fs[IJK];
        
        else
        if(d->Sn[IJK]<=0.0)
        d->Fx[IJK] = d->Fn[IJK];
        
        else
        {
        denom = d->Sn[IJK]-d->Ss[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        d->Fx[IJK] = (d->Sn[IJK]*d->Fs[IJK] - d->Ss[IJK]*d->Fn[IJK] + d->Sn[IJK]*d->Ss[IJK]*(Un[IJK] - Us[IJK]))/denom;
        }
    }
    
    // HLL flux y-dir
    if(p->j_dir==1)
    {
    VLOOP
    {
        if(d->Se[IJK]>=0.0)
        d->Fy[IJK] = d->Fe[IJK];
        
        else
        if(d->Sw[IJK]<=0.0)
        d->Fy[IJK] = d->Fw[IJK];
        
        else
        {
        denom = d->Sw[IJK]-d->Se[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        d->Fy[IJK] = (d->Sw[IJK]*d->Fe[IJK] - d->Se[IJK]*d->Fw[IJK] + d->Sw[IJK]*d->Se[IJK]*(Uw[IJK] - Ue[IJK]))/denom;
        }
    }
    }
}

double nhflow_HLL::HLL_E(lexer *p, fdm_nhf *&d)
{
    // HLL flux
    ULOOP
    {
        if(d->Ss[IJK]>=0.0)
        d->Fx[IJK] = d->Fs[IJK];
        
        else
        if(d->Sn[IJK]<=0.0)
        d->Fx[IJK] = d->Fn[IJK];
        
        else
        {
        denom = d->Sn[IJK]-d->Ss[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        d->Fx[IJK] = (d->Sn[IJK]*d->Fs[IJK] - d->Ss[IJK]*d->Fn[IJK] + d->Sn[IJK]*d->Ss[IJK]*(d->Dn(i,j) - d->Ds(i,j)))/denom;
        }
    }
    
    // HLL flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
        if(d->Se[IJK]>=0.0)
        d->Fy[IJK] = d->Fe[IJK];
        
        else
        if(d->Sw[IJK]<=0.0)
        d->Fy[IJK] = d->Fw[IJK];
        
        else
        {
        denom = d->Sw[IJK]-d->Se[IJK];
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        d->Fy[IJK] = (d->Sw[IJK]*d->Fe[IJK] - d->Se[IJK]*d->Fw[IJK] + d->Sw[IJK]*d->Se[IJK]*(d->Dw(i,j) - d->De(i,j)))/denom;
        }
    }
}
