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

#include"nhflow_FOU.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"slice.h"
#include"patchBC_interface.h"
#include"nhflow_reconstruct_hires.h"
#include"nhflow_signal_speed.h"
#include"nhflow_flux_build_f.h"

nhflow_FOU::nhflow_FOU (lexer *p, ghostcell *ppgc, patchBC_interface *ppBC) 
{
    pgc = ppgc;
    pBC = ppBC;
    
    pflux = new nhflow_flux_build_f(p,pgc,pBC);
}

nhflow_FOU::~nhflow_FOU()
{
}

void nhflow_FOU::precalc(lexer* p, fdm_nhf* d, int ipolL, slice &eta)
{
}

void nhflow_FOU::start(lexer *&p, fdm_nhf *&d, int ipol, slice &eta)
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

void nhflow_FOU::aij_U(lexer *&p,fdm_nhf *&d, int ipol)
{
    // FOU flux 
    pflux->start_U(p,d,pgc);
    FOU(p,d,d->UHs,d->UHn,d->UHe,d->UHw);
    
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

void nhflow_FOU::aij_V(lexer *&p, fdm_nhf *&d, int ipol)
{
    // FOU flux 
    pflux->start_V(p,d,pgc);
    FOU(p,d,d->VHs,d->VHn,d->VHe,d->VHw);
    
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

void nhflow_FOU::aij_W(lexer *&p,fdm_nhf *&d, int ipol)
{
    // FOU flux 
    pflux->start_W(p,d,pgc);
    FOU(p,d,d->WHs,d->WHn,d->WHe,d->WHw);
    
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

void nhflow_FOU::aij_E(lexer *&p, fdm_nhf *&d, int ipol)
{
    // FOU flux 
    pflux->start_E(p,d,pgc);
    
    FOU_E(p,d);  // -----
    
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

void nhflow_FOU::FOU(lexer *&p,fdm_nhf *&d, double *Us, double *Un, double *Ue, double *Uw)
{    
    // FOU flux
    ULOOP
    {
        if(d->U[IJK]>=0.0)
        d->Fx[IJK] = d->Fs[IJK];
        
        if(d->U[IJK]<0.0)
        d->Fx[IJK] = d->Fn[IJK];
    }
    
    // FOU flux y-dir
    if(p->j_dir==1)
    {
    VLOOP
    {
        if(d->V[IJK]>=0.0)
        d->Fy[IJK] = d->Fe[IJK];
        
        if(d->V[IJK]<0.0)
        d->Fy[IJK] = d->Fw[IJK];
    }
    }
}

void nhflow_FOU::FOU_E(lexer *&p, fdm_nhf *&d)
{
    // FOU flux
    ULOOP
    {
        if(d->U[IJK]>=0.0)
        d->Fx[IJK] = d->Fs[IJK];
        
        else
        if(d->U[IJK]<0.0)
        d->Fx[IJK] = d->Fn[IJK];
    }
    
    // FOU flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
        if(d->V[IJK]>=0.0)
        d->Fy[IJK] = d->Fe[IJK];
        
        if(d->V[IJK]<0.0)
        d->Fy[IJK] = d->Fw[IJK];
    }
}
