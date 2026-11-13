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

#include"sflow_flux_build_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm2D.h"
#include"patchBC_interface.h"
#include"sflow_curv.h"

sflow_flux_build_f::sflow_flux_build_f(lexer *p, ghostcell *ppgc, patchBC_interface *ppBC) 
{
    pBC = ppBC;
}

sflow_flux_build_f::~sflow_flux_build_f()
{
}

void sflow_flux_build_f::start_U(lexer* p, fdm2D *b, ghostcell *pgc)
{
    // curvilinear grid: normal fluxes of the faces
    if(p->geo_curv==1)
    {
    sflow_curv::flux(p,b,1);
    return;
    }
    
    // flux x-dir
    // X 402 1: pre-balanced form with the pressure-head surface zeta = eta - ps and the
    // depth below it d + ps (ps = 0 without a pressure patch: the original form)
    const double g = fabs(p->W22);
    const int wb = (p->X10==3 && p->X402==1)?1:0;
    double ps,zs,zn;
    
    SLICELOOP1
    {
    ps = wb==1?0.5*(b->pbh(i,j) + b->pbh(i+1,j)):0.0;
    zs = b->ETAs(i,j) - ps;
    zn = b->ETAn(i,j) - ps;
    
    b->Fs(i,j) = b->UHs(i,j)*b->Us(i,j)
            + 0.5*g*zs*zs
            + g*zs*(b->dfx(i,j) + ps);
    
    b->Fn(i,j) = b->UHn(i,j)*b->Un(i,j)
            + 0.5*g*zn*zn
            + g*zn*(b->dfx(i,j) + ps);
    }
    
    // flux y-dir
    if(p->j_dir==1)
    SLICELOOP2
    {
    b->Fe(i,j) = b->Ve(i,j)*b->UHe(i,j);
    
    b->Fw(i,j) = b->Vw(i,j)*b->UHw(i,j);
    }
}

void sflow_flux_build_f::start_V(lexer* p, fdm2D *b, ghostcell *pgc)
{
    // curvilinear grid: normal fluxes of the faces
    if(p->geo_curv==1)
    {
    sflow_curv::flux(p,b,2);
    return;
    }
    
    if(p->j_dir==1)
    {
    // flux x-dir
    SLICELOOP1
    {
    b->Fs(i,j) = b->Us(i,j)*b->VHs(i,j);
    
    b->Fn(i,j) = b->Un(i,j)*b->VHn(i,j);
    }
    
    // flux y-dir
    const double g = fabs(p->W22);
    const int wb = (p->X10==3 && p->X402==1)?1:0;
    double ps,ze,zw;
    
    SLICELOOP2
    {
    ps = wb==1?0.5*(b->pbh(i,j) + b->pbh(i,j+1)):0.0;
    ze = b->ETAe(i,j) - ps;
    zw = b->ETAw(i,j) - ps;
    
    b->Fe(i,j) = b->VHe(i,j)*b->Ve(i,j) 
            + 0.5*g*ze*ze
            + g*ze*(b->dfy(i,j) + ps);
    
    b->Fw(i,j) = b->VHw(i,j)*b->Vw(i,j) 
            + 0.5*g*zw*zw 
            + g*zw*(b->dfy(i,j) + ps);
    }
    }
}

void sflow_flux_build_f::start_W(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // curvilinear grid: normal fluxes of the faces
    if(p->geo_curv==1)
    {
    sflow_curv::flux(p,b,3);
    return;
    }
    
    // flux x-dir
    SLICELOOP1
    {
    b->Fs(i,j) = b->Us(i,j)*b->WHs(i,j);
    
    b->Fn(i,j) = b->Un(i,j)*b->WHn(i,j);
    }
    
    // flux y-dir
    if(p->j_dir==1)
    SLICELOOP2
    {
    b->Fe(i,j) = b->Ve(i,j)*b->WHe(i,j);
    
    b->Fw(i,j) = b->Vw(i,j)*b->WHw(i,j);
    }
}

void sflow_flux_build_f::start_E(lexer* p, fdm2D *b, ghostcell *pgc)
{
    // curvilinear grid: normal fluxes of the faces
    if(p->geo_curv==1)
    {
    sflow_curv::flux(p,b,4);
    return;
    }
    
    // flux x-dir
    SLICELOOP1
    {
    b->Fs(i,j) = b->UHs(i,j);
    
    b->Fn(i,j) = b->UHn(i,j);
    }
    
    // flux y-dir
    if(p->j_dir==1)
    SLICELOOP2
    {
    b->Fe(i,j) = b->VHe(i,j);
    
    b->Fw(i,j) = b->VHw(i,j);
    }
}





