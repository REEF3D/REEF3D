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

#include"nhflow_flux_build_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"slice.h"
#include"patchBC_interface.h"
#include"vrans.h"
#include"nhflow_flux_face.h"

nhflow_flux_build_f::nhflow_flux_build_f(lexer *p, ghostcell *ppgc, patchBC_interface *ppBC) 
{
    pBC = ppBC;
}

nhflow_flux_build_f::~nhflow_flux_build_f()
{
}

void nhflow_flux_build_f::start_U(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
    // flux x-dir
    ULOOP
    {
    d->Fs[IJK] = nhflow_face::U_s(p,d,i,j,k);
    d->Fn[IJK] = nhflow_face::U_n(p,d,i,j,k);
    }
    
    // flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
    d->Fe[IJK] = nhflow_face::U_e(p,d,i,j,k);
    d->Fw[IJK] = nhflow_face::U_w(p,d,i,j,k);
    }
    
    // flux z-dir
    nhflow_face::zflux(p,d,d->Ub,d->Ut);
}

void nhflow_flux_build_f::start_V(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->j_dir==1)
    {
    // flux x-dir
    ULOOP
    {
    d->Fs[IJK] = nhflow_face::V_s(p,d,i,j,k);
    d->Fn[IJK] = nhflow_face::V_n(p,d,i,j,k);
    }
    
    // flux y-dir
    VLOOP
    {
    d->Fe[IJK] = nhflow_face::V_e(p,d,i,j,k);
    d->Fw[IJK] = nhflow_face::V_w(p,d,i,j,k);
    }
    
    // flux z-dir
    nhflow_face::zflux(p,d,d->Vb,d->Vt);
    }
}

void nhflow_flux_build_f::start_W(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // flux x-dir
    ULOOP
    {
    d->Fs[IJK] = nhflow_face::W_s(p,d,i,j,k);
    d->Fn[IJK] = nhflow_face::W_n(p,d,i,j,k);
    }
    
    // flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
    d->Fe[IJK] = nhflow_face::W_e(p,d,i,j,k);
    d->Fw[IJK] = nhflow_face::W_w(p,d,i,j,k);
    }
    
    // flux z-dir
    nhflow_face::zflux(p,d,d->Wb,d->Wt);
}

void nhflow_flux_build_f::start_E(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
    // flux x-dir
    ULOOP
    {
    d->Fs[IJK] = d->UHs[IJK];
    
    d->Fn[IJK] = d->UHn[IJK];
    }
    
    // flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
    d->Fe[IJK] = d->VHe[IJK];
    
    d->Fw[IJK] = d->VHw[IJK];
    }
}

void nhflow_flux_build_f::start_U_yl(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
    // flux x-dir
    ULOOP
    {
    d->Fs[IJK] = d->UHs[IJK]*d->Us[IJK]/(PORVALNH1m*PORVALNH1m);
    
    d->Fn[IJK] = d->UHn[IJK]*d->Un[IJK]/(PORVALNH1*PORVALNH1);
    }
    
    // flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
    d->Fe[IJK] = d->Ve[IJK]*d->UHe[IJK]/(PORVALNH2m*PORVALNH2m);
    
    d->Fw[IJK] = d->Vw[IJK]*d->UHw[IJK]/(PORVALNH2*PORVALNH2);
    }
    
    // flux z-dir
    WLOOP
    d->Fz[IJK] = 0.5*(d->omegaF[FIJKp1]*(d->Ub[IJK] + d->Ut[IJK]))/(PORVALNH*PORVALNH) - 0.5*fabs(d->omegaF[FIJKp1])*(d->Ut[IJK] - d->Ub[IJK])/(PORVALNH*PORVALNH);
}

void nhflow_flux_build_f::start_V_yl(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->j_dir==1)
    {
    // flux x-dir
    ULOOP
    {
    d->Fs[IJK] = d->Us[IJK]*d->VHs[IJK]/(PORVALNH1m*PORVALNH1m);
    
    d->Fn[IJK] = d->Un[IJK]*d->VHn[IJK]/(PORVALNH1*PORVALNH1);
    }
    
    // flux y-dir
    VLOOP
    {
    d->Fe[IJK] = d->VHe[IJK]*d->Ve[IJK]/(PORVALNH2m*PORVALNH2m);
    
    d->Fw[IJK] = d->VHw[IJK]*d->Vw[IJK]/(PORVALNH2*PORVALNH2);
    }
    
    // flux z-dir
    WLOOP
    d->Fz[IJK] = 0.5*(d->omegaF[FIJKp1]*(d->Vb[IJK] + d->Vt[IJK]))/(PORVALNH*PORVALNH) - 0.5*fabs(d->omegaF[FIJKp1])*(d->Vt[IJK] - d->Vb[IJK])/(PORVALNH*PORVALNH);
    }
}
