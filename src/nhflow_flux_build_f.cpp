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

#include"nhflow_flux_build_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"slice.h"
#include"patchBC_interface.h"

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
    d->Fs[IJK] = d->UHs[IJK]*d->Us[IJK]
            + 0.5*fabs(p->W22)*d->ETAs(i,j)*d->ETAs(i,j)
            + fabs(p->W22)*d->ETAs(i,j)*d->dfx(i,j);
    
    d->Fn[IJK] = d->UHn[IJK]*d->Un[IJK]
            + 0.5*fabs(p->W22)*d->ETAn(i,j)*d->ETAn(i,j)
            + fabs(p->W22)*d->ETAn(i,j)*d->dfx(i,j);
    }
    
    // flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
    d->Fe[IJK] = d->Ve[IJK]*d->UHe[IJK];
    
    d->Fw[IJK] = d->Vw[IJK]*d->UHw[IJK];
    }
    
    // flux z-dir
    WLOOP
    d->Fz[IJK] = 0.5*(d->omegaF[FIJKp1]*(d->Ub[IJK] + d->Ut[IJK])) - 0.5*fabs(d->omegaF[FIJKp1])*(d->Ut[IJK] - d->Ub[IJK]);
}

void nhflow_flux_build_f::start_V(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->j_dir==1)
    {
    // flux x-dir
    ULOOP
    {
    d->Fs[IJK] = d->Us[IJK]*d->VHs[IJK];
    
    d->Fn[IJK] = d->Un[IJK]*d->VHn[IJK];
    }
    
    // flux y-dir
    VLOOP
    {
    d->Fe[IJK] = d->VHe[IJK]*d->Ve[IJK] 
            + 0.5*fabs(p->W22)*d->ETAe(i,j)*d->ETAe(i,j)
            + fabs(p->W22)*d->ETAe(i,j)*d->dfy(i,j);
    
    d->Fw[IJK] = d->VHw[IJK]*d->Vw[IJK] 
            + 0.5*fabs(p->W22)*d->ETAw(i,j)*d->ETAw(i,j) 
            + fabs(p->W22)*d->ETAw(i,j)*d->dfy(i,j);
    }
    
    // flux z-dir
    WLOOP
    d->Fz[IJK] = 0.5*(d->omegaF[FIJKp1]*(d->Vb[IJK] + d->Vt[IJK])) - 0.5*fabs(d->omegaF[FIJKp1])*(d->Vt[IJK] - d->Vb[IJK]);
    }
}

void nhflow_flux_build_f::start_W(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // flux x-dir
    ULOOP
    {
    d->Fs[IJK] = d->Us[IJK]*d->WHs[IJK];
    
    d->Fn[IJK] = d->Un[IJK]*d->WHn[IJK];
    }
    
    // flux y-dir
    if(p->j_dir==1)
    VLOOP
    {
    d->Fe[IJK] = d->Ve[IJK]*d->WHe[IJK];
    
    d->Fw[IJK] = d->Vw[IJK]*d->WHw[IJK];
    }
    
    // flux z-dir
    WLOOP
    d->Fz[IJK] = 0.5*(d->omegaF[FIJKp1]*(d->Wb[IJK] + d->Wt[IJK])) - 0.5*fabs(d->omegaF[FIJKp1])*(d->Wt[IJK] - d->Wb[IJK]);
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





