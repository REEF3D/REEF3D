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

#include"nhflow_flux_build_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"slice.h"
#include"patchBC_interface.h"

#define WLVL (fabs(d->WL_n(i,j))>1.0e-20?d->WL_n(i,j):1.0e20)

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
            + p->A523*(0.5*fabs(p->W22)*d->ETAs(i,j)*d->ETAs(i,j) + fabs(p->W22)*d->ETAs(i,j)*0.5*(d->depth(i,j) + d->depth(i+1,j)))
            + (1.0-p->A523)*(0.5*fabs(p->W22)*d->ETAs_n(i,j)*d->ETAs_n(i,j) + fabs(p->W22)*d->ETAs_n(i,j)*0.5*(d->depth(i,j) + d->depth(i+1,j)));
    
    d->Fn[IJK] = d->UHn[IJK]*d->Un[IJK]
            + p->A523*(0.5*fabs(p->W22)*d->ETAn(i,j)*d->ETAn(i,j) + fabs(p->W22)*d->ETAn(i,j)*0.5*(d->depth(i,j) + d->depth(i+1,j)))
            + (1.0-p->A523)*(0.5*fabs(p->W22)*d->ETAn_n(i,j)*d->ETAn_n(i,j) + fabs(p->W22)*d->ETAn_n(i,j)*0.5*(d->depth(i,j) + d->depth(i+1,j)));
    }
    
    // flux y-dir
    VLOOP
    {
    d->Fe[IJK] = d->UHe[IJK]*d->Ve[IJK];
    
    d->Fw[IJK] = d->UHw[IJK]*d->Vw[IJK];
    }
    
    // flux z-dir
    WLOOP
    d->Fz[IJK] = 0.5*(d->omegaF[FIJKp1]*(d->Ub[IJK] + d->Ut[IJK]) - fabs(d->omegaF[FIJKp1])*(d->Ub[IJK] - d->Ut[IJK]));
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
            + 0.5*fabs(p->W22)*d->ETAe(i,j)*d->ETAe(i,j) + fabs(p->W22)*d->ETAe(i,j)*0.5*(d->depth(i,j) + d->depth(i,j+1));
    
    d->Fw[IJK] = d->VHw[IJK]*d->Vw[IJK] 
            + 0.5*fabs(p->W22)*d->ETAw(i,j)*d->ETAw(i,j) + fabs(p->W22)*d->ETAw(i,j)*0.5*(d->depth(i,j) + d->depth(i,j+1));
    }
    
    // flux z-dir
    WLOOP
    d->Fz[IJK] = 0.5*(d->omegaF[FIJKp1]*(d->Vb[IJK] + d->Vt[IJK]) - fabs(d->omegaF[FIJKp1])*(d->Vb[IJK] - d->Vt[IJK]));
    }
}

void nhflow_flux_build_f::start_W(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // flux x-dir
    ULOOP
    {
    d->Fs[IJK] = d->Ds(i,j)*d->Us[IJK]*d->Ws[IJK];
    
    d->Fn[IJK] = d->Dn(i,j)*d->Un[IJK]*d->Wn[IJK];
    }
    
    // flux y-dir
    VLOOP
    {
    d->Fe[IJK] = d->De(i,j)*d->Ve[IJK]*d->We[IJK];
    
    d->Fw[IJK] = d->Dw(i,j)*d->Vw[IJK]*d->Ww[IJK];
    }
    
    // flux z-dir
    WLOOP
    d->Fz[IJK] = 0.5*(d->omegaF[FIJKp1]*(d->Wb[IJK] + d->Wt[IJK]) - fabs(d->omegaF[FIJKp1])*(d->Wb[IJK] - d->Wt[IJK]));

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
    VLOOP
    {
    d->Fe[IJK] = d->VHe[IJK];
    
    d->Fw[IJK] = d->VHw[IJK];
    }
    
}





