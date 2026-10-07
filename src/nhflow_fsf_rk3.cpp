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

#include"nhflow_fsf_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"seastate_nhflow.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"patchBC_interface.h"
#include"vrans.h"

void nhflow_fsf_f::rk3_step1(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& WLRK1, slice &WLRK2, double alpha)
{
    pflow->fsfinflow_flux_nhflow(p,d,pgc);
    
    SLICELOOP4
    K(i,j) = 0.0;
    
    LOOP
    WETDRY
    K(i,j) += -p->DZN[KP]*((d->FEx[IJK] - d->FEx[Im1JK])/p->DXN[IP]  + (d->FEy[IJK] - d->FEy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    // VRANS: d/dt int_{-h}^{eta} n dz = n(eta) deta/dt -> divide the column divergence by the
    // porosity of the surface cell, not layer by layer (mass error for vertically varying n)
    SLICELOOP4
    {
    k = p->knoz-1;
    K(i,j) /= PORVALNH;
    }
    
    // REEF3D::SEASTATE: Stokes transport in the continuity (A 751 2)
    if(d->wave!=nullptr)
    d->wave->mass_source(p,d,K);
    
    SLICELOOP4
    WLRK1(i,j) = d->WL(i,j) + p->dt*K(i,j);
     
    pflow->WL_relax(p,pgc,WLRK1,d->depth);
    pflow->fsfinflow_nhflow(p,d,pgc,WLRK1);
    if(d->wave!=nullptr)
    d->wave->wl_ghostcells(p,d,WLRK1);
    pgc->gcsl_start4(p,WLRK1,gcval_eta);
    
    SLICELOOP4
    d->eta_n(i,j) = d->eta(i,j);
    
    SLICELOOP4
    d->eta(i,j) = WLRK1(i,j) - d->depth(i,j);
    
    SLICELOOP4
    d->detadt(i,j) = K(i,j);
    
    pgc->gcsl_start4(p,d->eta,gcval_eta);
    pgc->gcsl_start4(p,d->detadt,1);
    
    wetdry(p,d,pgc,U,V,W,WLRK1);
}

void nhflow_fsf_f::rk3_step2(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& WLRK1, slice &WLRK2, double alpha)
{
    pflow->fsfinflow_flux_nhflow(p,d,pgc);
    
    SLICELOOP4
    K(i,j) = 0.0;
    
    LOOP
    WETDRY
    K(i,j) += -p->DZN[KP]*((d->FEx[IJK] - d->FEx[Im1JK])/p->DXN[IP]  + (d->FEy[IJK] - d->FEy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    // VRANS: d/dt int_{-h}^{eta} n dz = n(eta) deta/dt -> divide the column divergence by the
    // porosity of the surface cell, not layer by layer (mass error for vertically varying n)
    SLICELOOP4
    {
    k = p->knoz-1;
    K(i,j) /= PORVALNH;
    }
    
    // REEF3D::SEASTATE: Stokes transport in the continuity (A 751 2)
    if(d->wave!=nullptr)
    d->wave->mass_source(p,d,K);
    
    SLICELOOP4
    WLRK2(i,j) = 0.75*d->WL(i,j) + 0.25*WLRK1(i,j) + 0.25*p->dt*K(i,j);

    pflow->WL_relax(p,pgc,WLRK2,d->depth);
    pflow->fsfinflow_nhflow(p,d,pgc,WLRK2);
    if(d->wave!=nullptr)
    d->wave->wl_ghostcells(p,d,WLRK2);
    pgc->gcsl_start4(p,WLRK2,gcval_eta);
    
    SLICELOOP4
    d->eta(i,j) = WLRK2(i,j) - d->depth(i,j);
    
    SLICELOOP4
    d->detadt(i,j) = K(i,j);
    
    pgc->gcsl_start4(p,d->eta,gcval_eta);
    pgc->gcsl_start4(p,d->detadt,1);
    
    wetdry(p,d,pgc,U,V,W,WLRK2);
}

void nhflow_fsf_f::rk3_step3(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& WLRK1, slice &WLRK2, double alpha)
{
    pflow->fsfinflow_flux_nhflow(p,d,pgc);
    
    SLICELOOP4
    K(i,j) = 0.0;
    
    LOOP
    WETDRY
    K(i,j) += - p->DZN[KP]*((d->FEx[IJK] - d->FEx[Im1JK])/p->DXN[IP]  + (d->FEy[IJK] - d->FEy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    // VRANS: d/dt int_{-h}^{eta} n dz = n(eta) deta/dt -> divide the column divergence by the
    // porosity of the surface cell, not layer by layer (mass error for vertically varying n)
    SLICELOOP4
    {
    k = p->knoz-1;
    K(i,j) /= PORVALNH;
    }
    
    // REEF3D::SEASTATE: Stokes transport in the continuity (A 751 2)
    if(d->wave!=nullptr)
    d->wave->mass_source(p,d,K);
    
    SLICELOOP4
    d->WL(i,j) = (1.0/3.0)*d->WL(i,j) + (2.0/3.0)*WLRK2(i,j) + (2.0/3.0)*p->dt*K(i,j);


    pflow->WL_relax(p,pgc,d->WL,d->depth);
    pflow->fsfinflow_nhflow(p,d,pgc,d->WL);
    if(d->wave!=nullptr)
    d->wave->wl_ghostcells(p,d,d->WL);
    pgc->gcsl_start4(p,d->WL,gcval_eta);
    
    //SLICELOOP4
    //d->eta_n(i,j) = d->eta(i,j);

    SLICELOOP4
    d->eta(i,j) = d->WL(i,j) - d->depth(i,j);
    
    SLICELOOP4
    d->detadt(i,j) = K(i,j);
    
    pgc->gcsl_start4(p,d->eta,gcval_eta);
    pgc->gcsl_start4(p,d->detadt,1);

    wetdry(p,d,pgc,U,V,W,d->WL);
}


