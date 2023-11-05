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

#include"nhflow_fsf_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"patchBC_interface.h"

nhflow_fsf_f::nhflow_fsf_f(lexer *p, fdm_nhf* d, ghostcell *pgc, ioflow *pflow, patchBC_interface *ppBC) : eps(1.0e-6),P(p),Q(p),K(p)
{
    pBC = ppBC;
    
    p->Iarray(temp,p->imax*p->jmax);
    
    if(p->F50==1)
	gcval_eta = 51;
    
    if(p->F50==2)
	gcval_eta = 52;
    
    if(p->F50==3)
	gcval_eta = 53;
    
    if(p->F50==4)
	gcval_eta = 54;
}

nhflow_fsf_f::~nhflow_fsf_f()
{
}

void nhflow_fsf_f::start(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow)
{
}

void nhflow_fsf_f::update(lexer *p, fdm_nhf* d, ghostcell *pgc, slice &f)
{
}

void nhflow_fsf_f::rk2_step1(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice &WLRK1, slice &WLRK2, double alpha)
{
    SLICELOOP4
    K(i,j) = 0.0;
    
    LOOP
    WETDRY
    K(i,j) += -p->DZN[KP]*((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP]  + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    SLICELOOP4
    WLRK1(i,j) = d->WL(i,j) + p->dt*K(i,j);
     
    pflow->WL_relax(p,pgc,WLRK1,d->depth);
    pflow->fsfinflow_nhflow(p,d,pgc,WLRK1);
    pgc->gcsl_start4(p,WLRK1,gcval_eta);
    
    SLICELOOP4
    {
    d->eta_n(i,j) = d->eta(i,j);
    d->detadt_n(i,j) = d->detadt(i,j);
    }
    
    SLICELOOP4
    d->eta(i,j) = WLRK1(i,j) - d->depth(i,j);
    
    SLICELOOP4
    d->detadt(i,j) = K(i,j);
    
    pgc->gcsl_start4(p,d->eta,gcval_eta);
    pgc->gcsl_start4(p,d->detadt,1);
    
    wetdry(p,d,pgc,U,V,W,WLRK1);
    //breaking(p,d,pgc,d->eta,d->eta_n,1.0);
}

void nhflow_fsf_f::rk2_step2(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice &WLRK1, slice &WLRK2, double alpha)
{
    SLICELOOP4
    K(i,j) = 0.0;
    
    LOOP
    WETDRY
    K(i,j) += -p->DZN[KP]*((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP]  + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    SLICELOOP4
    d->WL(i,j) = 0.5*d->WL(i,j) + 0.5*WLRK1(i,j) + 0.5*p->dt*K(i,j);

    pflow->WL_relax(p,pgc,d->WL,d->depth);
    pflow->fsfinflow_nhflow(p,d,pgc,d->WL);
    pgc->gcsl_start4(p,d->WL,gcval_eta);
    
    SLICELOOP4
    d->eta_n(i,j) = d->eta(i,j);
    
    SLICELOOP4
    d->eta(i,j) = d->WL(i,j) - d->depth(i,j);
    
    SLICELOOP4
    d->detadt_n(i,j) = d->detadt(i,j);
    
    SLICELOOP4
    d->detadt(i,j) = K(i,j);
    
    pgc->gcsl_start4(p,d->eta,gcval_eta);
    pgc->gcsl_start4(p,d->detadt,1);
    
    wetdry(p,d,pgc,U,V,W,d->WL);
    //breaking(p,d,pgc,d->eta,d->eta_n,1.0);
}

void nhflow_fsf_f::flux_update(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& etark1, slice &etark2, double alpha)
{    
    SLICELOOP4
    K(i,j) = 0.0;
    
    LOOP
    K(i,j) += - p->DZN[KP]*((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP]  + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir);
}


