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

void nhflow_fsf_f::rk3_step1(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& WLRK1, slice &WLRK2, double alpha)
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
    d->eta_n(i,j) = d->eta(i,j);
    
    SLICELOOP4
    d->eta(i,j) = WLRK1(i,j) - d->depth(i,j);
    
    SLICELOOP4
    d->detadt(i,j) = K(i,j);
    
    pgc->gcsl_start4(p,d->eta,gcval_eta);
    pgc->gcsl_start4(p,d->detadt,1);
    
    wetdry(p,d,pgc,U,V,W,WLRK1);
    //breaking(p,d,pgc,etark1,d->eta,1.0);
}

void nhflow_fsf_f::rk3_step2(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& WLRK1, slice &WLRK2, double alpha)
{
    SLICELOOP4
    K(i,j) = 0.0;
    
    LOOP
    WETDRY
    K(i,j) += -p->DZN[KP]*((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP]  + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    SLICELOOP4
    WLRK2(i,j) = 0.75*d->WL(i,j) + 0.25*WLRK1(i,j) + 0.25*p->dt*K(i,j);

    pflow->WL_relax(p,pgc,WLRK2,d->depth);
    pflow->fsfinflow_nhflow(p,d,pgc,WLRK2);
    pgc->gcsl_start4(p,WLRK2,gcval_eta);
    
    SLICELOOP4
    d->eta(i,j) = WLRK2(i,j) - d->depth(i,j);
    
    SLICELOOP4
    d->detadt(i,j) = K(i,j);
    
    pgc->gcsl_start4(p,d->eta,gcval_eta);
    pgc->gcsl_start4(p,d->detadt,1);
    
    wetdry(p,d,pgc,U,V,W,WLRK2);
    //breaking(p,d,pgc,etark2,etark1,0.25);
}

void nhflow_fsf_f::rk3_step3(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& WLRK1, slice &WLRK2, double alpha)
{
    SLICELOOP4
    K(i,j) = 0.0;
    
    LOOP
    WETDRY
    K(i,j) += - p->DZN[KP]*((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP]  + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    SLICELOOP4
    d->WL(i,j) = (1.0/3.0)*d->WL(i,j) + (2.0/3.0)*WLRK2(i,j) + (2.0/3.0)*p->dt*K(i,j);


    pflow->WL_relax(p,pgc,d->WL,d->depth);
    pflow->fsfinflow_nhflow(p,d,pgc,d->WL);
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
    //breaking(p,d,pgc,d->eta,etark2,2.0/3.0);
}


