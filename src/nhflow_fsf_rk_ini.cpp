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

#include"nhflow_fsf_rk.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"sflow_hxy_weno.h"
#include"nhflow_flux_HLL.h"

void nhflow_fsf_rk::ini(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W)
{    
    wetdry(p,d,pgc,U,V,W,d->eta);
    
    pfluxfsf->face_flux_3D(p,pgc,d,d->eta,U,V,Fx,Fy);
    
    SLICELOOP1
    P(i,j)=0.0;
    
    SLICELOOP2
    Q(i,j)=0.0;

    ULOOP
    P(i,j) += 0.5*(U[IJK] + U[Ip1JK]);

    VLOOP
	Q(i,j) += 0.5*(V[IJK] + V[IJp1K]);
    
    pgc->gcsl_start1(p,P,10);
    pgc->gcsl_start2(p,Q,11);
    
    phxy->start(p,d->hx,d->hy,d->depth,p->wet,d->eta,P,Q);
    
    pgc->gcsl_start1(p,d->hx,10);
    pgc->gcsl_start2(p,d->hy,11);
    
    SLICELOOP4
    d->detadt(i,j) = 0.0;
    
    pgc->gcsl_start4(p,d->detadt,1);
    pgc->start1V(p,Fx,10);
    pgc->start2V(p,Fy,10);
    
    LOOP
    d->detadt(i,j) += -p->DZN[KP]*((Fx[IJK] - Fx[Im1JK])/p->DXN[IP]  + (Fy[IJK] - Fy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    pgc->gcsl_start4(p,d->detadt,1);
    
    pgc->start4V(p,Fx,10);
    
    //LOOP    
    //d->test[IJK] = Fx[IJK] - Fx[Im1JK];
    
    pgc->start4V(p,d->test,10);
     
    /*
    breaking(p,d,pgc,d->eta,d->eta,1.0);
    p->sigma_update(p,d,pgc,d->eta,d->eta,1.0);
    p->omega_update(p,d,pgc,U,V,W,d->eta,d->eta,1.0);*/
}