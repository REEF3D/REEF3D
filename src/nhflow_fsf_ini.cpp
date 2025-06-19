/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

void nhflow_fsf_f::ini(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W)
{   
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    wetdry(p,d,pgc,U,V,W,d->WL);
    
    SLICELOOP4
    d->detadt(i,j) = 0.0;
    
    pgc->gcsl_start4(p,d->detadt,1);
    
    LOOP
    d->detadt(i,j) += -p->DZN[KP]*((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP]  + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    pgc->gcsl_start4(p,d->detadt,1);
    
    pgc->start1V(p,d->Fx,10);
    
    pgc->start4V(p,d->test,1);
    
    // fsf guard
    if(p->A580==1)
    {
    guard_is = p->posc_i(p->A580_xs);
    guard_ie = p->posc_i(p->A580_xe);
    
    guard_js = p->posc_j(p->A580_ys);
    guard_je = p->posc_j(p->A580_ye);
    }
    
    if(p->A580==1)
    SLICELOOP4
    if(i>=guard_is && i<=guard_ie && j>=guard_js && j<=guard_je) 
    {
    
    p->flagfsf[IJ]=0;
        
    }
    
    pgc->gcslflagx(p,p->flagfsf);
    
    // FSF Box
    if(p->F72>0)
    {
    int istart, iend, jstart, jend, kstart, kend;
    
        for(int qn=0;qn<p->F72;++qn)
        {
            istart = p->posc_i(p->F72_xs[qn]);
            iend = p->posc_i(p->F72_xe[qn]);
            
            jstart = p->posc_j(p->F72_ys[qn]);
            jend = p->posc_j(p->F72_ye[qn]);

            SLICELOOP4
            if(i>=istart && i<iend && j>=jstart && j<jend)
            d->eta(i,j)= p->F72_h[qn] - p->F60;

        }
    }
    
    
    wetdry(p,d,pgc,U,V,W,d->WL);
    
    SLICELOOP4
    d->WL(i,j) = MAX(p->A544,d->eta(i,j) + d->depth(i,j));
    
    SLICELOOP4
    d->eta_n(i,j) = d->eta(i,j);
    
    pgc->gcsl_start4(p,d->eta,50);
    pgc->gcsl_start4(p,d->WL,50);
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    
    pgc->start1V(p,d->Fx,10);
    pgc->start2V(p,d->Fy,10);
}