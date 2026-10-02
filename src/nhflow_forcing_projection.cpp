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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"6DOF.h"
#include"nhflow_pressure.h"
#include"nhflow_momentum_func.h"

void nhflow_forcing::projection(lexer *p, fdm_nhf *d, ghostcell *pgc, sixdof *p6dof, nhflow_pressure *ppress, solver *psolv,
                                ioflow *pflow, nhflow_momentum_func *pmom, int iter, double alpha,
                                double *UH, double *VH, double *WH, slice &WL)
{
    // strongly coupled flexible membranes (X 330, membrane.dat 'coupling iterated'), called by
    // nhflow_momentum_func::phase_P in place of ppress->start: the projection is repeated until the membrane velocities
    // the forcing used and the ones the structure returns for the resulting loads agree (net_membrane_coupling.cpp).
    // Each iteration starts from the predicted, forced state of the stage and the pressure of the stage start, so the
    // fluid response is a deterministic function of the membrane velocity; the membrane module adds the forcing
    // response to the change of its velocities. The pressure solver is called unchanged. phase_P2 (velocities,
    // reforcing) follows as without membranes.
    const int n = p->imax*p->jmax*(p->kmax+2);
    
    if(CUH==nullptr)
    {
        p->Darray(CUH,n);
        p->Darray(CVH,n);
        p->Darray(CWH,n);
        p->Darray(CU,n);
        p->Darray(CV,n);
        p->Darray(CW,n);
        p->Darray(CP,n);
    }
    
    for(int q=0; q<n; ++q)
    {
        CUH[q] = UH[q];
        CVH[q] = VH[q];
        CWH[q] = WH[q];
        CU[q] = d->U[q];
        CV[q] = d->V[q];
        CW[q] = d->W[q];
        CP[q] = d->P[q];
    }
    
    for(int it=0; ; ++it)
    {
        if(it>0)
        {
            for(int q=0; q<n; ++q)
            {
                UH[q] = CUH[q];
                VH[q] = CVH[q];
                WH[q] = CWH[q];
                d->U[q] = CU[q];
                d->V[q] = CV[q];
                d->W[q] = CW[q];
                d->P[q] = CP[q];
            }
            
            p6dof->membrane_reforce_nhflow(p,d,pgc,alpha,UH,VH,WH,WL);
            
            pgc->start4V(p,d->U,gcval_u);
            pgc->start4V(p,d->V,gcval_v);
            pgc->start4V(p,d->W,gcval_w);
            
            pgc->start4V(p,UH,gcval_uh);
            pgc->start4V(p,VH,gcval_vh);
            pgc->start4V(p,WH,gcval_wh);
            
            pgc->gciobc_update(p,d);
        }
        
        ppress->start(p,d,psolv,pgc,pflow,WL,UH,VH,WH,alpha);
        
        // velocities for the membrane loads
        pmom->velcalc(p,d,pgc,UH,VH,WH,WL,alpha);
        
        if(p6dof->membrane_couple_nhflow(p,d,pgc,iter,alpha,WL,it))
        break;
    }
    
    // phase_P2 computes the velocities again (DWDT from the W before the projection)
    for(int q=0; q<n; ++q)
    d->W[q] = CW[q];
}
