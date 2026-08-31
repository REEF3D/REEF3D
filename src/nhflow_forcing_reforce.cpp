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
#include"nhflow_geometry.h"
#include"6DOF.h"
#include"nhflow_reinidisc_fsf.h"
#include"vrans.h"

void nhflow_forcing::reforcing(lexer *p, fdm_nhf *d, ghostcell *pgc, sixdof *p6dof, 
                             int iter, double alpha, double *UH, double *VH, double *WH, slice &WL, bool finalize)
{
    starttime=pgc->timer();
    
    // ini forcing terms
    reset(p,d,pgc);
    
    if(solid_flag==1)
    {
    // solid forcing
    solid_forcing(p,d,pgc,alpha,d->U,d->V,d->W,WL);
    }
    
    // 6DOF forcing
    p6dof->reforce_nhflow(p,d,pgc,iter,d->U,d->V,d->W,FX,FY,FZ,WL,fe,finalize);

    if(forcing_flag==1)
    {
    // add forcing term to RHS
    LOOP
    {
        UH[IJK]   += alpha*p->dt*CPORNH*FX[IJK]*WL(i,j);
        
        d->U[IJK] += alpha*p->dt*CPORNH*FX[IJK];
        
        d->maxF = MAX(fabs(d->maxF),alpha*p->dt*CPORNH*FX[IJK]);
    }
    
    LOOP
    {
        VH[IJK]   += alpha*p->dt*CPORNH*FY[IJK]*WL(i,j);
        
        d->V[IJK] += alpha*p->dt*CPORNH*FY[IJK];
        
        d->maxH = MAX(fabs(d->maxH),alpha*p->dt*CPORNH*FY[IJK]);
    }
    
    LOOP
    {
        WH[IJK]   += alpha*p->dt*CPORNH*FZ[IJK]*WL(i,j);
        
        d->W[IJK] += alpha*p->dt*CPORNH*FZ[IJK];
        
        d->maxG = MAX(fabs(d->maxG),alpha*p->dt*CPORNH*FZ[IJK]);
    }
    }

    // DLM
    if(dlm_flag==1)
    {
        
        LOOP
        {
            UH[IJK] += alpha*p->dt*CPORNH*FX[IJK]*WL(i,j);
            
            d->U[IJK] += alpha*p->dt*CPORNH*FX[IJK];
        }
        
        LOOP
        {
            VH[IJK] += alpha*p->dt*CPORNH*FY[IJK]*WL(i,j);
            
            d->V[IJK] += alpha*p->dt*CPORNH*FY[IJK];
        }
        
        LOOP
        {
            WH[IJK] += alpha*p->dt*CPORNH*FZ[IJK]*WL(i,j);
            
            d->W[IJK] += alpha*p->dt*CPORNH*FZ[IJK];
        }
    }
    
    pgc->gcsl_start4(p,d->eta,gcval_eta);
    pgc->gcsl_start4(p,WL,gcval_eta);
    pgc->gcsl_start4(p,d->bed,1);
    
    pgc->start4V(p,d->U,gcval_u);
    pgc->start4V(p,d->V,gcval_v);
    pgc->start4V(p,d->W,gcval_w);
    
    pgc->start4V(p,UH,gcval_uh);
    pgc->start4V(p,VH,gcval_vh);
    pgc->start4V(p,WH,gcval_wh);
    
    pgc->gciobc_update(p,d);
    
    p->dftime+=pgc->timer()-starttime;
}

