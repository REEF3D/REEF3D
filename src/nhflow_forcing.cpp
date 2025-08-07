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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"6DOF.h"
#include"nhflow_reinidisc_fsf.h"

nhflow_forcing::nhflow_forcing(lexer *p, fdm_nhf *d, ghostcell *pgc) : epsi(1.6), fe(p)
{
    forcing_flag=0;
    solid_flag=0;
    floating_flag=0;
    dlm_flag=0;
    
    if(p->A581>0 || p->A583>0 || p->A584>0   || p->A585>0  || p->A586>0 || p->A587>0 || p->A588>0 || p->A589>0 || p->A590>0)
    {
    forcing_flag=1;
    solid_flag=1;
    }
    
    if(p->X10>0)
    {
    forcing_flag=1;
    floating_flag=1;
    }
    
    if(p->A599==1)
    {
    dlm_flag=1;
    forcing_flag=0;
    solid_flag=0;
    floating_flag=0;
    }
    
    // ----
    if(forcing_flag==1)
    {
    p->Iarray(IO,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(CL,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(CR,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(FRK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(dt,p->imax*p->jmax*(p->kmax+2));
    p->Darray(L,p->imax*p->jmax*(p->kmax+2));

    prdisc = new nhflow_reinidisc_fsf(p);
    }
    
    p->Darray(FX,p->imax*p->jmax*(p->kmax+2));
    p->Darray(FY,p->imax*p->jmax*(p->kmax+2));
    p->Darray(FZ,p->imax*p->jmax*(p->kmax+2));
    
    if(dlm_flag==1)
    dlm_forcing_ini(p,pgc);
    
    if(p->F50==1)
	gcval_eta = 51;
    
    if(p->F50==2)
	gcval_eta = 52;
    
    if(p->F50==3)
	gcval_eta = 53;
    
    if(p->F50==4)
	gcval_eta = 54;
    
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_uh=14;
	gcval_vh=15;
	gcval_wh=16;
}

nhflow_forcing::~nhflow_forcing()
{
}

void nhflow_forcing::forcing(lexer *p, fdm_nhf *d, ghostcell *pgc, sixdof *p6dof, 
                             int iter, double alpha, double *UH, double *VH, double *WH, slice &WL, bool finalize)
{
    starttime=pgc->timer();
    
    // ini forcing terms
    reset(p,d,pgc);
    
    if(solid_flag==1)
    {
    // update direct forcing function
    ray_cast(p, d, pgc);
    reini_RK2(p, d, pgc, d->SOLID);
    
    // solid forcing
    solid_forcing(p,d,pgc,alpha,d->U,d->V,d->W,WL);
    }
    
    // 6DOF forcing
    p6dof->start_nhflow(p,d,pgc,iter,d->U,d->V,d->W,FX,FY,FZ,WL,fe,finalize);


    if(forcing_flag==1)
    {
    // add forcing term to RHS
    LOOP
    {
        UH[IJK]   += alpha*p->dt*CPORNH*FX[IJK]*WL(i,j);
        
        d->U[IJK] += alpha*p->dt*CPORNH*FX[IJK];
    }
    
    LOOP
    {
        VH[IJK]   += alpha*p->dt*CPORNH*FY[IJK]*WL(i,j);
        
        d->V[IJK] += alpha*p->dt*CPORNH*FY[IJK];
    }
    
    LOOP
    {
        WH[IJK]   += alpha*p->dt*CPORNH*FZ[IJK]*WL(i,j);
        
        d->W[IJK] += alpha*p->dt*CPORNH*FZ[IJK];
    }
    
    
    // DF
    LOOP
    p->DF[IJK]=1;
    
    if(solid_flag==1)
    LOOP
    if(d->SOLID[IJK]<0.0)
    p->DF[IJK]=-1;

    if(floating_flag==1)
    LOOP
    if(d->FB[IJK]<0.0)
    p->DF[IJK]=-1;
    
    pgc->startintV(p,p->DF,1);
    
    // DFSL slice
    pgc->gcsldf_update(p);
    pgc->solid_forcing_eta(p,WL);
    pgc->solid_forcing_eta(p,d->eta);
    pgc->solid_forcing_bed(p,d->bed);
    }

    // DLM
    if(dlm_flag==1)
    {
        dlm_forcecalc(p,d,pgc,alpha,d->U,d->V,d->W,WL);
        dlm_forcing(p,d,pgc,alpha,d->U,d->V,d->W,WL);
        
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
