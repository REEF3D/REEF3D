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

#include"nhflow_momentum_func.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_bcmom.h"
#include"nhflow_reconstruct.h"
#include"nhflow_convection.h"
#include"nhflow_signal_speed.h"
#include"nhflow_diffusion.h"
#include"nhflow_pressure.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"nhflow_fsf.h"

#define WLVL (fabs(WL(i,j))>(1.0*p->A544)?WL(i,j):1.0e20)

nhflow_momentum_func::nhflow_momentum_func(lexer *p, fdm_nhf *d, ghostcell *pgc)
                                                    : nhflow_bcmom(p), nhflow_sigma(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_uh=14;
	gcval_vh=15;
	gcval_wh=16;
   
}

nhflow_momentum_func::~nhflow_momentum_func()
{
}

void nhflow_momentum_func::reconstruct(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_fsf *pfsf, nhflow_signal_speed *pss, 
                                     nhflow_reconstruct *precon, slice &WL, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    starttime=pgc->timer();
    
    // reconstruct eta
    precon->reconstruct_2D_x(p, pgc, d, d->eta, d->ETAs, d->ETAn);
    precon->reconstruct_2D_y(p, pgc, d, d->eta, d->ETAe, d->ETAw);
    precon->reconstruct_2D_WL(p, pgc, d);

    // reconstruct U 
    precon->reconstruct_3D_x(p, pgc, d, U, d->Us, d->Un);
    precon->reconstruct_3D_y(p, pgc, d, U, d->Ue, d->Uw);
    precon->reconstruct_3D_z(p, pgc, d, U, d->Ub, d->Ut);
    
    // reconstruct  V
    precon->reconstruct_3D_x(p, pgc, d, V, d->Vs, d->Vn);
    precon->reconstruct_3D_y(p, pgc, d, V, d->Ve, d->Vw);
    precon->reconstruct_3D_z(p, pgc, d, V, d->Vb, d->Vt);
    
    // reconstruct  W
    precon->reconstruct_3D_x(p, pgc, d, W, d->Ws, d->Wn);
    precon->reconstruct_3D_y(p, pgc, d, W, d->We, d->Ww);
    precon->reconstruct_3D_z(p, pgc, d, W, d->Wb, d->Wt);
    
    // reconstruct UH
    precon->reconstruct_3D_x(p, pgc, d, UH, d->UHs, d->UHn);
    precon->reconstruct_3D_y(p, pgc, d, UH, d->UHe, d->UHw);
    
    // reconstruct  VH
    precon->reconstruct_3D_x(p, pgc, d, VH, d->VHs, d->VHn);
    precon->reconstruct_3D_y(p, pgc, d, VH, d->VHe, d->VHw);
    
    // reconstruct  WH
    precon->reconstruct_3D_x(p, pgc, d, WH, d->WHs, d->WHn);
    precon->reconstruct_3D_y(p, pgc, d, WH, d->WHe, d->WHw);
    
    // wetdry
    pfsf->wetdry_fluxes(p,d,pgc,WL,U,V,W,UH,VH,WH);
    
    pss->signal_speed_update(p, pgc, d, d->Us, d->Un, d->Ve, d->Vw, d->Ds, d->Dn, d->De, d->Dw);
    
    p->recontime+=pgc->timer()-starttime;
}

void nhflow_momentum_func::velcalc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *UH, double *VH, double *WH, slice &WL)
{
    // Fr nuber limiter
    LOOP
    WETDRY
    {
    UH[IJK] = MIN(UH[IJK], p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    VH[IJK] = MIN(VH[IJK], p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    WH[IJK] = MIN(WH[IJK], p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));      
    
    UH[IJK] = MAX(UH[IJK], -p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    VH[IJK] = MAX(VH[IJK], -p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    WH[IJK] = MAX(WH[IJK], -p->A531*WL(i,j)*sqrt(9.81*WL(i,j))); 
    }
    
    
    if(p->B60==1)
    LOOP
    if(p->wet[Ip1J]==0 || p->wet[Im1J]==0 || p->wet[IJp1]==0 || p->wet[IJm1]==0 || p->deep[IJ]==0)
    {
    UH[IJK] = MIN(UH[IJK], 0.1*p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    VH[IJK] = MIN(VH[IJK], 0.1*p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    WH[IJK] = MIN(WH[IJK], 0.1*p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));      
    
    UH[IJK] = MAX(UH[IJK], -0.1*p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    VH[IJK] = MAX(VH[IJK], -0.1*p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    WH[IJK] = MAX(WH[IJK], -0.1*p->A531*WL(i,j)*sqrt(9.81*WL(i,j))); 
    }
    
    
    LOOP
    WETDRY
    {
    d->U[IJK] = UH[IJK]/WLVL;
    d->V[IJK] = VH[IJK]/WLVL;
    d->W[IJK] = WH[IJK]/WLVL;       
    }
    

    LOOP
    if(p->wet[IJ]==0)
    {
    d->U[IJK] = 0.0;
    d->V[IJK] = 0.0;
    d->W[IJK] = 0.0;
    }
    
    pgc->start4V(p,d->U,gcval_u);
    pgc->start4V(p,d->V,gcval_v);
    pgc->start4V(p,d->W,gcval_w);
}

void nhflow_momentum_func::irhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    /*
	n=0;
	LOOP
	{
    d->maxF=MAX(fabs(d->rhsvec.V[n]),d->maxF);
	d->F[IJK] += (d->rhsvec.V[n])*PORVALNH;
	d->rhsvec.V[n]=0.0;
	++n;
	}*/
}

void nhflow_momentum_func::jrhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    /*
	n=0;
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n]),a->maxG);
	a->G[IJK] += (a->rhsvec.V[n])*PORVAL2;
	a->rhsvec.V[n]=0.0;
	++n;
	}*/
}

void nhflow_momentum_func::krhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
}

void nhflow_momentum_func::clearrhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
	n=0;
	LOOP
	{
	d->rhsvec.V[n]=0.0;
	++n;
	}
}

void nhflow_momentum_func::inidisc(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_fsf *pfsf)
{
    pfsf->wetdry(p,d,pgc,d->U,d->V,d->W,d->WL);
    sigma_update(p,d,pgc,d->WL);
    pfsf->kinematic_fsf(p,d,d->U,d->V,d->W,d->eta);
    pfsf->kinematic_bed(p,d,d->U,d->V,d->W);
    velcalc(p,d,pgc,d->UH,d->VH,d->WH,d->WL);
}
     



