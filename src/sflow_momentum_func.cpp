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

#include"sflow_momentum_func.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sflow_bcmom.h"
#include"sflow_reconstruct.h"
#include"sflow_convection.h"
#include"sflow_signal_speed.h"
#include"sflow_diffusion.h"
#include"sflow_pressure.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"sflow_fsf.h"

#define WLVL (fabs(WL(i,j))>(1.0*p->A544)?WL(i,j):1.0e20)

sflow_momentum_func::sflow_momentum_func(lexer *p, fdm2D *b, ghostcell *pgc)
                                                    : sflow_bcmom(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_uh=14;
	gcval_vh=15;
	gcval_wh=16;
   
}

sflow_momentum_func::~sflow_momentum_func()
{
}

void sflow_momentum_func::reconstruct(lexer *p, fdm2D *b, ghostcell *pgc, sflow_fsf *pfsf, sflow_signal_speed *pss, 
                                     sflow_reconstruct *precon, slice &WL, double *U, double *V, double *W, double *uh, double *vh, double *wh)
{
    starttime=pgc->timer();
    
    // reconstruct eta
    precon->reconstruct_x(p, pgc, b, b->eta, b->ETAs, b->ETAn);
    precon->reconstruct_y(p, pgc, b, b->eta, b->ETAe, b->ETAw);
    precon->reconstruct_WL(p, pgc, d);

    // reconstruct U 
    precon->reconstruct_x(p, pgc, b, U, b->Us, b->Un);
    precon->reconstruct_y(p, pgc, b, U, b->Ue, b->Uw);

    // reconstruct  V
    precon->reconstruct_x(p, pgc, b, V, b->Vs, b->Vn);
    precon->reconstruct_y(p, pgc, b, V, b->Ve, b->Vw);

    // reconstruct  W
    precon->reconstruct_x(p, pgc, b, W, b->Ws, b->Wn);
    precon->reconstruct_y(p, pgc, b, W, b->We, b->Ww);
    
    // reconstruct uh
    precon->reconstruct_x(p, pgc, b, uh, b->uhs, b->uhn);
    precon->reconstruct_y(p, pgc, b, uh, b->uhe, b->uhw);
    
    // reconstruct  vh
    precon->reconstruct_x(p, pgc, b, vh, b->vhs, b->vhn);
    precon->reconstruct_y(p, pgc, b, vh, b->vhe, b->vhw);
    
    // reconstruct  wh
    precon->reconstruct_x(p, pgc, b, wh, b->whs, b->whn);
    precon->reconstruct_y(p, pgc, b, wh, b->whe, b->whw);
    
    // wetdry
    pfsf->wetdry_fluxes(p,d,pgc,WL,U,V,W,uh,vh,wh);
    
    pss->signal_speed_update(p, pgc, b, b->Us, b->Un, b->Ve, b->Vw, b->Ds, b->Dn, b->De, b->Dw);
    
    p->recontime+=pgc->timer()-starttime;
}

void sflow_momentum_func::velcalc(lexer *p, fdm2D *b, ghostcell *pgc, double *uh, double *vh, double *wh, slice &WL)
{
    // Fr nuber limiter
    SLICELOOP4
    WETDRY
    {
    uh(i,j) = MIN(uh(i,j), p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    vh(i,j) = MIN(vh(i,j), p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    wh(i,j) = MIN(wh(i,j), p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));      
    
    uh(i,j) = MAX(uh(i,j), -p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    vh(i,j) = MAX(vh(i,j), -p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    wh(i,j) = MAX(wh(i,j), -p->A531*WL(i,j)*sqrt(9.81*WL(i,j))); 
    }
    
    
    if(p->B60==1)
    SLICELOOP4
    if(p->wet[Ip1J]==0 || p->wet[Im1J]==0 || p->wet[IJp1]==0 || p->wet[IJm1]==0 || p->deep[IJ]==0)
    {
    uh(i,j) = MIN(uh(i,j), 0.1*p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    vh(i,j) = MIN(vh(i,j), 0.1*p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    wh(i,j) = MIN(wh(i,j), 0.1*p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));      
    
    uh(i,j) = MAX(uh(i,j), -0.1*p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    vh(i,j) = MAX(vh(i,j), -0.1*p->A531*WL(i,j)*sqrt(9.81*WL(i,j)));
    wh(i,j) = MAX(wh(i,j), -0.1*p->A531*WL(i,j)*sqrt(9.81*WL(i,j))); 
    }
    
    
    SLICELOOP4
    WETDRY
    {
    b->U(i,j) = uh(i,j)/WLVL;
    b->V(i,j) = vh(i,j)/WLVL;
    b->W(i,j) = wh(i,j)/WLVL;       
    }
    

    SLICELOOP4
    if(p->wet[IJ]==0)
    {
    b->U(i,j) = 0.0;
    b->V(i,j) = 0.0;
    b->W(i,j) = 0.0;
    }
    
    pgc->gcsl_start4(p,UHRK1,gcval_u);
    pgc->gcsl_start4(p,VHRK1,gcval_u);
    pgc->gcsl_start4(p,WHRK1,gcval_u);
}

void sflow_momentum_func::irhs(lexer *p, fdm2D *b, ghostcell *pgc)
{
	SLICELOOP4
	{
	b->F(i,j) += (b->Fext(i,j))*PORVALNH;
	b->Fext(i,j)=0.0;
	}
}

void sflow_momentum_func::jrhs(lexer *p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
	{
	b->G(i,j) += (b->Gext(i,j))*PORVALNH;
	b->Gext(i,j)=0.0;
	}
}

void sflow_momentum_func::krhs(lexer *p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
	{
	b->H(i,j) += (b->Hext(i,j))*PORVALNH;
	b->Hext(i,j)=0.0;
	}
}

void sflow_momentum_func::clearrhs(lexer *p, fdm2D *b, ghostcell *pgc)
{
	n=0;
	SLICELOOP4
	{
	b->rhsvec.V[n]=0.0;
	++n;
	}
}

void sflow_momentum_func::inidisc(lexer *p, fdm2D *b, ghostcell *pgc, sflow_fsf *pfsf)
{
    /*
    pfsf->wetdry(p,d,pgc,b->U,b->V,b->W,b->WL);
    sigma_update(p,d,pgc,b->WL);
    pfsf->kinematic_fsf(p,d,b->U,b->V,b->W,b->eta);
    pfsf->kinematic_bed(p,d,b->U,b->V,b->W);
    velcalc(p,d,pgc,b->uh,b->vh,b->wh,b->WL);*/
}
     



