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

#include"sflow_momentum_RK2.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sflow_fsf.h"

sflow_momentum_RK2::sflow_momentum_RK2(lexer *p, fdm2D *b, ghostcell *pgc, sflow_HLL *pphll, sflow_signal_speed *ppss, 
                                       sflow_reconstruct *pprecon, sflow_diffusion *ppdiff, sflow_pressure *pppress, 
                                       solver2D *ppsolv, solver2D *pppoissonsolv, ioflow *ppflow, sflow_fsf *ppfsf, 
                                       sflow_forcing *ppsfdf, sixdof *pp6dof)
                                      : sflow_momentum_func(p,b,pgc,pphll,ppss,pprecon,ppdiff,pppress,ppsolv,pppoissonsolv,ppflow,ppfsf,ppsfdf,pp6dof),
                                        WLRK1(p),UHRK1(p),VHRK1(p),WHRK1(p)
{
}

sflow_momentum_RK2::~sflow_momentum_RK2()
{
}

void sflow_momentum_RK2::start(lexer *p, fdm2D* b, ghostcell* pgc)
{
    // in- and outflow ghost cells
    inflow(p,b,pgc,pflow);
    
//Step 1
//--------------------------------------------------------
    stage(p,b,pgc, b->WL,b->UH,b->VH,b->WH, WLRK1,UHRK1,VHRK1,WHRK1, 0.0, 0, 0);
    
//Step 2
//--------------------------------------------------------
    stage(p,b,pgc, WLRK1,UHRK1,VHRK1,WHRK1, b->WL,b->UH,b->VH,b->WH, 0.5, 1, 1);
    
    pfsf->breaking_persist(p,b,pgc,b->eta,b->eta_n,1.0);
    
    SLICELOOP4
    b->eta_n(i,j) = b->eta(i,j);
    
    pgc->gcsl_start4(p,b->eta_n,gcval_eta);
}
