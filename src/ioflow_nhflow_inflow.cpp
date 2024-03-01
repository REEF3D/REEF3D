/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void ioflow_f::inflow_nhflow(lexer *p, fdm_nhf *d,ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        U[Im1JK]=p->Ui;
        U[Im2JK]=p->Ui;
        U[Im3JK]=p->Ui;
		
        V[Im1JK]=0.0;
        V[Im2JK]=0.0;
        V[Im3JK]=0.0;
		
        W[Im1JK]=0.0;
        W[Im2JK]=0.0;
        W[Im3JK]=0.0;
        
        UH[Im1JK]=(d->eta(i,j)+d->depth(i,j))*p->Ui;
        UH[Im2JK]=(d->eta(i,j)+d->depth(i,j))*p->Ui;
        UH[Im3JK]=(d->eta(i,j)+d->depth(i,j))*p->Ui;
		
        VH[Im1JK]=0.0;
        VH[Im2JK]=0.0;
        VH[Im3JK]=0.0;
		
        WH[Im1JK]=0.0;
        WH[Im2JK]=0.0;
        WH[Im3JK]=0.0;
    }
}

void ioflow_f::rkinflow_nhflow(lexer *p, fdm_nhf *d,ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        U[Im3JK]=U[Im2JK]=U[Im1JK]=d->U[Im1JK];
        V[Im3JK]=V[Im2JK]=V[Im1JK]=d->V[Im1JK];
        W[Im3JK]=W[Im2JK]=W[Im1JK]=d->W[Im1JK];
        
        UH[Im3JK]=UH[Im2JK]=UH[Im1JK]=d->UH[Im1JK];
        VH[Im3JK]=VH[Im2JK]=VH[Im1JK]=d->VH[Im1JK];
        WH[Im3JK]=WH[Im2JK]=WH[Im1JK]=d->WH[Im1JK];
    }

}

void ioflow_f::wavegen_precalc_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    
}

void ioflow_f::wavegen_precalc_ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    
}