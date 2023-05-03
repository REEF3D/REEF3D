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

#include"iowave.h"
#include"lexer.h"
#include"ghostcell.h"
 
void iowave::inflow_nhflow(lexer *p, fdm_nhf *d, ghostcell* pgc, double *U, double *V, double *W)
{
    if(p->I230==0)
    {
    if(p->B98==0)
    nhflow_inflow_plain(p,d,pgc,U,V,W);
    
	if(p->B98==3)
	nhflow_dirichlet_wavegen(p,d,pgc,U,V,W);
	
	if(p->B98==4)
	nhflow_active_wavegen(p,d,pgc,U,V,W);
	}
    
	if(p->B99==3||p->B99==4||p->B99==5)
	nhflow_active_beach(p,d,pgc,U,V,W);
    
    //if(p->I230>0)
    //ff_inflow(p,d,pgc,U,V,W);
}

void iowave::rkinflow_nhflow(lexer *p, fdm_nhf *d, ghostcell* pgc, double *U, double *V, double *W)
{
}

void iowave::discharge_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{

}

void iowave::nhflow_inflow_plain(lexer *p, fdm_nhf *d, ghostcell* pgc, double *U, double *V, double *W)
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
    }
}
