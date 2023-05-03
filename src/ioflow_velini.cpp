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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void ioflow_f::velini(lexer *p, fdm *a, ghostcell *pgc)
{
    double Ai,area;

    area=0.0;
    Ai=0.0;
    p->Qi=0.0;
    p->Qo=0.0;

    // in
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        if(a->phi(i-1,j,k)>=0.0)
        {

            if(a->phi(i-1,j,k)>=0.5*p->DXM)
            area=p->DXM*p->DXM;

            if(a->phi(i-1,j,k)<0.5*p->DXM)
            area=p->DXM*(p->DXM*0.5 + a->phi(i-1,j,k));

            Ai+=area;
            p->Qi+=area*a->u(i-1,j,k);
        }
    }

    if(p->B60==1)
    {
    p->Ua=p->Ui=p->W10/(Ai>1.0e-20?Ai:1.0e20); 
    p->Qi=p->W10;
    }
    
    if(p->B60==2||p->B60==4)
    {
    p->Ua=p->Ui=hydrograph_ipol(p,pgc,hydro_in,hydro_in_count)/(Ai>1.0e-20?Ai:1.0e20); 
    p->Qi=hydrograph_ipol(p,pgc,hydro_in,hydro_in_count);
    }
	
	if(p->B60==3||p->B60==4)
    {
    p->Uo=p->Ui=hydrograph_ipol(p,pgc,hydro_out,hydro_out_count)/(Ao>1.0e-20?Ao:1.0e20); 
    p->Qo=hydrograph_ipol(p,pgc,hydro_out,hydro_out_count);
    }
    
    
    ULOOP
    a->u(i,j,k)=p->Ui;

    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];


            a->u(i-1,j,k)=p->Ui;
            a->u(i-2,j,k)=p->Ui;
            a->u(i-3,j,k)=p->Ui;

            a->v(i-1,j,k)=0.0;
            a->v(i-2,j,k)=0.0;
            a->v(i-3,j,k)=0.0;

            a->w(i-1,j,k)=0.0;
            a->w(i-2,j,k)=0.0;
            a->w(i-3,j,k)=0.0;
    }

}

