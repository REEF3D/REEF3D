/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void partres::bedchange(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, int mode)
{
    // topo
    ALOOP
    a->topo(i,j,k) -= 0.1*(1.0/P.ParcelFactor)*bedch(i,j)*1.0/6.0*PI*pow(P.d50,3.0)/(p->DXN[IP]*p->DYN[JP]*p->S24);
    
    // bedzh
    double h;
    ILOOP
    JLOOP
	{
		KLOOP
		PBASECHECK
		{
        if(a->topo(i,j,k-1)<0.0 && a->topo(i,j,k)>0.0)
        h = -(a->topo(i,j,k-1)*p->DZP[KP])/(a->topo(i,j,k)-a->topo(i,j,k-1)) + p->pos_z()-p->DZP[KP];
        }
		s->bedzh(i,j)=h;
	}
}

void partres::bedchange_update(lexer *p, ghostcell *pgc, sediment_fdm *s, int mode)
{
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        // step 1
        if(mode==1)
        {
        i=p->posc_i(P.X[n]);
        j=p->posc_j(P.Y[n]);
        k=p->posc_k(P.Z[n]);
        }
        
        if(mode==2)
        {
        i=p->posc_i(P.XRK1[n]);
        j=p->posc_j(P.YRK1[n]);
        k=p->posc_k(P.ZRK1[n]);
        }
        
        bedch(i,j) -= P.ParcelFactor;
        
        
        // step 2
        if(mode==1)
        {
        i=p->posc_i(P.XRK1[n]);
        j=p->posc_j(P.YRK1[n]);
        k=p->posc_k(P.ZRK1[n]);
        }
        
        if(mode==2)
        {
        i=p->posc_i(P.X[n]);
        j=p->posc_j(P.Y[n]);
        k=p->posc_k(P.Z[n]);
        }
        
        bedch(i,j) += P.ParcelFactor;
    }
    
    pgc->gcsl_start4(p,bedch,1);
}