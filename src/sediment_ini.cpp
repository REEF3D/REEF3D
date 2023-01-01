/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void sediment_f::ini_cfd(lexer *p, fdm *a,ghostcell *pgc)
{
	double h,h1;

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
        s->bedzh0(i,j)=h;
	}
	
	pgc->gcsl_start4(p,s->bedzh,1);
	
    active_ini_cfd(p,a,pgc);
    
    topo_zh_update(p,a,pgc,s);
    
    log_ini(p);
}

void sediment_f::ini_sflow(lexer *p, fdm2D *b, ghostcell *pgc)
{
    //relax(p,b,pgc);
    
    SLICELOOP4
    {
    s->ks(i,j) = p->S20;
    
    s->bedzh(i,j)=b->topobed(i,j);
    s->bedzh0(i,j)=b->topobed(i,j);
    }
    
    SLICELOOP4
    b->bed(i,j) = MAX(b->topobed(i,j),b->solidbed(i,j));
    
    active_ini_sflow(p,b,pgc);
    
    log_ini(p);
}