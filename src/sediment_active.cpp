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

void sediment_f::active_cfd(lexer *p, fdm *a,ghostcell *pgc)
{
    
    SLICELOOP4
    s->active(i,j)=0;
    
    // #define ALOOP ILOOP JLOOP KLOOP PSOLIDCHECK
    ILOOP
    JLOOP
    {
        KWLOOP
        PSOLIDCHECK
        {
            
        if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
        s->active(i,j)=1;
        }
    }
}

void sediment_f::active_ini_cfd(lexer *p, fdm *a,ghostcell *pgc)
{
    SLICELOOP4
    s->active(i,j)=1;
}

void sediment_f::active_sflow(lexer *p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
    s->active(i,j)=1;
    
    SLICELOOP4
    if(b->solidbed(i,j) >= s->bedzh(i,j))
    {
    b->test(i,j) = b->solidbed(i,j);

    s->active(i,j)=0;
    }
    
}

void sediment_f::active_ini_sflow(lexer *p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
    s->active(i,j)=1;
    
  
  
    
}