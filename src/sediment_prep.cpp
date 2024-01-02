/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

void sediment_f::prep_cfd(lexer *p, fdm *a,ghostcell *pgc)
{    
    
    // vel prep --------
    pgc->start1(p,a->u,14);
	pgc->start2(p,a->v,15);
	pgc->start3(p,a->w,16);
    
    // find bedk -------
    fill_bedk(p,a,pgc);
    
    fill_PQ_cfd(p,a,pgc);
    
    waterlevel(p,a,pgc);
    
}

void sediment_f::prep_sflow(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q)
{
    
    fill_PQ_sflow(p,b,pgc,P,Q);
    
}