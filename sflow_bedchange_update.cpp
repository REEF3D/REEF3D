/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"sflow_sediment_f.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"fnpf_weno.h"

void sflow_sediment_f::bedchange_update(lexer *p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
    b->zb(i,j) = b->bed(i,j) - b->bed0(i,j);
    
    pgc->gcsl_start4(p,b->zb,50);
    b->zb.ggcpol(p);
}