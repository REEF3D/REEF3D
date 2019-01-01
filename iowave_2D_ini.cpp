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

#include"iowave.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"vrans_v.h"
#include"vrans_f.h"

void iowave::ini2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    if((p->B98==1 || p->B98==2) && p->A10==2)
    {
    wavegen_2D_precalc_ini(p,pgc);
    
    wavegen_2D_precalc(p,b,pgc);
    }
    
    if(p->I30==1)
	full_initialize2D(p,b,pgc);
}