/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
#include"ghostcell.h"

void iowave::wavegen_precalc_relax_func(lexer *p, ghostcell *pgc)
{

    SLICELOOP1
    {
		dg = distgen(p);
		db = distbeach(p);
		
		// Wave Generation
        if(p->B98==2)
        {
                relax1_wg(i,j) = rb1(p,dg);
		}
        
        // Numerical Beach
        if(p->B98==1 || p->B98==2)
        {
                relax1_nb(i,j) = rb3(p,db);
		}
    }
    pgc->gcsl_start4(p,relax1_wg,50);
}