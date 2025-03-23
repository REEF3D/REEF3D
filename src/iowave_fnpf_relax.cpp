/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

void iowave::fifsf_relax(lexer *p, ghostcell *pgc, slice& f)
{
    starttime=pgc->timer();
    
    count=0;
    SLICELOOP4
    {
        dg = distgen(p);    
        db = distbeach(p);
        
		// Wave Generation
		if(p->B98==2 && f_switch==1)
        {
            if(dg<1.0e20)
            {
            f(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*Fifsfval[count]  + relax4_wg(i,j)*f(i,j);
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->A10!=3 || p->A348==1 || p->A348==3)
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            f(i,j) = relax4_nb(i,j)*f(i,j);
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}