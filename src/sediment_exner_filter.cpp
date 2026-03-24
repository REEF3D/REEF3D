/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"sediment_exner.h"
#include"lexer.h"
#include"ghostcell.h"
#include"bedconc_VR.h"
#include"topo_relax.h"
#include"sediment_exnerdisc.h"
#include"sediment_fdm.h"

void sediment_exner::filter(lexer *p,ghostcell *pgc, slice &f, int outer_iter, int inner_iter)
{
    slice4 h(p),dh(p); 
	
	for(int qn=0;qn<outer_iter;++qn)
	{
		SLICEBASELOOP
        if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
		h(i,j) = f(i,j);
		
		pgc->gcsl_start4(p,h,1);
	
        // predictor
		SLICEBASELOOP
        if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
		f(i,j) = p->S102*h(i,j) + (1.0-p->S102)*0.25*(h(i-1,j) + h(i+1,j) + h(i,j-1) + h(i,j+1));
		
        // corrector
		for(int qqn=0;qqn<inner_iter;++qqn)
		{
            SLICEBASELOOP
            if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
            dh(i,j) = h(i,j) - f(i,j);
            
            
            SLICEBASELOOP
            if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
            dh(i,j) = p->S102*dh(i,j) + (1.0-p->S102)*0.25*(dh(i-1,j) + dh(i+1,j) + dh(i,j-1) + dh(i,j+1));
            
            SLICEBASELOOP
            if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
            f(i,j) += dh(i,j);
		}
    }
}