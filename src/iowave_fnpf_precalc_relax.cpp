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

#include"iowave.h"
#include"lexer.h"
#include"ghostcell.h"

void iowave::fnpf_precalc_relax(lexer *p, ghostcell *pgc)
{
    p->wavetime = p->simtime + p->RK_alpha*p->dt;
    
    if(!gen_built) genzone4_build(p,pgc);
    
    // eta and Fifsf of a cell together, so the cell's cached phases are read
    // once (they do not fit in cache for many components). Fifsf is evaluated
    // at the surface eta(i,j) just computed, which the halo exchange below
    // does not change at interior cells; the Fifsfval count sequence is that
    // of the relaxation functions (SLICELOOP4 order).
    count=0;
    
    if(p->B98==2)
    for(size_t q=0; q<gen_i.size(); ++q)
    {
        i = gen_i[q];
        j = gen_j[q];
        
        PSLICECHECK4
        {
        eta(i,j) = wave_eta_c(p,pgc,int(q));
        
        if(f_switch==1)
        {
        Fifsfval[count] = wave_fi_c(p,pgc,int(q),eta(i,j));
        ++count;
        }
        }
    }
    
    pgc->gcsl_start4(p,eta,50);
}
    
