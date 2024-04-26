/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::ini_nhflow(lexer *p, fdm_nhf *d, ghostcell* pgc)
{
    // relax_ini or dirichlet_ini
    wavegen_precalc_ini_nhflow(p,d,pgc);
    wavegen_precalc_relax_func_nhflow(p,pgc);
    
    if(p->B89==1 && p->B98==2)
    nhflow_wavegen_precalc_decomp_space(p,pgc);
    
    if(p->B89==1 && p->B98>=3)
    nhflow_wavegen_precalc_decomp_space_dirichlet(p,pgc);
    
    wavegen_precalc(p,pgc);
    
    U_relax(p,pgc,d->U,d->UH);
    V_relax(p,pgc,d->V,d->VH);
    W_relax(p,pgc,d->W,d->WH);

    if(p->I30==1)
	full_initialize_nhflow(p,d,pgc);
}