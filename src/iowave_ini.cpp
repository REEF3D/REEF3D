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
#include"fdm.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::ini(lexer *p, fdm* a, ghostcell* pgc)
{
    // relax_ini OR dirichlet_ini
    if(p->A10==6)
    {
    wavegen_precalc_ini(p,pgc);
    wavegen_precalc_relax_func(p,pgc);
    
    if(p->B89==1 && p->B98==2)
    wavegen_precalc_space(p,pgc);
    
    if(p->B89==1 && p->B98>=3)
    wavegen_precalc_space_dirichlet(p,pgc);

    wavegen_precalc(p,pgc);
    
    u_relax(p,a,pgc,a->u);
	v_relax(p,a,pgc,a->v);
	w_relax(p,a,pgc,a->w);
    }
    
    if(p->I30==1)
	full_initialize(p,a,pgc);
}

void iowave::ini_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    wavegen_precalc_ini(p,pgc);
    wavegen_precalc_relax_func_fnpf(p,pgc);
    
    if(p->B89==1 && p->B98==2)
    wavegen_precalc_decomp_space_fnpf(p,pgc);
    
    if(p->B89==1 && p->B98==3)
    wavegen_precalc_decomp_space_dirichlet_fnpf(p,pgc);

    wavegen_precalc(p,pgc);

    if(p->I30==1)
	full_initialize_fnpf(p,c,pgc);
}

void iowave::ini_ptf(lexer *p, fdm *a, ghostcell *pgc)
{
    wavegen_precalc_ini(p,pgc);
    wavegen_precalc_relax_func_nhflow(p,pgc);
    
    //if(p->B89==1 && p->B98==2)
    //wavegen_precalc_decomp_space_fnpf(p,pgc);

    wavegen_precalc(p,pgc);
    
    if(p->I30==1)
	full_initialize_ptf(p,a,pgc);
}

