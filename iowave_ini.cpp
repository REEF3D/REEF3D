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
#include"fdm.h"
#include"ghostcell.h"
#include"vrans_v.h"
#include"vrans_f.h"
#include"vrans_veg.h"

void iowave::ini(lexer *p, fdm* a, ghostcell* pgc)
{
    if(p->B269==0)
	pvrans = new vrans_v(p,a,pgc);
	
	if(p->B269==1 || p->S10==2)
	pvrans = new vrans_f(p,a,pgc);
    
    if(p->B269==2)
	pvrans = new vrans_veg(p,a,pgc);
    
    // relax_ini OR dirichlet_ini
    if(p->A10==3 || p->A10==4 || p->A10==5)
    {
    wavegen_precalc_ini(p,pgc);
    
    if(p->B89==1 && p->B98==2)
    wavegen_precalc_space(p,pgc);

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
    
    if(p->B89==1 && p->B98==2)
    wavegen_precalc_decomp_space_fnpf(p,pgc);

    wavegen_precalc(p,pgc);

    
    if(p->I30==1)
	full_initialize_fnpf(p,c,pgc);
}

