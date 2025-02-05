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
#include"fdm_fnpf.h"
#include"ghostcell.h"

void iowave::wavegen_precalc_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    starttime=pgc->timer();
    
    // prestep
    wave_prestep(p,pgc);
    

        if(p->B89==0 )
        {
            if(p->B98==2)
            fnpf_precalc_relax(p,pgc);
            
            if(p->B98==3 || p->B98==4)
            fnpf_precalc_dirichlet(p,pgc);
        }
        
        if(p->B89==1)
        {
            if(p->B98==2)
            {
            wavegen_precalc_decomp_time_fnpf(p,pgc);
            wavegen_precalc_decomp_relax_fnpf(p,pgc);
            }
            
            if(p->B98==3 || p->B98==4)
            {
            wavegen_precalc_decomp_time_fnpf(p,pgc);
            wavegen_precalc_decomp_dirichlet_fnpf(p,pgc);
            }
        }

    p->wavecalctime+=pgc->timer()-starttime;
}
    
