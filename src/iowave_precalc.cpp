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
#include"ghostcell.h"

void iowave::wavegen_precalc(lexer *p, ghostcell *pgc)
{
    
    // prestep
    wave_prestep(p,pgc);
    
    if(p->A10!=3 && p->A10!=5)
    {
        if(p->B89==0 )
        {
            if(p->B98==2)
            wavegen_precalc_relax(p,pgc);
            
            if(p->B98==3 || p->B98==4)
            wavegen_precalc_dirichlet(p,pgc);
        }
        
        if(p->B89==1)
        {
            if(p->B98==2)
            {
            wavegen_precalc_time(p,pgc);
            wavegen_precalc_decomp_relax(p,pgc);
            }
            
            if(p->B98==3 || p->B98==4)
            wavegen_precalc_dirichlet(p,pgc);
        }
    }
    
    if(p->A10==3)
    {
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
    }
    
    
    
}
    
