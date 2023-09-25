/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#include"fdm_nhf.h"
#include"lexer.h"
#include"ghostcell.h"

void iowave::wavegen_precalc_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    wave_prestep(p,pgc);
    
    if(p->B89==0)
    {
        if(p->B98==2)
        nhflow_precalc_relax(p,d,pgc);
                
        if(p->B98==3 || p->B98==4)
        nhflow_precalc_dirichlet(p,d,pgc);
    }
    
    if(p->B89==1)
    {
        if(p->B98==2)
        {
        nhflow_wavegen_precalc_decomp_time(p,pgc);
        nhflow_wavegen_precalc_decomp_relax(p,d,pgc);
        }
            
        if(p->B98==3 || p->B98==4)
        nhflow_wavegen_precalc_decomp_dirichlet(p,pgc);
    }
}

void iowave::wavegen_precalc_ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    wave_prestep(p,pgc);
    
    if(p->B98==2)
    nhflow_precalc_relax_ini(p,d,pgc);
        
    if(p->B98==3 || p->B98==4)
    nhflow_precalc_dirichlet_ini(p,d,pgc);
}
