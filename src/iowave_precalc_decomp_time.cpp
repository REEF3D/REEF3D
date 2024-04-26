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

void iowave::wavegen_precalc_time(lexer *p, ghostcell *pgc)
{
    int qn;
    
    for(qn=0;qn<wave_comp;++qn)
    {
    etaval_T_sin[qn] = wave_eta_time_sin(p,pgc,qn);
    etaval_T_cos[qn] = wave_eta_time_cos(p,pgc,qn);
    }
    
    for(qn=0;qn<wave_comp;++qn)
    {
    uval_T_sin[qn] = wave_u_time_sin(p,pgc,qn);
    uval_T_cos[qn] = wave_u_time_cos(p,pgc,qn);
    }

    for(qn=0;qn<wave_comp;++qn)
    {
    vval_T_sin[qn] = wave_v_time_sin(p,pgc,qn);
    vval_T_cos[qn] = wave_v_time_cos(p,pgc,qn);
    }

    for(qn=0;qn<wave_comp;++qn)
    {
    wval_T_sin[qn] = wave_w_time_sin(p,pgc,qn);
    wval_T_cos[qn] = wave_w_time_cos(p,pgc,qn);
    }

    for(qn=0;qn<wave_comp;++qn)
    {
    Fival_T_sin[qn] = wave_fi_time_sin(p,pgc,qn);
    Fival_T_cos[qn] = wave_fi_time_cos(p,pgc,qn);
    }
    
}
