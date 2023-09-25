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
#include"lexer.h"
#include"ghostcell.h"

void iowave::nhflow_wavegen_precalc_decomp_dirichlet(lexer *p, ghostcell *pgc)
{
    int qn;
        
    // ETA
        count=0;
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];

        eta(i,j) = 0.0;
        etaval[count] = 0.0;
            
        for(qn=0;qn<wave_comp;++qn)
        {
        eta(i,j) += etaval_S_cos[count][qn]*etaval_T_cos[qn] - etaval_S_sin[count][qn]*etaval_T_sin[qn];
        etaval[count] = eta(i,j);
        }
        
        ++count;
        }
        
    // -------------
        count=0;
		for(n=0;n<p->gcin_count;n++)
		{
        uval[count] = 0.0;
        vval[count] = 0.0;
        wval[count] = 0.0;
        
    // U
        for(qn=0;qn<wave_comp;++qn)
        uval[count] += uval_S_cos[count][qn]*uval_T_cos[qn] - uval_S_sin[count][qn]*uval_T_sin[qn];

    // V
        for(qn=0;qn<wave_comp;++qn)
        vval[count] += vval_S_cos[count][qn]*vval_T_cos[qn] - vval_S_sin[count][qn]*vval_T_sin[qn];
                    
    // W
        for(qn=0;qn<wave_comp;++qn)
        wval[count] += wval_S_cos[count][qn]*wval_T_sin[qn] + wval_S_sin[count][qn]*wval_T_cos[qn];
        
        ++count;
        }
}




