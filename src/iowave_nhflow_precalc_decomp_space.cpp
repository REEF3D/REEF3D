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

void iowave::nhflow_wavegen_precalc_decomp_space(lexer *p, ghostcell *pgc)
{
    int qn;

// ETA
    count=0;
    SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
		
		// Wave Generation
        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
                for(qn=0;qn<wave_comp;++qn)
                {
                etaval_S_sin[count][qn] = wave_eta_space_sin(p,pgc,xg,yg,qn);
                etaval_S_cos[count][qn] = wave_eta_space_cos(p,pgc,xg,yg,qn);
                }
            ++count;
            }
		}
    }
    pgc->gcsl_start4(p,eta,50);
    
// U
    count=0;
    LOOP
    {
		xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
        
        z=p->ZSP[IJK]-p->phimean;

		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
                for(qn=0;qn<wave_comp;++qn)
                {
                uval_S_sin[count][qn] = wave_u_space_sin(p,pgc,xg,yg,z,qn);
                uval_S_cos[count][qn] = wave_u_space_cos(p,pgc,xg,yg,z,qn);
                }
            ++count;
            }
		}
    }

// V
    count=0;    
    LOOP
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
        
        z=p->ZSP[IJK]-p->phimean;
        
		// Wave Generation		
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
                for(qn=0;qn<wave_comp;++qn)
                {
                vval_S_sin[count][qn] = wave_v_space_sin(p,pgc,xg,yg,z,qn);
                vval_S_cos[count][qn] = wave_v_space_cos(p,pgc,xg,yg,z,qn);
                }
            ++count;
            }
		}
    }

// W
    count=0;
    LOOP
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
        
        z=p->ZSP[IJK]-p->phimean;
        
		// Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
                for(qn=0;qn<wave_comp;++qn)
                {
                wval_S_sin[count][qn] = wave_w_space_sin(p,pgc,xg,yg,z,qn);
                wval_S_cos[count][qn] = wave_w_space_cos(p,pgc,xg,yg,z,qn);
                }
            ++count;
            }
		}
    }	
    
}
