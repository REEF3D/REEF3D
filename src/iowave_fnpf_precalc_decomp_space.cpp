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

void iowave::wavegen_precalc_decomp_space_fnpf(lexer *p, ghostcell *pgc)
{
    double fsfloc;
    int qn;


    // pre-calc every iteration
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
    

    //FI
    count=0;
    FLOOP
    {
		
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        z=p->ZSN[FIJK]-p->phimean;
		
		// Wave Generation
		if(p->B98==2)
        {
            // Zone 1
            if(dg<dist1)
            {
                for(qn=0;qn<wave_comp;++qn)
                {
                Fival_S_sin[count][qn] = wave_fi_space_sin(p,pgc,xg,yg,z,qn);
                Fival_S_cos[count][qn] = wave_fi_space_cos(p,pgc,xg,yg,z,qn);
                }
            ++count;
            }
		}
    }
    
    
    
    // Fifsf
    count=0;
    SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
        
        z = eta(i,j);
		
		// Wave Generation
        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
                for(qn=0;qn<wave_comp;++qn)
                {
                Fifsfval_S_sin[count][qn] = wave_fi_space_sin(p,pgc,xg,yg,z,qn);
                Fifsfval_S_cos[count][qn] = wave_fi_space_cos(p,pgc,xg,yg,z,qn);
                }
            ++count;
            }
		}
    }



}
