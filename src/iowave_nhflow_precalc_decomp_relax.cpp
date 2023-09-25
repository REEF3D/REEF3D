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
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::nhflow_wavegen_precalc_decomp_relax(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
        /*
        a: space
        b: time
        
        sin(a + b) = sin(a) cos(b) + cos(a) sin(b)
        cos(a + b) = cos(a) cos(b) - sin(a) sin(b)
        */

        /*
         U: cos()
         V: cos()
         W: sin()
         ETA: cos()
        */

    int qn;


// ETA
    count=0;
    SLICELOOP4
    {
		dg = distgen(p);

        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
                eta(i,j) = 0.0;
                etaval[count] = 0.0;
                
                for(qn=0;qn<wave_comp;++qn)
                {
                eta(i,j) += etaval_S_cos[count][qn]*etaval_T_cos[qn] - etaval_S_sin[count][qn]*etaval_T_sin[qn];
                etaval[count] = eta(i,j);
                }
            ++count;
            }
		}
    }
    pgc->gcsl_start4(p,eta,50);

    
    int dbcount=0;
    
// U
    count=0;
    LOOP
    {
        dg = distgen(p);
		
		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            
            // Zone 1
            if(dg<1.0e20)
            {
            uval[count]=0.0;
            
            for(qn=0;qn<wave_comp;++qn)
            uval[count] += uval_S_cos[count][qn]*uval_T_cos[qn] - uval_S_sin[count][qn]*uval_T_sin[qn];
            
            UHval[count] = (eta(i,j) + d->depth(i,j))*uval[count];
            
            ++count;
            }
		}
    }


    count=0;
    LOOP
    {
        dg = distgen(p);
        
		// Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            vval[count]=0.0;
            
            for(qn=0;qn<wave_comp;++qn)
            vval[count] += vval_S_cos[count][qn]*vval_T_cos[qn] - vval_S_sin[count][qn]*vval_T_sin[qn];
            
            VHval[count] = (eta(i,j) + d->depth(i,j))*vval[count];
            
            ++count;
            }
		}
    }


    count=0;
    LOOP
    {
        dg = distgen(p);

		// Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            wval[count]=0.0;

            for(qn=0;qn<wave_comp;++qn)
            wval[count] += wval_S_cos[count][qn]*wval_T_sin[qn] + wval_S_sin[count][qn]*wval_T_cos[qn];
            
            WHval[count] = (eta(i,j) + d->depth(i,j))*wval[count];

            ++count;
            }
		}
    }	
    


    
}
    
