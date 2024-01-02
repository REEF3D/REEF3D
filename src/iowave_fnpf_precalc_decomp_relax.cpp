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

void iowave::wavegen_precalc_decomp_relax_fnpf(lexer *p, ghostcell *pgc)
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

    double fsfloc;
    int qn;


    // pre-calc every iteration
    count=0;
    SLICELOOP4
    {
		dg = distgen(p);
		db = distbeach(p);


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

    
    count=0;
    int dbcount=0;
    
    FILOOP 
    FJLOOP 
    {
        dg = distgen(p);
		db = distbeach(p);
        
        FKLOOP 
        FPCHECK
        {
        
        z=p->ZSN[FIJK]-p->phimean;

		
		if(p->B98==2 && f_switch==1)
        {  
            // Zone 1
            if(dg<dist1)
            {
            Fival[count]=0.0;
                        
            for(qn=0;qn<wave_comp;++qn)
            Fival[count] += Fival_S_cos[count][qn]*Fival_T_sin[qn] + Fival_S_sin[count][qn]*Fival_T_cos[qn];
            
            rb1val[count] = rb1(p,dg);
            
            ++count;
            }
		}
        
        if(p->B99==1||p->B99==2)
        {
                // Zone 2
                if(db<dist2)
                {
                rb3val[dbcount] = rb3(p,db);
                ++dbcount;
                }
        }
            
        }
    }
     
    
    // pre-calc every iteration
    count=0;
    SLICELOOP4
    {
		dg = distgen(p);
		db = distbeach(p);
		
		// Wave Generation
        if(p->B98==2 && f_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
                Fifsfval[count] = 0.0;
                
                for(qn=0;qn<wave_comp;++qn)
                Fifsfval[count] += Fifsfval_S_cos[count][qn]*Fifsfval_T_sin[qn] + Fifsfval_S_sin[count][qn]*Fifsfval_T_cos[qn];

            ++count;
            }
		}
    }

    
}
    
