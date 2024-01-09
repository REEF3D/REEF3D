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

void iowave::wavegen_precalc_decomp_dirichlet_fnpf(lexer *p, ghostcell *pgc)
{
    int qn;
        
        // eta and Fifsfval
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
        
        z = eta(i,j);
        
        time_1=time_0;
        time_0=p->simtime;
        time_n=p->simtime+p->dt;
        Fifsfval1[count] = Fifsfval0[count];
        Fifsfval0[count] = Fifsfval[count];
        
        
        Fifsfval[count]=0.0;
        
        
        for(qn=0;qn<wave_comp;++qn)
        Fifsfval[count] += Fifsfval_S_cos[count][qn]*Fifsfval_T_cos[qn] - Fifsfval_S_sin[count][qn]*Fifsfval_T_sin[qn];
        
        ++count;
        }
        
        // Uin
        count=0;
		for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
            Uinval[count]=0.0;
        
            FKLOOP
            FPCHECK
            {
            for(qn=0;qn<wave_comp;++qn)
            Uinval[count] += uval_S_cos[count][qn]*uval_T_cos[qn] - uval_S_sin[count][qn]*uval_T_sin[qn];
            
            ++count;
            }
        }
        
    
    // beach
    count=0;
    FILOOP 
    FJLOOP 
    {

		db = distbeach(p);
        
        FKLOOP 
        FPCHECK
        {
                    
            if(p->B99==1||p->B99==2)
            {
                // Zone 2
                if(db<dist2)
                {
                rb3val[count] = rb3(p,db);
                ++count;
                }
            }
        }
    }
   

}
