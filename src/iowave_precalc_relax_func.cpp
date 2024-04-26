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

void iowave::wavegen_precalc_relax_func(lexer *p, ghostcell *pgc)
{
    // ini fill
    SLICELOOP1
    {
    relax1_wg(i,j) = 1.0;
    relax1_nb(i,j) = 1.0;
    }
    
    SLICELOOP2
    {
    relax2_wg(i,j) = 1.0;
    relax2_nb(i,j) = 1.0;
    }
    
    SLICELOOP4
    {
    relax4_wg(i,j) = 1.0;
    relax4_nb(i,j) = 1.0;
    }
    
    // 1
    SLICELOOP1
    {
		// Wave Generation
        if(p->B98==2)
        {
                relax1_wg(i,j) = rb1_ext(p,1);
		}
        
        // Numerical Beach
        if(p->B99==1 || p->B99==2 || p->B107>0)
        {
                relax1_nb(i,j) = rb3_ext(p,1);
		}
    }
    pgc->gcsl_start1(p,relax1_wg,50);
    pgc->gcsl_start1(p,relax1_nb,50);
    
    // 2
    SLICELOOP2
    {
		// Wave Generation
        if(p->B98==2)
        {
                relax2_wg(i,j) = rb1_ext(p,2);
		}
        
        // Numerical Beach
        if(p->B99==1 || p->B99==2 || p->B107>0)
        {
                relax2_nb(i,j) = rb3_ext(p,2);
		}
    }
    pgc->gcsl_start2(p,relax2_wg,50);
    pgc->gcsl_start2(p,relax2_nb,50);
    
    // 4
    SLICELOOP4
    {
		// Wave Generation
        if(p->B98==2)
        {
                relax4_wg(i,j) = rb1_ext(p,4);
		}
        
        // Numerical Beach
        if(p->B99==1 || p->B99==2 || p->B107>0)
        {
        
                relax4_nb(i,j) = rb3_ext(p,4);
		}
    }
    pgc->gcsl_start4(p,relax4_wg,50);
    pgc->gcsl_start4(p,relax4_nb,50);
}

void iowave::wavegen_precalc_relax_func_fnpf(lexer *p, ghostcell *pgc)
{
    // ini fill
    SLICELOOP4
    {
    relax4_wg(i,j) = 1.0;
    relax4_nb(i,j) = 1.0;
    }
    
    
    // 4
    SLICELOOP4
    {
		// Wave Generation
        if(p->B98==2)
        {
                relax4_wg(i,j) = rb1_ext(p,4);
		}
        
        // Numerical Beach
        if(p->B99==1 || p->B99==2 || p->B107>0)
        {
                relax4_nb(i,j) = rb3_ext(p,4);
		}
    }
    pgc->gcsl_start4(p,relax4_wg,50);
    pgc->gcsl_start4(p,relax4_nb,50);
}

void iowave::wavegen_precalc_relax_func_nhflow(lexer *p, ghostcell *pgc)
{
    // ini fill
    SLICELOOP4
    {
    relax4_wg(i,j) = 1.0;
    relax4_nb(i,j) = 1.0;
    }
    
    
    // 4
    SLICELOOP4
    {
		// Wave Generation
        if(p->B98==2)
        {
                relax4_wg(i,j) = rb1_ext(p,4);
		}
        
        // Numerical Beach
        if(p->B99==1 || p->B99==2 || p->B107>0)
        {
                relax4_nb(i,j) = rb3_ext(p,4);
		}
    }
    pgc->gcsl_start4(p,relax4_wg,50);
    pgc->gcsl_start4(p,relax4_nb,50);
}