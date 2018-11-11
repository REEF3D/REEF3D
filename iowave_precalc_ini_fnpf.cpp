/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"iowave.h"
#include"lexer.h"
#include"ghostcell.h"


void iowave::fnpf_precalc_relax_ini(lexer *p, ghostcell *pgc)
{
    // count number of relax points
    // allocate double* array
    
    upt_count=vpt_count=wpt_count=ppt_count=ept_count=0;
    
    if(p->B89==1)
    {
        if(p->B92==5)
        wave_comp = 5;
        
        if(p->B92==31 || p->B92==41 || p->B92==51)
        wave_comp = p->wN;
        
    }
    
    // FI ------------------------------------------------
    FLOOP
    {
		dg = distgen(p); 

		// Wave Generation
        if(p->B98==1)
        {
            // Zone 1
            if(dg<dist1)
            ++ppt_count;

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            ++ppt_count;
		}
		
		if(p->B98==2)
        {
            // Zone 1
            if(dg<dist1)
            ++ppt_count;

		}
    }	

// ETA ------------------------------------------------
    SLICEBASELOOP
    {
		dg = distgen(p); 

		// Wave Generation
        if(p->B98==1)
        {
            // Zone 1
            if(dg<dist1)
            ++ept_count;

            // Zone 2
            if(dg>=dist1 && dg<dist2+3.0*p->dx)
            ++ept_count;
		}
		
		if(p->B98==2)
        {
            // Zone 1
            if(dg<dist1+3.0*p->dx)
            ++ept_count;

		}
    }	
    
    // precalc array alloc
    p->Darray(etaval,ept_count);
    p->Darray(Fival,ppt_count);
    p->Darray(Fifsfval,ept_count);
    
    if(p->B89==1) 
    {
    p->Darray(etaval_S_sin,ept_count,wave_comp);
    p->Darray(Fival_S_sin,ppt_count,wave_comp);
    
    p->Darray(etaval_S_cos,ept_count,wave_comp);
    p->Darray(Fival_S_cos,ppt_count,wave_comp);
    
    p->Darray(etaval_T_sin,wave_comp);
    p->Darray(Fival_T_sin,wave_comp);
    
    p->Darray(etaval_T_cos,wave_comp);
    p->Darray(Fival_T_cos,wave_comp);
    }
}
