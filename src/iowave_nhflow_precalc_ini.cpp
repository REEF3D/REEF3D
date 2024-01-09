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

void iowave::nhflow_precalc_relax_ini(lexer *p,fdm_nhf *d, ghostcell *pgc)
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
    
    // U ------------------------------------------------
    BASELOOP
    {
        dg = distgen(p);
        
        // Wave Generation
		if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            ++upt_count;
		}
    }
    
    
    // V ------------------------------------------------
    BASELOOP
    {
		dg = distgen(p);

        
		// Wave Generation
		if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            ++vpt_count;
		}
    }
    
    // W ------------------------------------------------
    BASELOOP
    {
		dg = distgen(p); 

		// Wave Generation
		if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            ++wpt_count;

		}
    }

    // FI ------------------------------------------------
    FBASELOOP
    {
		dg = distgen(p); 

		// Wave Generation
		if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            ++ppt_count;

		}
    }	

// ETA ------------------------------------------------
    SLICEBASELOOP
    {
		dg = distgen(p); 

		// Wave Generation
		if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            ++ept_count;

		}
    }	
    
    // precalc array alloc
    p->Darray(uval,upt_count);
    p->Darray(vval,vpt_count);
    p->Darray(wval,wpt_count);
    p->Darray(UHval,upt_count);
    p->Darray(VHval,vpt_count);
    p->Darray(WHval,wpt_count);
    p->Darray(etaval,ept_count);
    
    
    
    if(p->B89==1) 
    {
    p->Darray(uval_S_sin,upt_count,wave_comp);
    p->Darray(vval_S_sin,vpt_count,wave_comp);
    p->Darray(wval_S_sin,wpt_count,wave_comp);
    p->Darray(etaval_S_sin,ept_count,wave_comp);

    p->Darray(uval_S_cos,upt_count,wave_comp);
    p->Darray(vval_S_cos,vpt_count,wave_comp);
    p->Darray(wval_S_cos,wpt_count,wave_comp);
    p->Darray(etaval_S_cos,ept_count,wave_comp);

    p->Darray(uval_T_sin,wave_comp);
    p->Darray(vval_T_sin,wave_comp);
    p->Darray(wval_T_sin,wave_comp);
    p->Darray(etaval_T_sin,wave_comp);
    
    p->Darray(uval_T_cos,wave_comp);
    p->Darray(vval_T_cos,wave_comp);
    p->Darray(wval_T_cos,wave_comp);
    p->Darray(etaval_T_cos,wave_comp);
    }
}

void iowave::nhflow_precalc_dirichlet_ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // count number of relax points
    // allocate double* array
    
    upt_count=vpt_count=wpt_count=ppt_count=ept_count = p->gcin_count;
    
    if(p->B89==1)
    {
        if(p->B92==5)
        wave_comp = 5;
        
        if(p->B92==31 || p->B92==41 || p->B92==51)
        wave_comp = p->wN;
    }
  
    // precalc array alloc
    p->Darray(uval,upt_count);
    p->Darray(vval,vpt_count);
    p->Darray(wval,wpt_count);
    p->Darray(UHval,upt_count);
    p->Darray(VHval,vpt_count);
    p->Darray(WHval,wpt_count);
    p->Darray(etaval,ept_count);
}
