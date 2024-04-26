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
#include"fdm.h"
#include"ghostcell.h"

void iowave::wavegen_2D_precalc_ini(lexer *p, ghostcell *pgc)
{
    // count number of relax points
    // allocate double* array
    upt_count=vpt_count=wpt_count=ppt_count=ept_count=0;
    
    // U ------------------------------------------------
    SLICEBASELOOP
    {
        dg = distgen(p);
        
        // Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            ++upt_count;
		}
    }
    
    
    // V ------------------------------------------------
    SLICEBASELOOP
    {
		dg = distgen(p);

        
		// Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            ++vpt_count;
		}
    }
    
    // W ------------------------------------------------
    SLICEBASELOOP
    {
		dg = distgen(p); 

		// Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            ++wpt_count;

		}
    }

// ETA ------------------------------------------------
    SLICEBASELOOP
    {
		dg = distgen(p); 

		// Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            ++ept_count;

		}
    }	
    
    //cout<<"EPT_COUNT: "<<ept_count;
    
    // precalc array alloc
    upt_count *= (p->B160+1); 
    vpt_count *= (p->B160+1);
    wpt_count *= (p->B160+1);
    
    p->Darray(uval,upt_count);
    p->Darray(vval,upt_count);
    p->Darray(wval,upt_count);
    p->Darray(etaval,ept_count);

    
    if((p->B92==31 || p->B92==41 || p->B92==51 ) && p->B89==1) 
    {
    p->Darray(uval_S_sin,upt_count,p->wN);
    p->Darray(vval_S_sin,vpt_count,p->wN);
    p->Darray(wval_S_sin,wpt_count,p->wN);
    p->Darray(etaval_S_sin,ept_count,p->wN);
    p->Darray(Fival_S_sin,ppt_count,p->wN);
    
    p->Darray(uval_S_cos,upt_count,p->wN);
    p->Darray(vval_S_cos,vpt_count,p->wN);
    p->Darray(wval_S_cos,wpt_count,p->wN);
    p->Darray(etaval_S_cos,ept_count,p->wN);
    p->Darray(Fival_S_cos,ppt_count,p->wN);
    
    p->Darray(uval_T_sin,p->wN);
    p->Darray(vval_T_sin,p->wN);
    p->Darray(wval_T_sin,p->wN);
    p->Darray(etaval_T_sin,p->wN);
    p->Darray(Fival_T_sin,p->wN);
    
    p->Darray(uval_T_cos,p->wN);
    p->Darray(vval_T_cos,p->wN);
    p->Darray(wval_T_cos,p->wN);
    p->Darray(etaval_T_cos,p->wN);
    p->Darray(Fival_T_cos,p->wN);
    }

}
