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
#include"looping.h"

void iowave::wavegen_precalc_ini(lexer *p, ghostcell *pgc)
{
    if(p->A10!=3)
    {
        if(p->B98==2)
        wavegen_precalc_relax_ini(p,pgc);
        
        if(p->B98==3 || p->B98==4)
        wavegen_precalc_dirichlet_ini(p,pgc);
    }
    
    if(p->A10==3) // FNPF
    {
        if(p->B98==2)
        fnpf_precalc_relax_ini(p,pgc);
        
        if(p->B98==3 || p->B98==4)
        fnpf_precalc_dirichlet_ini(p,pgc);
    }
    
    if(p->A10==55) // NHFLOW
    {
        if(p->B98==2)
        nhflow_precalc_relax_ini(p,pgc);
        
        if(p->B98==3 || p->B98==4)
        nhflow_precalc_dirichlet_ini(p,pgc);
    }
}

void iowave::wavegen_precalc_relax_ini(lexer *p, ghostcell *pgc)
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
    UBASELOOP
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
    
    
    // U ------------------------------------------------
    VBASELOOP
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
    WBASELOOP
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
    p->Darray(etaval,ept_count);
    p->Darray(lsval,ppt_count);
    p->Darray(Fival,ppt_count);
    p->Darray(Fifsfval,ept_count);
    
    if(p->B89==1) 
    {
    p->Darray(uval_S_sin,upt_count,wave_comp);
    p->Darray(vval_S_sin,vpt_count,wave_comp);
    p->Darray(wval_S_sin,wpt_count,wave_comp);
    p->Darray(etaval_S_sin,ept_count,wave_comp);
    p->Darray(Fival_S_sin,ppt_count,wave_comp);
    
    p->Darray(uval_S_cos,upt_count,wave_comp);
    p->Darray(vval_S_cos,vpt_count,wave_comp);
    p->Darray(wval_S_cos,wpt_count,wave_comp);
    p->Darray(etaval_S_cos,ept_count,wave_comp);
    p->Darray(Fival_S_cos,ppt_count,wave_comp);
    
    p->Darray(uval_T_sin,wave_comp);
    p->Darray(vval_T_sin,wave_comp);
    p->Darray(wval_T_sin,wave_comp);
    p->Darray(etaval_T_sin,wave_comp);
    p->Darray(Fival_T_sin,wave_comp);
    
    p->Darray(uval_T_cos,wave_comp);
    p->Darray(vval_T_cos,wave_comp);
    p->Darray(wval_T_cos,wave_comp);
    p->Darray(etaval_T_cos,wave_comp);
    p->Darray(Fival_T_cos,wave_comp);
    }
}
void iowave::wavegen_precalc_dirichlet_ini(lexer *p, ghostcell *pgc)
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
    p->Darray(etaval,ept_count);
    p->Darray(lsval,ppt_count);
    p->Darray(Fival,ppt_count);
    p->Darray(Fifsfval,ept_count);
    
    if(p->B89==1) 
    {
    p->Darray(uval_S_sin,upt_count,wave_comp);
    p->Darray(vval_S_sin,vpt_count,wave_comp);
    p->Darray(wval_S_sin,wpt_count,wave_comp);
    p->Darray(etaval_S_sin,ept_count,wave_comp);
    p->Darray(Fival_S_sin,ppt_count,wave_comp);
    
    p->Darray(uval_S_cos,upt_count,wave_comp);
    p->Darray(vval_S_cos,vpt_count,wave_comp);
    p->Darray(wval_S_cos,wpt_count,wave_comp);
    p->Darray(etaval_S_cos,ept_count,wave_comp);
    p->Darray(Fival_S_cos,ppt_count,wave_comp);
    
    p->Darray(uval_T_sin,wave_comp);
    p->Darray(vval_T_sin,wave_comp);
    p->Darray(wval_T_sin,wave_comp);
    p->Darray(etaval_T_sin,wave_comp);
    p->Darray(Fival_T_sin,wave_comp);
    
    p->Darray(uval_T_cos,wave_comp);
    p->Darray(vval_T_cos,wave_comp);
    p->Darray(wval_T_cos,wave_comp);
    p->Darray(etaval_T_cos,wave_comp);
    p->Darray(Fival_T_cos,wave_comp);
    }
}
