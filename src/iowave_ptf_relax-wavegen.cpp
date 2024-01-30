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


void iowave::ptf_precalc_relax_ini(lexer *p, ghostcell *pgc)
{
    // count number of relax points
    // allocate double* array

    int dbcount=0;
    
    upt_count=vpt_count=wpt_count=ppt_count=ept_count=0;
    
    if(p->B89==1)
    {
        if(p->B92==5)
        wave_comp = 5;
        
        if(p->B92==31 || p->B92==41 || p->B92==51)
        wave_comp = p->wN;
        
    }
    
    // FI ------------------------------------------------

    FLUIDLOOP
    {
		dg = distgen(p); 
        db = distbeach(p); 

		// Wave Generation
		if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            ++ppt_count;

		}
        
        if(p->B99==1||p->B99==2)
		{
            if(db<1.0e20)
            ++dbcount;
        }
    }	

// ETA ------------------------------------------------
    SLICELOOP4
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
    p->Darray(etaval,ept_count);
    p->Darray(Fival,ppt_count);
    p->Darray(rb1val,ppt_count);
    p->Darray(rb3val,dbcount);
    p->Darray(Fifsfval,ept_count);
    p->Darray(Fifsfval0,ept_count);
    p->Darray(Fifsfval1,ept_count);
    
    
    if(p->B89==1) 
    {
    p->Darray(etaval_S_sin,ept_count,wave_comp);
    p->Darray(Fival_S_sin,ppt_count,wave_comp);
    p->Darray(Fifsfval_S_sin,ept_count,wave_comp);
    p->Darray(uval_S_sin,ppt_count,wave_comp);
    
    p->Darray(etaval_S_cos,ept_count,wave_comp);
    p->Darray(Fival_S_cos,ppt_count,wave_comp);
    p->Darray(Fifsfval_S_cos,ept_count,wave_comp);
    p->Darray(uval_S_cos,ppt_count,wave_comp);
    
    p->Darray(etaval_T_sin,wave_comp);
    p->Darray(Fival_T_sin,wave_comp);
    p->Darray(Fifsfval_T_sin,wave_comp);
    p->Darray(uval_T_sin,wave_comp);
    
    p->Darray(etaval_T_cos,wave_comp);
    p->Darray(Fival_T_cos,wave_comp);
    p->Darray(Fifsfval_T_cos,wave_comp);
    p->Darray(uval_T_cos,wave_comp);
    }

}

void iowave::ptf_precalc_relax(lexer *p, ghostcell *pgc)
{
    double fsfloc;
    int dbcount;
    
    
    // pre-calc every iteration
    // eta
    count=0;
    SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
		
		// Wave Generation
        if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            eta(i,j) = wave_eta(p,pgc,xg,yg);
            etaval[count] = eta(i,j);
            ++count;
            }
		}
    }
    pgc->gcsl_start4(p,eta,50);
    
    
    // Fi
    count=0;
    dbcount=0;

    ILOOP 
    JLOOP 
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        KLOOP 
        PFLUIDCHECK
        {
        
            z=p->ZSN[IJK]-p->phimean;

            
            // Wave Generation
            if(p->B98==2 && f_switch==1)
            {
                // Zone 1
                if(dg<1.0e19)
                { 
                //Fival[count] = wave_fi(p,pgc,xg,yg,z);
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


    count=0;
    SLICELOOP4
    {
		
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        z = eta(i,j);
		
		// Wave Generation
		if(p->B98==2 && f_switch==1)
        {
            // Zone 1
            if(dg<1.0e19)
            { 
            Fifsfval[count] = wave_fi(p,pgc,xg,yg,z);
            
            ++count;
            }
		}
    }
    
}

void iowave::wavegen_precalc_decomp_relax_ptf(lexer *p, ghostcell *pgc)
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
    
    ILOOP 
    JLOOP 
    {
        dg = distgen(p);
		db = distbeach(p);
        
        KLOOP 
        PFLUIDCHECK
        {
        
        z=p->ZSN[IJK]-p->phimean;

		
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
    
