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

void iowave::wavegen_precalc_decomp_relax(lexer *p, ghostcell *pgc)
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
		
		// Wave Generation
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
    ULOOP
    {
        dg = distgen(p);
        
        zloc1 = p->pos1_z();
        fsfloc = 0.5*(eta(i,j)+eta(i+1,j)) + p->phimean;
    

        if(zloc1<=fsfloc)
        {
        if(zloc1<=p->phimean)
        z=-(fabs(p->phimean-zloc1));
		
		if(zloc1>p->phimean)
        z=(fabs(p->phimean-zloc1));
        }
        
        if(zloc1>fsfloc)
        z = 0.5*(eta(i,j)+eta(i+1,j));
		
		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            
            // Zone 1
            if(dg<dist1)
            {
            uval[count]=0.0;
            
            if(zloc1<=fsfloc+epsi)
            for(qn=0;qn<wave_comp;++qn)
            uval[count] += uval_S_cos[count][qn]*uval_T_cos[qn] - uval_S_sin[count][qn]*uval_T_sin[qn];
            
            if(zloc1>fsfloc+epsi)
            uval[count] = 0.0;
            
            ++count;
            }
            
		}
    }


    count=0;
    VLOOP
    {
        dg = distgen(p);
        
        zloc2 = p->pos2_z();
        fsfloc = 0.5*(eta(i,j)+eta(i,j+1)) + p->phimean;
    

        if(zloc2<=fsfloc)
        {
        if(zloc2<=p->phimean)
        z=-(fabs(p->phimean-zloc2));
		
		if(zloc2>p->phimean)
        z=(fabs(p->phimean-zloc2));
        }
        
        if(zloc2>fsfloc)
        z = 0.5*(eta(i,j)+eta(i,j+1));

		// Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            vval[count]=0.0;
            
            if(zloc2<=fsfloc+epsi)
            for(qn=0;qn<wave_comp;++qn)
            vval[count] += vval_S_cos[count][qn]*vval_T_cos[qn] - vval_S_sin[count][qn]*vval_T_sin[qn];
            
            if(zloc2>fsfloc+epsi)
            vval[count] = 0.0;
            
            ++count;
            }
		}
    }


    count=0;
    WLOOP
    {
        dg = distgen(p);
        
        zloc3 = p->pos3_z();
        fsfloc = eta(i,j) + p->phimean;
    
        if(zloc3<=fsfloc)
        {
        if(zloc3<=p->phimean)
        z=-(fabs(p->phimean-zloc3));
		
		if(zloc3>p->phimean)
        z=(fabs(p->phimean-zloc3));
        }
        
        if(zloc3>fsfloc)
        z = eta(i,j);


		// Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            wval[count]=0.0;
            
            if(zloc3<=fsfloc+epsi)
            for(qn=0;qn<wave_comp;++qn)
            wval[count] += wval_S_cos[count][qn]*wval_T_sin[qn] + wval_S_sin[count][qn]*wval_T_cos[qn];
            
            if(zloc3>fsfloc+epsi)
            wval[count] = 0.0;
            
            ++count;
            }
		}
    }	
    
    count=0;
    LOOP
    {
		dg = distgen(p);
        
		if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos_z()));
		
		if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos_z()));
		
		// Wave Generation
        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            lsval[count] = eta(i,j)+p->phimean-p->pos_z();
            
            ++count;
            }
		}
    }

    
    count=0;
    LOOP
    {
        dg = distgen(p);
        
        zloc4 = p->pos_z();
        fsfloc = eta(i,j) + p->phimean;
    

        if(zloc4<=fsfloc)
        {
        if(zloc4<=p->phimean)
        z=-(fabs(p->phimean-zloc4));
		
		if(zloc4>p->phimean)
        z=(fabs(p->phimean-zloc4));
        }
        
        if(zloc4>fsfloc)
        z = eta(i,j);
		
		// Wave Generation		
		if(p->B98==2 && u_switch==1)
        {  
            // Zone 1
            if(dg<dist1)
            {
            Fival[count]=0.0;
            
            if(zloc4<=fsfloc+epsi)
            for(qn=0;qn<wave_comp;++qn)
            Fival[count] += Fival_S_cos[count][qn]*Fival_T_cos[qn] - Fival_S_sin[count][qn]*Fival_T_sin[qn];
            
            if(zloc4>fsfloc+epsi)
            Fival[count] = 0.0;
            
            ++count;
            }
		}
    }
    
}
    
