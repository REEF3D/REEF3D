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

void iowave::fnpf_precalc_relax(lexer *p, ghostcell *pgc)
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

    FILOOP 
    FJLOOP 
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        FKLOOP 
        FPCHECK
        {
        
            z=p->ZSN[FIJK]-p->phimean;

            
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
    
