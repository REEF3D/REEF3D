/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

void iowave::nhflow_precalc_relax(lexer *p, ghostcell *pgc)
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
    
    
    // U
    count=0;
    ULOOP
    {
		xg = xgen1(p);
        yg = ygen1(p);
        dg = distgen(p);
		db = distbeach(p);
        
        z=p->ZSP[IJK]-p->phimean;
		
		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            uval[count] = wave_u(p,pgc,xg,yg,z);
            ++count;
            }
		}
    }
	
    // V	
    count=0;
    VLOOP
    {
        xg = xgen2(p);
        yg = ygen2(p);
        dg = distgen(p);
		db = distbeach(p);
        
        z=p->ZSP[IJK]-p->phimean;

		// Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            vval[count] = wave_v(p,pgc,xg,yg,z);
            ++count;
            }
		}
    }
    
    // W
    count=0;
    WLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        z=p->ZSN[(i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + (k+1)-p->kmin]-p->phimean;


		// Wave Generation		
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            if(zloc3<=fsfloc+epsi)
            wval[count] = wave_w(p,pgc,xg,yg,z);
            
            if(zloc3>fsfloc+epsi)
            wval[count] = 0.0;
            
            ++count;
            }
		}
    }	
    
}
    
