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

void iowave::fnpf_precalc_relax(lexer *p, ghostcell *pgc)
{
    double fsfloc;
    
    
    // pre-calc every iteration
    count=0;
    SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
		
		// Wave Generation
        if(p->B98==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            eta(i,j) = wave_eta(p,pgc,xg,yg);
            etaval[count] = eta(i,j);
            ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2+3.0*p->dx)
            {
            eta(i,j) = wave_eta(p,pgc,xg,yg);
            etaval[count] = eta(i,j);
            ++count;
            }
		}

        if(p->B98==2)
        {
            // Zone 1
            if(dg<dist1+3.0*p->DXM)
            {
            eta(i,j) = wave_eta(p,pgc,xg,yg);
            etaval[count] = eta(i,j);
            ++count;
            }
		}
    }
    pgc->gcsl_start4(p,eta,50);
    
    
    count=0;
    FLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        z=p->ZSN[FIJK]-p->phimean;

        
		// Wave Generation
        if(p->B98==1 && f_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            Fival[count] = wave_fi(p,pgc,xg,yg,z);
            ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            { 
            Fival[count] = wave_fi(p,pgc,xg,yg,z);
            ++count;
            }
		}
		
		if(p->B98==2 && f_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            { 
            Fival[count] = wave_fi(p,pgc,xg,yg,z);
            ++count;
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
        if(p->B98==1 && f_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            Fifsfval[count] = wave_fi(p,pgc,xg,yg,z);
            
            ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            { 
            Fifsfval[count] = wave_fi(p,pgc,xg,yg,z);
            
            ++count;
            }
		}
		
		if(p->B98==2 && f_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            { 
            Fifsfval[count] = wave_fi(p,pgc,xg,yg,z);
            
            ++count;
            }
		}
    }
    
}
    
