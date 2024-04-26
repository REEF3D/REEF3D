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
#include"fdm2D.h"
#include"ghostcell.h"

void iowave::wavegen_2D_precalc(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double fsfloc;
    double u_val,v_val,w_val;
    double deltaz;
    
    // pre-calc every iteration
    count=0;
    SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
		
		// Wave Generation
        if(p->B98==2 && h_switch==1)
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
    
    count=0;
    SLICELOOP1
    {
		xg = xgen1(p);
        yg = ygen1(p);
        dg = distgen(p);
		db = distbeach(p);
        
        
        deltaz = (0.5*(eta(i,j)+eta(i+1,j)) + p->wd - 0.5*(b->bed(i,j)+b->bed(i+1,j)))/(double(p->B160));
        
        u_val=0.0;
        z=-p->wd;
        for(int qn=0;qn<=p->B160;++qn)
        {
        u_val += wave_u(p,pgc,xg,yg,z);
        
        z+=deltaz;
        }
        u_val/=double(p->B160+1);
		
		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            uval[count] = u_val;
            ++count;
            }
		}
    }
		
    count=0;
    SLICELOOP2
    {
        xg = xgen2(p);
        yg = ygen2(p);
        dg = distgen(p);
		db = distbeach(p);
        
        deltaz = (0.5*(eta(i,j)+eta(i,j+1)) + p->wd - 0.5*(b->bed(i,j)+b->bed(i,j+1)))/(double(p->B160));
        
        v_val=0.0;
        z=-p->wd;
        for(int qn=0;qn<=p->B160;++qn)
        {
        v_val += wave_v(p,pgc,xg,yg,z);
        
        z+=deltaz;
        }
        v_val/=double(p->B160+1);
        
		// Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            vval[count] = v_val;
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

        deltaz = (eta(i,j) + p->wd - b->bed(i,j))/(double(p->B160));
        
        w_val=0.0;
        z=-p->wd;
        for(int qn=0;qn<=p->B160;++qn)
        {
        w_val += wave_w(p,pgc,xg,yg,z);
        
        z+=deltaz;
        }
        w_val/=double(p->B160+1);
        
        
		// Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            wval[count] = w_val;
            ++count;
            }
		}
    }
    
}
    
