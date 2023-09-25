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
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::nhflow_precalc_relax(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    double fsfloc;
    
// ETA
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
    LOOP
    {
		xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        z=p->ZSP[IJK]-p->phimean;
		
		// Wave Generation
		if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            uval[count] = wave_u(p,pgc,xg,yg,z) + p->Ui;
            UHval[count] = (eta(i,j) + d->depth(i,j))*uval[count];
            ++count;
            }
		}
    }
		
// V
    count=0;
    if(p->j_dir==1)
    LOOP
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        z=p->ZSP[IJK]-p->phimean;
        
		// Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            vval[count] = wave_v(p,pgc,xg,yg,z);
            VHval[count] = (eta(i,j) + d->depth(i,j))*vval[count];
            ++count;
            }
		}
    }
    
// W
    count=0;
    LOOP
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        zloc3 = p->pos3_z();
        fsfloc = eta(i,j) + p->phimean;

        z=p->ZSP[IJK]-p->phimean;

		// Wave Generation		
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            wval[count] = wave_w(p,pgc,xg,yg,z);
            WHval[count] = (eta(i,j) + d->depth(i,j))*wval[count];
            ++count;
            }
		}
    }	
}
    
