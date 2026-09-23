/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
    
    p->wavetime = p->simtime + p->RK_alpha*p->dt;
    
    // pre-calc every iteration
    // eta
    // only relaxation-zone cells, SLICELOOP4 order (cached geometry)
    if(!rz4_built) relaxzone4_build(p);
    for(size_t rzq=0; rzq<rz4_i.size(); ++rzq)
    {
        i = rz4_i[rzq];
        j = rz4_j[rzq];
        xg = rz4_xg[rzq];
        yg = rz4_yg[rzq];
        dg = rz4_dg[rzq];
        db = rz4_db[rzq];
		
		// Wave Generation
        if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            eta(i,j) = wave_eta(p,pgc,xg,yg);
		}
    }
    pgc->gcsl_start4(p,eta,50);
    


    count=0;
    // only relaxation-zone cells, SLICELOOP4 order (cached geometry)
    if(!rz4_built) relaxzone4_build(p);
    for(size_t rzq=0; rzq<rz4_i.size(); ++rzq)
    {
        i = rz4_i[rzq];
        j = rz4_j[rzq];
        xg = rz4_xg[rzq];
        yg = rz4_yg[rzq];
        dg = rz4_dg[rzq];
        
        z = eta(i,j);
		
		// Wave Generation
		if(p->B98==2 && f_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            { 
            Fifsfval[count] = wave_fi(p,pgc,xg,yg,z);
            
            ++count;
            }
		}
    }
    
}
    
