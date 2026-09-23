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
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::nhflow_precalc_relax(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    p->wavetime = p->simtime;
    
    if(!gen_built) genzone4_build(p,pgc);
    
// ETA (SLICELOOP4 order)
    if(p->B98==2)
    for(size_t q=0; q<gen_i.size(); ++q)
    {
        i = gen_i[q];
        j = gen_j[q];
        
        PSLICECHECK4
        eta(i,j) = wave_eta_c(p,pgc,int(q));
    }
    pgc->gcsl_start4(p,eta,50);
    
// U, V, W in one pass over the generation-zone cells (LOOP order), with the
// three velocities of a cell from one evaluation of the phase and depth
// functions. uval/vval/wval share the count sequence of the former three loops.
    count=0;
    
    if(p->B98==2)
    for(size_t q=0; q<gen_i.size(); ++q)
    {
        i = gen_i[q];
        j = gen_j[q];
        
        const double hval = eta(i,j) + d->depth(i,j);
        
        KLOOP
        PCHECK
        {
            z=p->ZSP[IJK]-p->phimean;
            
            double uw,vw,ww;
            wave_uvw_c(p,pgc,int(q),z,uw,vw,ww);
            
            uval[count] = uw + p->Ui;
            UHval[count] = hval*uval[count];
            
            if(p->j_dir==1 && v_switch==1)
            {
            vval[count] = vw;
            VHval[count] = hval*vval[count];
            }
            
            if(w_switch==1)
            {
            wval[count] = ww;
            WHval[count] = hval*wval[count];
            }
            
            ++count;
        }
    }
}
    
