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
#include"fdm_nhf.h"
#include"lexer.h"
#include"ghostcell.h"

void iowave::wavegen_precalc_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    starttime=pgc->timer();
    
    // time of this step for wave_prestep (HDC reads its time levels there)
    p->wavetime = p->simtime;
    wave_prestep(p,pgc);
    
    // piston/flap wavemakers: 2nd-order set-up (B113, B119) with the depth at the paddle
    nhflow_paddle_setup(p,d,pgc);
    
    // tidal / current background (B 510-514)
    nhflow_bg_update(p,d,pgc);
    
    if(p->B89==0)
    {
        if(p->B98==2)
        nhflow_precalc_relax(p,d,pgc);
                
        if(p->B98==3 || p->B98==4)
        nhflow_precalc_dirichlet(p,d,pgc);
    }
    
    if(p->B89==1)
    {
        if(p->B98==2)
        {
        nhflow_wavegen_precalc_decomp_time(p,pgc);
        nhflow_wavegen_precalc_decomp_relax(p,d,pgc);
        }
            
        if(p->B98==3 || p->B98==4)
        {
        nhflow_wavegen_precalc_decomp_time(p,pgc);
        nhflow_wavegen_precalc_decomp_dirichlet(p,pgc);
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
    
    pgc->gcsync();
}

void iowave::wavegen_precalc_ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    p->wavetime = p->simtime;
    wave_prestep(p,pgc);
    
    nhflow_paddle_setup(p,d,pgc);
    
    nhflow_bg_update(p,d,pgc);
    
    if(p->B98==2)
    nhflow_precalc_relax_ini(p,d,pgc);
        
    if(p->B98==3 || p->B98==4)
    nhflow_precalc_dirichlet_ini(p,d,pgc);
}

void iowave::nhflow_paddle_setup(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // the paddle sits at the inflow boundary: use the still-water depth there
    // (e.g. a wavemaker pit), unless the wave theory depth is given with B 94
    if(paddle_setup==1 || !(p->B92==21 || p->B92==22 || p->B92==23))
    return;
    
    double depth=0.0;
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
    depth = MAX(depth,d->depth(i,j));
    }
    
    depth = pgc->globalmax(depth);
    
    if(p->B94==1)
    depth = p->B94_wdt;
    
    // depth not available yet: try again at the next call
    if(depth<=0.0)
    return;
    
    wave_paddle_setup(p,pgc,depth);
    paddle_setup=1;
}
