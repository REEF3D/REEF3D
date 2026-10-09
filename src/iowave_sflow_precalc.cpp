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
#include"fdm2D.h"
#include"ghostcell.h"

void iowave::wavegen_2D_precalc(lexer *p, fdm2D *b, ghostcell *pgc)
{
    starttime=pgc->timer();
    p->wavetime = p->simtime;
    
    // only relaxation generation (B 98 2) uses the precalc values
    if(p->B98!=2)
    {
    p->wavecalctime+=pgc->timer()-starttime;
    return;
    }
    
    double u_val,v_val,w_val;
    double deltaz;
    
    // generation-zone columns, registered for the cached evaluation (iowave_dist.cpp); only these
    // are evaluated (before, the depth-averaged velocities were computed for every column of the
    // domain and kept only in the zone). Zone sources (B 524) per column.
    if(!gen_built) genzone4_build(p,pgc);
    
    const bool zsel = zones.has_sources();
    
    if(p->B98==2 && h_switch==1)
    for(size_t q=0; q<gen_i.size(); ++q)
    {
        i = gen_i[q];
        j = gen_j[q];
        
        if(zsel)
        select_sources(gen_src[q]);
        
        PSLICECHECK4
        eta(i,j) = wave_eta_c(p,pgc,int(q));
    }
    
    if(zsel)
    select_sources(nullptr);
    
    pgc->gcsl_start4(p,eta,50);
    
    // depth-averaged wave velocities at the cell centres (cell-centred HLL SFLOW), the count
    // sequence of the relaxation functions (SLICELOOP4 order, zone columns)
    const bool uv = p->B98==2 && (u_switch==1 || v_switch==1);
    const bool ww = p->B98==2 && w_switch==1;
    int cuv=0, cw=0;
    
    if(uv || ww)
    for(size_t q=0; q<gen_i.size(); ++q)
    {
        i = gen_i[q];
        j = gen_j[q];
        
        PSLICECHECK4
        {
        if(zsel)
        select_sources(gen_src[q]);
        
        deltaz = (eta(i,j) + p->wd - b->bed(i,j))/(double(p->B160));
        u_val=0.0;
        v_val=0.0;
        w_val=0.0;
        z=-p->wd;
        
        for(int qn=0;qn<=p->B160;++qn)
        {
        double uw, vw, wq;
        wave_uvw_c(p,pgc,int(q),z,uw,vw,wq);
        u_val += uw;
        if(p->j_dir==1)
        v_val += vw;
        w_val += wq;
        z+=deltaz;
        }
        
        u_val/=double(p->B160+1);
        v_val/=double(p->B160+1);
        w_val/=double(p->B160+1);
        
        // Boussinesq (A 220 4): velocity at the reference level z_a
        if(p->A220==4)
        {
        double uw, vw, wq;
        z = -0.53*(p->wd - b->bed(i,j)) + 0.47*eta(i,j);
        wave_uvw_c(p,pgc,int(q),z,uw,vw,wq);
        u_val = uw;
        v_val = (p->j_dir==1)?vw:0.0;
        }
        
        if(uv)
        {
        uval[cuv] = u_val;
        vval[cuv] = v_val;
        ++cuv;
        }
        
        if(ww)
        {
        wval[cw] = w_val;
        ++cw;
        }
        }
    }
    
    if(zsel)
    select_sources(nullptr);
    
    p->wavecalctime+=pgc->timer()-starttime;
}
