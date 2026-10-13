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

void iowave::fnpf_precalc_dirichlet(lexer *p, ghostcell *pgc)
{
    p->wavetime = p->simtime;
    
        count=0;
		for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        xg=xgen(p);
        yg=ygen(p);
        x1=xgen1(p);
        y2=ygen2(p);
        
        eta(i,j) = wave_eta(p,pgc,xg,yg);
        
        z = eta(i,j);
        
        time_1=time_0;
        time_0=time_n;
        time_n=p->simtime;
        Fifsfval1[count] = Fifsfval0[count];
        Fifsfval0[count] = Fifsfval[count];
        
        Fifsfval[count] = wave_u(p,pgc,xg,yg,z);
        ++count;
        }
        
        
        // Uin
        count=0;
		for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        xg=xgen(p);
        yg=ygen(p);
        x1=xgen1(p);
        y2=ygen2(p);
        
            FKLOOP
            FPCHECK
            {
            z=p->ZSN[FIJK]-p->phimean;

            Uinval[count] = wave_u(p,pgc,xg,yg,z);
            ++count;
            }
        }
}



// combined generation (B 98 6): targets of the absorbing generating edge at x- (surface elevation,
// surface velocity, inflow velocity per level), in their own arrays: the relaxation zone in front
// keeps eta / Fifsfval of the zone (fnpf_precalc_relax)
void iowave::fnpf_precalc_edge(lexer *p, ghostcell *pgc)
{
    const double wt = p->wavetime;
    
    p->wavetime = p->simtime;
    
    if((int)edge_eta.size()<p->gcslin_count)
    {
        edge_eta.assign(p->gcslin_count,0.0);
        edge_us.assign(p->gcslin_count,0.0);
    }
    
    int nlev=0;
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
        FKLOOP
        FPCHECK
        ++nlev;
    }
    
    if((int)edge_uin.size()<nlev)
    edge_uin.assign(nlev,0.0);
    
    count=0;
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
    xg=xgen(p);
    yg=ygen(p);
    
    edge_eta[n] = wave_eta(p,pgc,xg,yg);
    edge_us[n]  = wave_u(p,pgc,xg,yg,edge_eta[n]);
    
        FKLOOP
        FPCHECK
        {
        z=p->ZSN[FIJK]-p->phimean;
        
        edge_uin[count] = wave_u(p,pgc,xg,yg,z);
        ++count;
        }
    }
    
    p->wavetime = wt;
}
