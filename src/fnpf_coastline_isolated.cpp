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


#include"fnpf_coastline.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"slice.h"

// Wet cells are classified locally (wd - bed >= wd_criterion), so inland
// depressions below the still water level become separate wet pockets with
// their own coastline. Here only the wet area connected (through cell faces)
// to the deepest wet cell of the domain is kept; all other wet cells are set
// dry (coastline = -1, and wet = 0 when wetting-drying is on).
// Parallel flood fill: alternating sweeps with a halo swap after each sweep,
// until no cell changes on any rank.
//
// Input: coastline = +1 wet / -1 dry, as set in fnpf_coastline::start
//    L : 1 = wet and connected, 0 = not (yet) connected, -1 = dry / inactive

void fnpf_coastline::isolated(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f, int *wet)
{
    double dmax,changed;
    int sweep,iter,maxiter,nsweep,count;
    
    // deepest wet cell
    dmax = -1.0e20;
    
    SLICELOOP4
    if(f(i,j)>0.0)
    dmax = MAX(dmax, p->wd - c->bed(i,j));
    
    dmax = pgc->globalmax(dmax);
    
    if(dmax < 0.0)   // no wet cell at all
    return;
    
    SLICEBASELOOP
    L(i,j) = -1.0;
    
    SLICELOOP4
    {
    if(f(i,j)>0.0)
    L(i,j) = 0.0;
    
    if(f(i,j)>0.0 && p->wd - c->bed(i,j) >= dmax - 1.0e-10*MAX(fabs(dmax),1.0))
    L(i,j) = 1.0;
    }
    
    pgc->gcsl_start4(p,L,50);
    
    nsweep  = p->j_dir==1?4:2;
    maxiter = 100000;   // meandering channels need more outer iterations; exit on no change
    
    for(iter=0; iter<maxiter; ++iter)
    {
        changed=0.0;
        
        for(sweep=0; sweep<nsweep; ++sweep)
        {
            const int idir = (sweep%2==0)?1:-1;
            const int jdir = (sweep<2)?1:-1;
            
            for(int ii=0; ii<p->knox; ++ii)
            {
            i = idir>0?ii:p->knox-1-ii;
            
                for(int jj=0; jj<p->knoy; ++jj)
                {
                j = jdir>0?jj:p->knoy-1-jj;
                
                if(L(i,j)>-0.5 && L(i,j)<0.5)
                if(L(i-1,j)>0.5 || L(i+1,j)>0.5 || (p->j_dir==1 && (L(i,j-1)>0.5 || L(i,j+1)>0.5)))
                {
                L(i,j) = 1.0;
                changed = 1.0;
                }
                }
            }
            
            pgc->gcsl_start4(p,L,50);
        }
        
        changed = pgc->globalmax(changed);
        
        if(changed<0.5)
        break;
    }
    
    // dry out the isolated pockets
    count=0;
    
    SLICELOOP4
    if(L(i,j)>-0.5 && L(i,j)<0.5)
    {
    f(i,j) = -1.0;
    
    if(p->A343>=1)
    wet[IJ] = 0;
    
    ++count;
    }
    
    count = pgc->globalisum(count);
    
    pgc->gcsl_start4(p,f,50);
    pgc->gcsl_start4Vint(p,wet,50);
    
    if(p->mpirank==0)
    cout<<"coastline: "<<count<<" wet cells not connected to the main water body set dry ("<<iter+1<<" flood fill iterations)"<<endl<<endl;
}
