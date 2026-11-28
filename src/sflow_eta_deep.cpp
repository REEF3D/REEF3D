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

#include"sflow_eta.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

void sflow_eta::wetdrydeep(lexer* p, fdm2D* b, ghostcell* pgc, slice &WL)
{
    // one sweep over the arrays for the flag of the cells (the wet flags of the neighbours, then the depth), one
    // for the faces
    const int nj = p->jmax;
    const int *const flag = p->flagslice4, *const wet = p->wet;
    int *const deep = p->deep;
    const double *const wl = WL.V;
    const bool shore = (p->A221==0);
    const bool jd = (p->j_dir==1);
    
    for(int ii=0; ii<p->knox; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;
    
        for(int q=q0; q<q0+p->knoy; ++q)
        if(flag[q]>0)
        {
        int dp = wet[q];
        
            if(shore)
            {
            if(wet[q+nj]==0 || wet[q+2*nj]==0 || wet[q-nj]==0 || wet[q-2*nj]==0)
            dp=0;
            
            if(jd)
            // (two cells, as in x: with three the solution was not symmetric in x and y where cells are dry)
            if(wet[q+1]==0 || wet[q+2]==0 || wet[q-1]==0 || wet[q-2]==0)
            dp=0;
            
            if(jd)
            if(wet[q+nj+1]==0 || wet[q+nj-1]==0 || wet[q-nj+1]==0 || wet[q-nj-1]==0)
            dp=0;
            
            if(wl[q]<=10.0*wd_criterion)
            dp=0;
            }
        
        deep[q] = dp;
        }
    }
    
    {
    const int *const flag1 = p->flagslice1, *const flag2 = p->flagslice2;
    int *const deep1 = b->deep1.V, *const deep2 = b->deep2.V;
    const int ie = p->knox - p->ulast;
    const int je = p->knoy - p->vlast;
    
    for(int ii=0; ii<p->knox; ++ii)
    {
    const int q0 = (ii - p->imin)*nj - p->jmin;
    
        if(ii<ie)
        for(int q=q0; q<q0+p->knoy; ++q)
        if(flag1[q]>0)
        deep1[q] = (deep[q]==1 && deep[q+nj]==1)?1:0;
        
        for(int q=q0; q<q0+je; ++q)
        if(flag2[q]>0)
        deep2[q] = (deep[q]==1 && deep[q+1]==1)?1:0;
    }
    }
    
    pgc->gcsl_start4Vint(p,p->deep,50);
}
