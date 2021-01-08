/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"patchBC.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC::patchBC_gcb_count(lexer *p, ghostcell *pgc)
{
// count gcbs
    
    for(qn=0;qn<p->B221;++qn)
    {
        int count=0;
        {
            istart = p->posc_i(p->B221_xs[qn]);
            iend = p->posc_i(p->B221_xe[qn]);
            
            jstart = p->posc_j(p->B221_ys[qn]);
            jend = p->posc_j(p->B221_ye[qn]);
            
            kstart = p->posc_k(p->B221_zs[qn]);
            kend = p->posc_k(p->B221_ze[qn]);
            
            
            for(n=0;n<p->gcb4_count;++n)
            {
            i=p->gcb4[n][0];
            j=p->gcb4[n][1];
            k=p->gcb4[n][2];
            
                if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend && p->gcb4[n][3]==p->B221_face[qn] && (p->gcb4[n][4]==21||p->gcb4[n][4]==22))
                {
                ++count;
                }
            }
        }
        
        for(qq=0;qq<obj_count;++qq)
        {
        if(patch[qq]->ID == p->B221_ID[qn])
        patch[qq]->gcb_count += count;
        }
    
    }
    
    
    // xircle
    for(qn=0;qn<p->B221;++qn)
    {
        int count=0;
        {
            istart = p->posc_i(p->B221_xs[qn]);
            iend = p->posc_i(p->B221_xe[qn]);
            
            jstart = p->posc_j(p->B221_ys[qn]);
            jend = p->posc_j(p->B221_ye[qn]);
            
            kstart = p->posc_k(p->B221_zs[qn]);
            kend = p->posc_k(p->B221_ze[qn]);
            
            
            for(n=0;n<p->gcb4_count;++n)
            {
            i=p->gcb4[n][0];
            j=p->gcb4[n][1];
            k=p->gcb4[n][2];
            
                if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend && p->gcb4[n][3]==p->B221_face[qn] && (p->gcb4[n][4]==21||p->gcb4[n][4]==22))
                {
                ++count;
                }
            }
        }
        
        for(qq=0;qq<obj_count;++qq)
        {
        if(patch[qq]->ID == p->B221_ID[qn])
        patch[qq]->gcb_count += count;
        }
    
    }
    
    
    // allocate arrays
    for(q=0; q<obj_count;++q)
    patch[q]->patch_obj_gcb_generate(p,pgc);
    
    
    
    // gcin / gcout


} 

