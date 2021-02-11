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

#include"patchBC_2D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC_2D::patchBC_gcb_count(lexer *p, ghostcell *pgc)
{
// count gcbs
    for(qn=0;qn<p->B440;++qn)
    {
            int count=0;

            istart = p->posc_i(p->B440_xs[qn]);
            iend = p->posc_i(p->B440_xe[qn]);
            
            jstart = p->posc_j(p->B440_ys[qn]);
            jend = p->posc_j(p->B440_ye[qn]);
            
            //cout<<"p->B440_xs[qn]: "<<p->B440_xs[qn]<<" p->B440_xs[qn]: "<<p->B440_xe[qn]<<" p->B440_ys[qn]: "<<p->B440_ys[qn]<<" p->B440_ye[qn]: "<<p->B440_ye[qn]<<endl;
            //cout<<"istart: "<<istart<<" iend: "<<iend<<" jstart: "<<jstart<<" jend: "<<jend<<" p->gcbsl4_count: "<<p->gcbsl4_count<<endl;
            
            for(n=0;n<p->gcbsl4_count;++n)
            {
            i=p->gcbsl4[n][0];
            j=p->gcbsl4[n][1];
            
            //cout<<p->mpirank<<" i: "<<i<<" j: "<<j<<endl;
            
                if(i>=istart && i<iend && j>=jstart && j<jend && p->gcbsl4[n][3]==p->B440_face[qn] && (p->gcbsl4[n][4]==21||p->gcbsl4[n][4]==22))
                {
                ++count;
                }
            }
            
            //cout<<p->mpirank<<" COUNT: "<<count<<endl;

        
        for(qq=0;qq<obj_count;++qq)
        {
        if(patch[qq]->ID == p->B440_ID[qn])
        {
        patch[qq]->gcb_count += count;
        
        }
    
        }
    
    }
    
    // allocate arrays in patch_obj
    for(q=0; q<obj_count;++q)
    patch[q]->patch_obj_gcb_generate(p,pgc);
    
    
} 

