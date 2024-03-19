/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"patchBC.h"
#include"lexer.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC::patchBC_gcb_count(lexer *p, ghostcell *pgc)
{
// count gcbs
    for(qn=0;qn<p->B440;++qn)
    {
            int count=0;
        
            istart = p->posc_i(p->B440_xs[qn]);
            iend = p->posc_i(p->B440_xe[qn]);
            
            jstart = p->posc_j(p->B440_ys[qn]);
            jend = p->posc_j(p->B440_ye[qn]);
            
            for(n=0;n<p->gcb4_count;++n)
            {
            i=p->gcb4[n][0];
            j=p->gcb4[n][1];
            k=p->gcb4[n][2];
            
                if(i>=istart && i<iend && j>=jstart && j<jend && p->gcb4[n][3]==p->B440_face[qn] && (p->gcb4[n][4]==21||p->gcb4[n][4]==22))
                {
                ++count;
                }
            }
        
        
        for(qq=0;qq<obj_count;++qq)
        {
        if(patch[qq]->ID==p->B440_ID[qn])
        patch[qq]->gcb_count += count;
        }
    
    }
    
    for(qn=0;qn<p->B441;++qn)
    {
            int count=0;
        
            istart = p->posc_i(p->B441_xs[qn]);
            iend = p->posc_i(p->B441_xe[qn]);
            
            jstart = p->posc_j(p->B441_ys[qn]);
            jend = p->posc_j(p->B441_ye[qn]);
            
            kstart = p->posc_k(p->B441_zs[qn]);
            kend = p->posc_k(p->B441_ze[qn]);
            
            
            for(n=0;n<p->gcb4_count;++n)
            {
            i=p->gcb4[n][0];
            j=p->gcb4[n][1];
            k=p->gcb4[n][2];
            
                if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend && p->gcb4[n][3]==p->B441_face[qn] && (p->gcb4[n][4]==21||p->gcb4[n][4]==22))
                {
                ++count;
                }
            }
    
        
        for(qq=0;qq<obj_count;++qq)
        {
        if(patch[qq]->ID==p->B441_ID[qn])
        patch[qq]->gcb_count += count;
        }
    
    }
    
    
    // circle
    double r;
    for(qn=0;qn<p->B442;++qn)
    {
            int count=0;
          
            
            for(n=0;n<p->gcb4_count;++n)
            {
            i=p->gcb4[n][0];
            j=p->gcb4[n][1];
            k=p->gcb4[n][2];
            
                // x-dir
                if(p->B442_face[qn]==1 || p->B442_face[qn]==4)
                {
                    r = sqrt(pow(p->YP[JP]-p->B442_ym[qn],2.0)+pow(p->ZP[KP]-p->B442_zm[qn],2.0));
                
                    if(r<=p->B442_r[qn] && p->pos_x()>p->B442_xm[qn]-p->DXP[IP] && p->pos_x()<=p->B442_xm[qn]+p->DXP[IP] && p->gcb4[n][3]==p->B442_face[qn] && (p->gcb4[n][4]==21||p->gcb4[n][4]==22))
                    {
                    ++count;
                    }
                }
                
                // y-dir
                if(p->B442_face[qn]==2 || p->B442_face[qn]==3)
                {
                    r = sqrt(pow(p->XP[IP]-p->B442_xm[qn],2.0)+pow(p->ZP[KP]-p->B442_zm[qn],2.0));
                
                    if(r<=p->B442_r[qn] && p->pos_y()>p->B442_ym[qn]-p->DYP[JP] && p->pos_y()<=p->B442_ym[qn]+p->DYP[JP] && p->gcb4[n][3]==p->B442_face[qn] && (p->gcb4[n][4]==21||p->gcb4[n][4]==22))
                    {
                    ++count;
                    }
                }
                
                // z-dir
                if(p->B442_face[qn]==5 || p->B442_face[qn]==6)
                {
                    r = sqrt(pow(p->XP[IP]-p->B442_xm[qn],2.0)+pow(p->YP[JP]-p->B442_ym[qn],2.0));
                
                    if(r<=p->B442_r[qn] && p->pos_z()>p->B442_zm[qn]-p->DZP[KP] && p->pos_z()<=p->B442_zm[qn]+p->DZP[KP] && p->gcb4[n][3]==p->B442_face[qn] && (p->gcb4[n][4]==21||p->gcb4[n][4]==22))
                    {
                    ++count;
                    }
                }
            }

        
        
        for(qq=0;qq<obj_count;++qq)
        {
        if(patch[qq]->ID==p->B442_ID[qn])
        patch[qq]->gcb_count += count;
        }
    }
    
    // allocate arrays in patch_obj
    for(q=0; q<obj_count;++q)
    patch[q]->patch_obj_gcb_generate(p,pgc);
    
    
} 

