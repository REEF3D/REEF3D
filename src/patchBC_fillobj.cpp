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

void patchBC::patchBC_fillobj(lexer *p, ghostcell *pgc)
{
    

// fill BC options
  
    // inflow/outflow
    for(qn=0;qn<p->B210;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID == p->B210_ID[qn])
        {
        patch[qq]->io_flag=p->B210_io[qn];
        }
    }
    
    // discharge
    for(qn=0;qn<p->B211;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID == p->B211_ID[qn])
        {
        patch[qq]->Q_flag=1;
        patch[qq]->Q=p->B211_Q[qn];
        
        patch[qq]->gcb_uflag=2;
        }
    }
    
    // pressure
    for(qn=0;qn<p->B212;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID == p->B212_ID[qn])
        {
        patch[qq]->pressure_flag=1;
        patch[qq]->pressure=p->B212_pressBC[qn];
        
        patch[qq]->gcb_pressflag=2;
        }
    }
    
    // waterlevel
    for(qn=0;qn<p->B213;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID == p->B213_ID[qn])
        {
        patch[qq]->waterlevel_flag=1;
        patch[qq]->waterlevel=p->B213_h[qn];
        
        patch[qq]->gcb_phiflag=2;
        }
    }
    
    // perpendicular velocity
    for(qn=0;qn<p->B214;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID == p->B214_ID[qn])
        {
        patch[qq]->Uio_flag=1;
        patch[qq]->Uio=p->B214_Uio[qn];
        
        patch[qq]->gcb_uflag=2;
        }
    }
    
    // perpendicular velocity
    for(qn=0;qn<p->B215;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID == p->B215_ID[qn])
        {
        patch[qq]->velcomp_flag=1;
        patch[qq]->U=p->B215_U[qn];
        patch[qq]->V=p->B215_V[qn];
        patch[qq]->V=p->B215_W[qn];
        
        patch[qq]->gcb_uflag=2;
        }
    }
    
    // inflow angle
    for(qn=0;qn<p->B216;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID == p->B216_ID[qn])
        {
        patch[qq]->flowangle_flag=1;
        patch[qq]->alpha=p->B216_alpha[qn];
        }
    }
    
    // inflow normals
    for(qn=0;qn<p->B217;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID == p->B217_ID[qn])
        {
        patch[qq]->velcomp_flag=1;
        patch[qq]->Nx=p->B217_Nx[qn];
        patch[qq]->Ny=p->B217_Ny[qn];
        patch[qq]->Nz=p->B217_Nz[qn];
        }
    }
    
    /*
    111 - 222
    110 - 221
    100 - 211
    101 - 212
    011 - 122
    010 -121
    001 -112
    000 -111
    */
    
// fill gcb_flags
    for(qq=0;qq<obj_count;++qq)
    {
    if(patch[qq]->gcb_uflag==2 && patch[qq]->gcb_pressflag==2 && patch[qq]->gcb_phiflag==2)
    patch[qq]->gcb_flag = 222;
    
    if(patch[qq]->gcb_uflag==2 && patch[qq]->gcb_pressflag==2 && patch[qq]->gcb_phiflag==1)
    patch[qq]->gcb_flag = 221;
    
    if(patch[qq]->gcb_uflag==2 && patch[qq]->gcb_pressflag==1 && patch[qq]->gcb_phiflag==1)
    patch[qq]->gcb_flag = 211;
    
    if(patch[qq]->gcb_uflag==2 && patch[qq]->gcb_pressflag==1 && patch[qq]->gcb_phiflag==2)
    patch[qq]->gcb_flag = 212;
    
    if(patch[qq]->gcb_uflag==1 && patch[qq]->gcb_pressflag==2 && patch[qq]->gcb_phiflag==2)
    patch[qq]->gcb_flag = 122;
    
    if(patch[qq]->gcb_uflag==1 && patch[qq]->gcb_pressflag==2 && patch[qq]->gcb_phiflag==1)
    patch[qq]->gcb_flag = 121;
    
    if(patch[qq]->gcb_uflag==1 && patch[qq]->gcb_pressflag==1 && patch[qq]->gcb_phiflag==2)
    patch[qq]->gcb_flag = 112;
    
    if(patch[qq]->gcb_uflag==1 && patch[qq]->gcb_pressflag==1 && patch[qq]->gcb_phiflag==1)
    patch[qq]->gcb_flag = 111;
    }
    
    
    
// fill gcbs
    int count=0;
    for(qn=0;qn<p->B221;++qn)
    {
        
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
                    for(qq=0;qq<obj_count;++qq)
                    if(patch[qq]->ID == p->B221_ID[qn])
                    {
                    patch[qq]->gcb[count][0]=i;
                    patch[qq]->gcb[count][1]=j;
                    patch[qq]->gcb[count][2]=k;
                    ++count;
                    
                    // convert gcb
                    p->gcb4[n][4]=patch[qq]->gcb_flag;
                    }
                }
            }
        }
    }
    
}



