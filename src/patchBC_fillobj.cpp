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
                    p->gcb4[n][4]=21;
                    }
                }
            }
        }
    }
    
    /*
    B210=0;        // int patchBC inflow/outflow
    B211=0;        // int patchBC discharge
    B212=0;        // int patchBC pressure BC
    B213=0;        // int patchBC waterlevel
    B214=0;        // int patchBC perpendicular velocity
    B215=0;        // int patchBC velocity components
    B216=0;        // int patchBC horizontal inflow angle
    B217=0;        // int patchBC inflow normals
    */
    
    /*
      int pressure_flag;
    double pressure;
    
    int waterlevel_flag;
    double waterlevel;
    
    int Uin_flag;
    double Uin;
    
    int velcomp_flag;
    double U,V,W;
    
    int flowangle_flag;
    double alpha;
    
    int flownormal_flag;
    double Nx,Ny,Nz;
    */
    
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
        }
    }
    
}



