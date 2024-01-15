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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<math.h>

void particle_f::seed_ini(lexer* p, fdm* a, ghostcell* pgc)
{
    // ini
    LOOP
    {
        active_box(i,j,k) = 0.0;
        active_topo(i,j,k) = 0.0;
    }
    
    // Box
    size_t cellcount=0;
    for(int qn=0;qn<p->Q110;++qn)
        LOOP
            if(p->XN[IP]>=p->Q110_xs[qn] && p->XN[IP]<p->Q110_xe[qn]
            && p->YN[JP]>=p->Q110_ys[qn] && p->YN[JP]<p->Q110_ye[qn]
            && p->ZN[KP]>=p->Q110_zs[qn] && p->ZN[KP]<p->Q110_ze[qn])
            {
                active_box(i,j,k) = 1.0;
                ++cellcount;
            }

    // Topo
    PLAINLOOP
        if((abs(a->topo(i,j,k))<(p->DZN[KP]*ceil(p->Q102)))&&(a->topo(i,j,k)<=0.25*p->DZN[KP])) //find better comparison to fix numerical drifts
        {
            active_topo(i,j,k) = 1.0;
            if(1!=active_box(i,j,k))
                cellcount++;
        }

    // guess particle demand
    if(p->Q24>0)
        ppcell = p->Q24;
    else
        ppcell = 0;
    
    partnum = cellcount * ppcell;
}

void particle_f::seed(lexer* p, fdm* a, ghostcell* pgc)
{
    if(p->Q110>0)
        posseed_box(p,a,pgc);
	if(p->Q101>0)
        posseed_topo(p,a,pgc);
}


void particle_f::posseed_box(lexer* p, fdm* a, ghostcell* pgc)
{
    if(p->Q29>0)
        srand(p->Q29);

    if(p->Q29==0)
        srand((unsigned)time(0)*(p->mpirank+1));
	
    double x,y,z;
    int flag;

    LOOP
        if(active_box(i,j,k)>0.0)
            for(int qn=0;qn<ppcell;++qn)
                {
                    if(PP.size+1>0.9*PP.capacity)
                        PP.reserve();
                    
                    x = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
                    y = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
                    z = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;

                    // pos[pactive][RADIUS] = p->Q31/2*double(rand() % int(drand/2) + int(drand/2))/drand;

                    PP.add(x,y,z,1);
                }
}


void particle_f::posseed_topo(lexer* p, fdm* a, ghostcell* pgc)
{
    if(p->Q29>0)
        srand(p->Q29);

    if(p->Q29==0)
        srand((unsigned)time(0)*(p->mpirank+1));

    double tolerance = 5e-18;
    double x,y,z,ipolTopo,ipolSolid;
    int flag;

    PLAINLOOP
        if(active_topo(i,j,k)>0.0)
            {
                if(PP.size+ppcell>0.9*PP.capacity)
                    PP.reserve();
                for(int qn=0;qn<ppcell;++qn)
                {

                    x = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
                    y = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
                    z = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;

                    // tempPos[tempActive][RADIUS] = p->Q31/2*double(rand() % int(drand/2) + int(drand/2))/drand;

                    ipolTopo = p->ccipol4_b(a->topo,x,y,z);
                    ipolSolid = p->ccipol4_b(a->solid,x,y,z);

                    if (!(ipolTopo>tolerance||ipolTopo<-p->Q102*p->DZN[KP]||ipolSolid<0))
                        PP.add(x,y,z,1);
                }
            }
}

void particle_f::posseed_suspended(lexer* p, fdm* a, ghostcell* pgc)
{
    if(p->Q29>0)
        srand(p->Q29);

    if(p->Q29==0)
        srand((unsigned)time(0)*(p->mpirank+1));
    
    double x,y,z;
    for(int n=0;n<p->gcin_count;n++)
        if(p->gcin[n][3]>0)
        {
            i=p->gcin[n][0];
            j=p->gcin[n][1];
            k=p->gcin[n][2];
            if(a->topo(i,j,k)>=0.0)
            {
                if(PP.size+ppcell>0.9*PP.capacity)
                    PP.reserve();
                for(int qn=0;qn<ppcell;++qn)
                {
                    x = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
                    y = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
                    z = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;
                    PP.add(x,y,z,1);
                }
            }
        }
}
