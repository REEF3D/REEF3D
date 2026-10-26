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
Authors: Hans Bihs
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::solid_forcing_flag_update(lexer *p, fdm *a)
{   
    // Update DF
    set_DF(p, a);
    
    // a patch of the mesh refinement (no partition neighbours): the cells around it in the domain take the
    // flag from their own fields (filled from the next patch or the coarser level), as the direct forcing
    // flags in gcdf_update; they kept the 1 of the patch setup, so the implicit diffusion (D 22 1) coupled
    // the first fluid cell of a patch to a bed or solid cell around it as to a fluid cell with zero velocity:
    // a no-slip wall the grids without refinement do not have (it slowed the flow over a patch whose lower
    // side lies on the bed by 6-7 %)
    auto margin = [p](int ii, int jj, int kk)
    {
        if(ii>=0 && ii<p->knox && jj>=0 && jj<p->knoy && kk>=0 && kk<p->knoz)
        return false;
        return ii+p->origin_i>=0 && ii+p->origin_i<p->gknox && jj+p->origin_j>=0 && jj+p->origin_j<p->gknoy
            && kk+p->origin_k>=0 && kk+p->origin_k<p->gknoz;
    };
    auto mloop = [&](auto &&f)
    {
        if(!do_comms)
        for(i=p->imin; i<p->imin+p->imax; ++i)
        for(j=p->jmin; j<p->jmin+p->jmax; ++j)
        for(k=p->kmin; k<p->kmin+p->kmax; ++k)
        if(margin(i,j,k) && p->flag4[IJK]>OBJ_FLAG)
        f();
    };

    mloop([&]()
    {
        p->DF[IJK] = ((p->solidread>0 && a->solid(i,j,k)<0.0) || (p->toporead>0 && a->topo(i,j,k)<0.0)
                  || (p->X10>0 && a->fb(i,j,k)<0.0)) ? -1 : 1;
    });
    
    // 1
    ULOOP
    p->DF1[IJK]=p->DF[IJK];
    
    ULOOP
    if(p->DF[IJK]>0 && p->DF[Ip1JK]<0)
    p->DF1[IJK]=-1;
    
    // 2
    VLOOP
    p->DF2[IJK]=p->DF[IJK];
    
    VLOOP
    if(p->DF[IJK]>0 && p->DF[IJp1K]<0)
    p->DF2[IJK]=-1;
    
    // 3
    WLOOP
    p->DF3[IJK]=p->DF[IJK];
    
    WLOOP
    if(p->DF[IJK]>0 && p->DF[IJKp1]<0)
    p->DF3[IJK]=-1;
    
    
    // the cells around a patch
    mloop([&]()
    {
        p->DF1[IJK] = p->DF2[IJK] = p->DF3[IJK] = p->DF[IJK];
        if(p->DF[IJK]>0 && i+1<p->imin+p->imax && p->DF[Ip1JK]<0)
        p->DF1[IJK] = -1;
        if(p->DF[IJK]>0 && j+1<p->jmin+p->jmax && p->DF[IJp1K]<0)
        p->DF2[IJK] = -1;
        if(p->DF[IJK]>0 && k+1<p->kmin+p->kmax && p->DF[IJKp1]<0)
        p->DF3[IJK] = -1;
    });

    // MPI exchange (not for the patches of the mesh refinement, whose work is rank local)
    if(do_comms)
    {
    gcparaxintV(p, p->DF1);
    gcparaxintV(p, p->DF2);
    gcparaxintV(p, p->DF3);
    }

}
