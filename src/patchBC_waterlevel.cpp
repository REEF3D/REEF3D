/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#include"fdm.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC::patchBC_waterlevel(lexer *p, fdm *a, ghostcell *pgc, field &phi)
{
    // hydrograph interpolation
    // waterlevel
    for(qq=0;qq<obj_count;++qq)
    if(patch[qq]->hydroFSF_flag==1)
    {
    patch[qq]->waterlevel = patchBC_hydrograph_FSF_ipol(p,pgc,qq,patch[qq]->ID);
    }
        
    // waterlevel
    for(qq=0;qq<obj_count;++qq)
    if(patch[qq]->waterlevel_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i=patch[qq]->gcb[n][0];
    j=patch[qq]->gcb[n][1];
    k=patch[qq]->gcb[n][2];
    
        if(patch[qq]->gcb[n][3]==1)
        {
        a->phi(i-1,j,k)=patch[qq]->waterlevel-p->pos_z();
        a->phi(i-2,j,k)=patch[qq]->waterlevel-p->pos_z();
        a->phi(i-3,j,k)=patch[qq]->waterlevel-p->pos_z();
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        a->phi(i,j+1,k)=patch[qq]->waterlevel-p->pos_z();
        a->phi(i,j+2,k)=patch[qq]->waterlevel-p->pos_z();
        a->phi(i,j+3,k)=patch[qq]->waterlevel-p->pos_z();
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        a->phi(i,j-1,k)=patch[qq]->waterlevel-p->pos_z();
        a->phi(i,j-2,k)=patch[qq]->waterlevel-p->pos_z();
        a->phi(i,j-3,k)=patch[qq]->waterlevel-p->pos_z();
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        a->phi(i+1,j,k)=patch[qq]->waterlevel-p->pos_z();
        a->phi(i+2,j,k)=patch[qq]->waterlevel-p->pos_z();
        a->phi(i+3,j,k)=patch[qq]->waterlevel-p->pos_z();
        }
        
    
    }
} 

void patchBC::patchBC_waterlevel2D(lexer*, fdm2D*, ghostcell*, slice&)
{
}