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

#include"patchBC_2D.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC_2D::patchBC_waterlevel2D(lexer *p, fdm2D *b, ghostcell *pgc, slice &eta)
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
    
        if(patch[qq]->gcb[n][3]==1)
        {
        //if(p->count<=1)
        //eta(i,j)   =  patch[qq]->waterlevel-p->wd;
        
        eta(i-1,j) =  patch[qq]->waterlevel-p->wd;
        eta(i-2,j) =  patch[qq]->waterlevel-p->wd;
        eta(i-3,j) =  patch[qq]->waterlevel-p->wd;
        
        //if(p->count<=1)
        //b->hp(i,j) =  patch[qq]->waterlevel-b->bed(i,j);
        
        b->hp(i-1,j) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hp(i-2,j) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hp(i-3,j) =  patch[qq]->waterlevel-b->bed(i,j);
        
        /*
        b->hx(i-1,j) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hx(i-2,j) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hx(i-3,j) =  patch[qq]->waterlevel-b->bed(i,j);
        
        b->hy(i-1,j) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hy(i-2,j) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hy(i-3,j) =  patch[qq]->waterlevel-b->bed(i,j);*/
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        //eta(i,j) =  patch[qq]->waterlevel-p->wd;
        eta(i,j+1) =  patch[qq]->waterlevel-p->wd;
        eta(i,j+2) =  patch[qq]->waterlevel-p->wd;
        eta(i,j+3) =  patch[qq]->waterlevel-p->wd;
        
        b->hp(i,j+1) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hp(i,j+2) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hp(i,j+3) =  patch[qq]->waterlevel-b->bed(i,j);
        
        /*
        b->hx(i,j+1) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hx(i,j+2) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hx(i,j+3) =  patch[qq]->waterlevel-b->bed(i,j);
        
        b->hy(i,j)   =  patch[qq]->waterlevel-b->bed(i,j);
        b->hy(i,j+1) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hy(i,j+2) =  patch[qq]->waterlevel-b->bed(i,j);*/
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        //eta(i,j) =  patch[qq]->waterlevel-p->wd;
        eta(i,j-1) =  patch[qq]->waterlevel-p->wd;
        eta(i,j-2) =  patch[qq]->waterlevel-p->wd;
        eta(i,j-3) =  patch[qq]->waterlevel-p->wd;
        
        b->hp(i,j-1) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hp(i,j-2) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hp(i,j-3) =  patch[qq]->waterlevel-b->bed(i,j);
        
        /*b->hx(i,j-1) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hx(i,j-2) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hx(i,j-3) =  patch[qq]->waterlevel-b->bed(i,j);
        
        b->hy(i,j-1) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hy(i,j-2) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hy(i,j-3) =  patch[qq]->waterlevel-b->bed(i,j);*/
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        //eta(i,j) =  patch[qq]->waterlevel-p->wd;
        eta(i+1,j) =  patch[qq]->waterlevel-p->wd;
        eta(i+2,j) =  patch[qq]->waterlevel-p->wd;
        eta(i+3,j) =  patch[qq]->waterlevel-p->wd;
        
        b->hp(i+1,j) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hp(i+2,j) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hp(i+3,j) =  patch[qq]->waterlevel-b->bed(i,j);
        
        /*
        b->hx(i,j)   =  patch[qq]->waterlevel-b->bed(i,j);
        b->hx(i+1,j) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hx(i+2,j) =  patch[qq]->waterlevel-b->bed(i,j);
        
        b->hy(i+1,j) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hy(i+2,j) =  patch[qq]->waterlevel-b->bed(i,j);
        b->hy(i+3,j) =  patch[qq]->waterlevel-b->bed(i,j);*/
        }
    }
}

void patchBC_2D::patchBC_waterlevel(lexer *p, fdm *a, ghostcell *pgc, field &phi)
{
} 

