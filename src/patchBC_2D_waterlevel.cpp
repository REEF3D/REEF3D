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

void patchBC_2D::patchBC_waterlevel2D(lexer *p, ghostcell *pgc, slice &eta)
{
    // waterlevel
    for(qq=0;qq<obj_count;++qq)
    if(patch[qq]->waterlevel_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i=patch[qq]->gcb[n][0];
    j=patch[qq]->gcb[n][1];
    
        if(patch[qq]->gcb[n][3]==1)
        {
        eta(i-1,j) =  patch[qq]->waterlevel-p->wd;
        eta(i-2,j) =  patch[qq]->waterlevel-p->wd;
        eta(i-3,j) =  patch[qq]->waterlevel-p->wd;
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        eta(i,j)   =  patch[qq]->waterlevel-p->wd;
        eta(i,j+1) =  patch[qq]->waterlevel-p->wd;
        eta(i,j+2) =  patch[qq]->waterlevel-p->wd;
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        eta(i,j-1) =  patch[qq]->waterlevel-p->wd;
        eta(i,j-2) =  patch[qq]->waterlevel-p->wd;
        eta(i,j-3) =  patch[qq]->waterlevel-p->wd;
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        eta(i,j)   =  patch[qq]->waterlevel-p->wd;
        eta(i+1,j) =  patch[qq]->waterlevel-p->wd;
        eta(i+2,j) =  patch[qq]->waterlevel-p->wd;
        }
    }
}

void patchBC_2D::patchBC_waterlevel(lexer *p, fdm *a, ghostcell *pgc, field &phi)
{
} 

