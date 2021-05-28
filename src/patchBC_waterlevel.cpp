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

void patchBC::patchBC_waterlevel(lexer *p, fdm *a, ghostcell *pgc, field &phi)
{
    /*
    for(n=0;n<p->gcin_count;++n)
        {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];

        a->phi(i-1,j,k)=p->phimean-p->pos_z();
        a->phi(i-2,j,k)=p->phimean-p->pos_z();
        a->phi(i-3,j,k)=p->phimean-p->pos_z();
        }*/
} 

void patchBC::patchBC_waterlevel2D(lexer*, fdm2D*, ghostcell*, slice&)
{
}