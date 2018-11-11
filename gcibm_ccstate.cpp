/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::ccstateflag(lexer* p, fdm* a)
{
 for(n=0;n<p->facetnum;++n)
 if(p->ccstate[n]==1)
 {
     i=p->facet[n][0];
     j=p->facet[n][1];
     k=p->facet[n][2];


     if(p->flag5[IJK]==-10)
     p->flag5[IJK]=-20;



     if(p->flag5[IJK]!=-10)
     p->flag5[IJK]=-30;
 }

    p->tpcellnum=0;

    LOOP
    if(p->flag5[IJK]!=-20 && p->flag5[IJK]!=-30)
    ++p->tpcellnum;


}
