/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"geotopo.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void geotopo::box(lexer* p, fdm* a, ghostcell* pgc)
{


    int qn;


    for(qn=0;qn<p->G60;++qn)
    {
        istart = conv((p->G60_xs[qn]-p->originx)/p->dx);
        iend = conv((p->G60_xe[qn]-p->originx)/p->dx);

        jstart = conv((p->G60_ys[qn]-p->originy)/p->dx);
        jend = conv((p->G60_ye[qn]-p->originy)/p->dx);

        kstart = conv((p->G60_zs[qn]-p->originz)/p->dx);
        kend = conv((p->G60_ze[qn]-p->originz)/p->dx);

        //cout<<p->mpirank<<" G60: "<<p->G60<<" . "<<istart<<" "<<iend<<" "<<jstart<<" "<<jend<<" "<<kstart<<" "<<kend<<endl;

        ALOOP
        if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
        a->topo(i,j,k)=-1.0;
    }

}




