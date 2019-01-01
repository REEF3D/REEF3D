/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

void geotopo::wedge(lexer* p, fdm* a, ghostcell* pgc)
{
    int qn;
    double xmin,ymin,zmin,xs,xe;
    double slope=0.0;

    // X-DIR
    for(qn=0;qn<p->G61;++qn)
    {
            zmin=MIN(p->G61_zs[qn],p->G61_ze[qn]);
            
            if(p->G61_xs[qn]<=p->G61_xe[qn])
            {
            xs = p->G61_xs[qn];
            xe = p->G61_xe[qn];
            }
            
            if(p->G61_xs[qn]>p->G61_xe[qn])
            {
            xs = p->G61_xe[qn];
            xe = p->G61_xs[qn];
            }

            slope=(p->G61_ze[qn]-p->G61_zs[qn])/(p->G61_xe[qn]-p->G61_xs[qn]);

            ALOOP
            if(p->pos_x()>=xs && p->pos_x()<xe
            && p->pos_y()>=p->G61_ys[qn] && p->pos_y()<p->G61_ye[qn]
            && p->pos_z()>=zmin && p->pos_z()<slope*(p->pos_x()-p->G61_xs[qn])+p->G61_zs[qn] )
            a->topo(i,j,k)=-1.0;
        
    }

/*
    // TYPE 2

    if(p->G61_type[qn]==3 && p->G61_ori[qn]==1)
    for(qn=0;qn<p->G61;++qn)
    {
        zmin=MIN(p->G61_zs[qn],p->G61_ze[qn]);

        slope=(p->G61_ze[qn]-p->G61_zs[qn])/(p->G61_xe[qn]-p->G61_xs[qn]);

        ALOOP
        if(p->pos_x()>=p->G61_xs[qn] && p->pos_x()<p->G61_xe[qn]
        && p->pos_y()>=p->G61_ys[qn] && p->pos_y()<p->G61_ye[qn]
        && p->pos_z()>=zmin && p->pos_z()<slope*(p->pos_x()-p->G61_xs[qn])+p->G61_zs[qn] )
        a->topo(i,j,k)=-1.0;
    }

*/


}





