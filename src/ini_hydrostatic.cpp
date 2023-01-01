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

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

void initialize::hydrostatic(lexer *p, fdm *a, ghostcell *pgc)
{
    double maxh=0.0;
    LOOP
    maxh=MAX(maxh, p->pos_z());

    maxh=pgc->globalmax(maxh);
	
    if(p->F30>0)
    maxh=p->phimean;
	
	if(p->I12==1 && (p->I30==0||p->B90==0))
    LOOP
    a->press(i,j,k) = (p->phimean-p->pos_z())*a->ro(i,j,k)*fabs(p->W22);

	if(p->I12==2 && (p->I30==0||p->B90==0))
    LOOP
    a->press(i,j,k) = a->phi(i,j,k)*a->ro(i,j,k)*fabs(p->W22);
	
	if(p->I12==3 && (p->I30==0||p->B90==0))
    LOOP
    a->press(i,j,k) = (maxh-p->pos_z())*a->ro(i,j,k)*fabs(p->W22);

    if(p->I56==1)
    LOOP
    {
    if(a->phi(i,j,k)<0.0)
    a->press(i,j,k)=0.0;
    }
	
    LOOP
    a->press(i,j,k)+=p->I55;
	
	if(p->I12==2 && p->I30==0)
    GC4LOOP
	{
	i = p->gcb4[n][0];
	j = p->gcb4[n][1];
	k = p->gcb4[n][2];
	
    a->press(i,j,k) = a->phi(i,j,k)*a->ro(i,j,k)*fabs(p->W22) + p->I55;
	}

    
    
    
	

}





