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

#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void particle_pls::setradius(lexer* p, fdm* a)
{
    for(n=0;n<posactive;++n)
    if(posflag[n]>0)
    posradius(p,a,n);

    for(n=0;n<negactive;++n)
    if(negflag[n]>0)
    negradius(p,a,n);

}


void particle_pls::posradius(lexer *p, fdm *a, int qx)
{
	pos[qx][3]=phipol(p,a,pos[qx][0],pos[qx][1],pos[qx][2]);
	
	if(pos[qx][3]<rmin)
	pos[qx][4]=rmin;

	if(pos[qx][3]>rmax)
	pos[qx][4]=rmax;

	if(pos[qx][3]>=rmin && pos[qx][3]<=rmax)
	pos[qx][4]=fabs(pos[qx][3]);
}

void particle_pls::negradius(lexer *p, fdm *a, int qx)
{
	neg[qx][3]=phipol(p,a,neg[qx][0],neg[qx][1],neg[qx][2]);

	if(neg[qx][3]>-rmin)
	neg[qx][4]=rmin;

    if(neg[qx][3]<-rmax)
	neg[qx][4]=rmax;

	if(neg[qx][3]<=-rmin && neg[qx][3]>=-rmax)
	neg[qx][4]=fabs(neg[qx][3]);
}
