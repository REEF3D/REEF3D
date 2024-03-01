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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void particle_f::advect(lexer* p, fdm* a, ghostcell* pgc,double** f,int *flag,int active)
{
    for(n=0;n<active;++n)
    if(flag[n]>0)
    {
    u1=p->dt*p->ccipol1(a->u,f[n][0],f[n][1],f[n][2]);
    coord1=f[n][0]+u1;
	
	v1=p->dt*p->ccipol2(a->v,f[n][0],f[n][1],f[n][2]);
    coord2=f[n][1]+v1;
	
	w1=p->dt*p->ccipol3(a->w,f[n][0],f[n][1],f[n][2]);
    coord3=f[n][2]+w1;
	
	
	u2=0.25*u1 + 0.25*p->dt*p->ccipol1(a->u,coord1,coord2,coord3);
    coord1=f[n][0]+u2;
	
	v2=0.25*v1 + 0.25*p->dt*p->ccipol2(a->v,coord1,coord2,coord3);
    coord2=f[n][1]+v2;
	
	w2=0.25*w1 + 0.25*p->dt*p->ccipol3(a->w,coord1,coord2,coord3);
	coord3=f[n][2]+w2;
	
	
	f[n][0] = f[n][0] + (2.0/3.0)*u2 + (2.0/3.0)*p->dt*p->ccipol1(a->u,coord1,coord2,coord3);

    f[n][1] = f[n][1] + (2.0/3.0)*v2 + (2.0/3.0)*p->dt*p->ccipol2(a->v,coord1,coord2,coord3);
	
	f[n][2] = f[n][2] + (2.0/3.0)*w2 + (2.0/3.0)*p->dt*p->ccipol3(a->w,coord1,coord2,coord3);


    //f[n][3]=phipol(p,a,f[n][0],f[n][1],f[n][2]);
    }
}


