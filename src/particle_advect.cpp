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

void particle_f::advect(lexer* p, fdm* a, ghostcell* pgc)
{
    double coord1, coord2, coord3, u1, u2, v1, v2, w1, w2;
    PARTLOOP
        if(PP.Flag[n]>0)
        {
            u1=p->dt*p->ccipol1(a->u,PP.X[n],PP.Y[n],PP.Z[n]);
            coord1=PP.X[n]+u1;
            
            v1=p->dt*p->ccipol2(a->v,PP.X[n],PP.Y[n],PP.Z[n]);
            coord2=PP.Y[n]+v1;
            
            w1=p->dt*p->ccipol3(a->w,PP.X[n],PP.Y[n],PP.Z[n]);
            coord3=PP.Z[n]+w1;
            
            
            u2=0.25*u1 + 0.25*p->dt*p->ccipol1(a->u,coord1,coord2,coord3);
            coord1=PP.X[n]+u2;
            
            v2=0.25*v1 + 0.25*p->dt*p->ccipol2(a->v,coord1,coord2,coord3);
            coord2=PP.Y[n]+v2;
            
            w2=0.25*w1 + 0.25*p->dt*p->ccipol3(a->w,coord1,coord2,coord3);
            coord3=PP.Z[n]+w2;
            
            
            PP.X[n] += (2.0/3.0)*u2 + (2.0/3.0)*p->dt*p->ccipol1(a->u,coord1,coord2,coord3);

            PP.Y[n] += (2.0/3.0)*v2 + (2.0/3.0)*p->dt*p->ccipol2(a->v,coord1,coord2,coord3);
            
            PP.Z[n] += (2.0/3.0)*w2 + (2.0/3.0)*p->dt*p->ccipol3(a->w,coord1,coord2,coord3);
        }
}


