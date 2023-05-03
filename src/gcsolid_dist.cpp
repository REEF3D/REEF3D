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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::gcsolid_gcb_dist(lexer *p, fdm *a)
{
    double nx,ny,nz,norm;

    for(n=p->gcb_fix;n<p->gcb_solid;++n)
    {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];
        
        if(p->gcb4[n][4]==22)
        {
        nx=(a->solid(i+1,j,k)-a->solid(i-1,j,k))/(2.0*p->DXM);
        ny=(a->solid(i,j+1,k)-a->solid(i,j-1,k))/(2.0*p->DXM);
        nz=(a->solid(i,j,k+1)-a->solid(i,j,k-1))/(2.0*p->DXM);

        norm=sqrt(nx*nx + ny*ny + nz*nz);

        nx/=norm;
        ny/=norm;
        nz/=norm;

        if(p->gcb4[n][3]==1)
        p->gcd4[n]=fabs(a->solid(i,j,k))/(fabs(nx)>1.0e-15?fabs(nx):1.0e-15);

        if(p->gcb4[n][3]==4)
        p->gcd4[n]=fabs(a->solid(i,j,k))/(fabs(nx)>1.0e-15?fabs(nx):1.0e-15);

        if(p->gcb4[n][3]==3)
        p->gcd4[n]=fabs(a->solid(i,j,k))/(fabs(ny)>1.0e-15?fabs(ny):1.0e-15);

        if(p->gcb4[n][3]==2)
        p->gcd4[n]=fabs(a->solid(i,j,k))/(fabs(ny)>1.0e-15?fabs(ny):1.0e-15);

        if(p->gcb4[n][3]==5)
        p->gcd4[n]=fabs(a->solid(i,j,k))/(fabs(nz)>1.0e-15?fabs(nz):1.0e-15);

        if(p->gcb4[n][3]==6)
        p->gcd4[n]=fabs(a->solid(i,j,k))/(fabs(nz)>1.0e-15?fabs(nz):1.0e-15);
        
        
        if(p->gcb4[n][3]==1)
        p->gcd4[n]=fabs(a->solid(i,j,k));

        if(p->gcb4[n][3]==4)
        p->gcd4[n]=fabs(a->solid(i,j,k));

        if(p->gcb4[n][3]==3)
        p->gcd4[n]=fabs(a->solid(i,j,k));

        if(p->gcb4[n][3]==2)
        p->gcd4[n]=fabs(a->solid(i,j,k));

        if(p->gcb4[n][3]==5)
        p->gcd4[n]=fabs(a->solid(i,j,k));

        if(p->gcb4[n][3]==6)
        p->gcd4[n]=fabs(a->solid(i,j,k));
        }
		
		p->gcd4[n] = MIN(p->gcd4[n],p->DXM);
		p->gcd4[n] = MAX(p->gcd4[n],0.1*p->DXM);
    }
}



