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

#include"ghostcell.h"
#include"fdm.h"
#include"lexer.h"

void ghostcell::nse2(lexer *p, fdm *a, field &f, int gcv)
{
    double nx,ny,nz,dnorm;
    double xp, yp, zp;
    double lsv;
    
    
    VAIRLOOP
    f(i,j,k)=0.0;
    

    VAIRLOOP
    {
    lsv = 0.5*(a->phi(i,j,k)+a->phi(i,j+1,k));
    
        if(lsv>-4.1*(1.0/3.0)*(p->DXP[IP] + p->DYN[JP] + p->DZP[KP]))
        {         
        nx = (a->phi(i+1,j,k)-a->phi(i,j,k))/p->DXP[IP];
        ny = (a->phi(i,j+1,k)-a->phi(i,j,k))/p->DYP[JP];
        nz = (a->phi(i,j,k+1)-a->phi(i,j,k))/p->DZP[KP];

        dnorm = sqrt(nx*nx + ny*ny + nz*nz);
        
        nx/=dnorm;
        ny/=dnorm;
        nz/=dnorm;
        
        xp = p->pos2_x() + nx*(1.0*fabs(lsv)+0.0*p->DXP[IP]);
        yp = p->pos2_y() + ny*(1.0*fabs(lsv)+0.0*p->DYN[JP]);
        zp = p->pos2_z() + nz*(1.0*fabs(lsv)+0.0*p->DZP[KP]);
        
        // chk bounds
        f(i,j,k) = p->ccipol2_a(f, xp, yp, zp);
        }
    
    }
    
    
}
