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

void ghostcell::nse3(lexer *p, fdm *a, field &f, int gcv)
{
    double nx,ny,nz,dnorm;
    double xp, yp, zp;
    double lsv,val;
    
    WAIRLOOP
    f(i,j,k)=0.0;

    WAIRLOOP
    {
    lsv = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
    
        if(lsv>-4.1*(1.0/3.0)*(p->DXP[IP] + p->DYP[JP] + p->DZN[KP]))
        {
        nx = (0.5*(a->phi(i+1,j,k)+a->phi(i+1,j,k+1)) - 0.5*(a->phi(i-1,j,k)+a->phi(i-1,j,k+1)))/(p->DXP[IP]+p->DXP[IM1]);
        ny = (0.5*(a->phi(i,j+1,k)+a->phi(i,j+1,k+1)) - 0.5*(a->phi(i,j-1,k)+a->phi(i,j-1,k+1)))/(p->DYP[JP]+p->DYP[JM1]);
        nz = (a->phi(i,j,k+1)-a->phi(i,j,k))/p->DZP[KP]; 

        dnorm = sqrt(nx*nx + ny*ny + nz*nz);
        
        nx/=dnorm;
        ny/=dnorm;
        nz/=dnorm;
        
        xp = p->pos3_x() + nx*(1.0*fabs(lsv)+0.0*p->DXP[IP]);
        yp = p->pos3_y() + ny*(1.0*fabs(lsv)+0.0*p->DYP[JP]);
        zp = p->pos3_z() + nz*(1.0*fabs(lsv)+0.0*p->DZN[KP]);
  
        f(i,j,k) = p->ccipol3_a(f, xp, yp, zp);  
        }

    }
}
