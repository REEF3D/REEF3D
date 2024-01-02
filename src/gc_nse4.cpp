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

void ghostcell::nse4(lexer *p, fdm *a, field &f, int gcv)
{
    double nx,ny,nz,dnorm;
    double xp, yp, zp;
    double lsv,val;
    
    if(gcv>=40 && gcv<=45)
    {
    AIRLOOP
    f(i,j,k)=0.0; 
    /*
    AIRLOOP
    {
    lsv = a->phi(i,j,k);
    
        if(lsv>-3.1*(1.0/3.0)*(p->DXP[IP] + p->DYP[JP] + p->DZP[KP]))
        {
         nx = (a->phi(i+1,j,k)-a->phi(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]);
         ny = (a->phi(i,j+1,k)-a->phi(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP]);
         nz = (a->phi(i,j,k+1)-a->phi(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]);  

        dnorm = sqrt(nx*nx + ny*ny + nz*nz);
        
        nx/=dnorm;
        ny/=dnorm;
        nz/=dnorm;
        
        xp = p->pos_x() + nx*(1.0*fabs(lsv)+0.5*p->DXM);
        yp = p->pos_y() + ny*(1.0*fabs(lsv)+0.5*p->DXM);
        zp = p->pos_z() + nz*(1.0*fabs(lsv)+0.5*p->DXM);

        
        val = 0.0;p->ccipol4(f, xp, yp, zp);
        
        f(i,j,k) =  -val;
        }
    }*/
    
    
    /*
    AIRLOOP
    {
    val = 0.0;//(a->eta(i,j))*9.81;
     
    if(p->flag4[Im1JK]>0)
    f(i,j,k) =  val*(1.0-fabs(a->phi(i,j,k))/p->DXM) - f(i-1,j,k)*fabs(a->phi(i,j,k))/p->DXM;
    
    if(p->flag4[Ip1JK]>0)
    f(i,j,k) =  val*(1.0-fabs(a->phi(i,j,k))/p->DXM) - f(i+1,j,k)*fabs(a->phi(i,j,k))/p->DXM;
    
    if(p->flag4[IJKm1]>0)
    f(i,j,k) =  val*(1.0-fabs(a->phi(i,j,k))/p->DXM) - f(i,j,k-1)*fabs(a->phi(i,j,k))/p->DXM;
    
    if(p->flag4[IJKp1]>0)
    f(i,j,k) =  val*(1.0-fabs(a->phi(i,j,k))/p->DXM) - f(i,j,k+1)*fabs(a->phi(i,j,k))/p->DXM;
    }*/
        
        /*
    AIRLOOP
    {
    if(p->flag4[Im1JK]>0)
    f(i,j,k) =  f(i-1,j,k);
    
    if(p->flag4[Ip1JK]>0)
    f(i,j,k) =  f(i+1,j,k);
    
    if(p->flag4[IJKm1]>0)
    f(i,j,k) =  f(i,j,k-1);
    
    if(p->flag4[IJKp1]>0)
    f(i,j,k) =  f(i,j,k+1);
    }*/
    
    /*
    AIRLOOP
    {
    double val = (a->eta(i,j))*fabs(p->W22);
    
    if(p->flag4[Im1JK]>0)
    f(i,j,k) =  val;
    
    if(p->flag4[Ip1JK]>0)
    f(i,j,k) =  val;
    
    if(p->flag4[IJKm1]>0)
    f(i,j,k) =  val;
    
    if(p->flag4[IJKp1]>0)
    f(i,j,k) =  val;
    }
    */

    }
}




