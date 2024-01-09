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

void ghostcell::nse1(lexer *p, fdm *a, field &f, int gcv)
{
    
    double nx,ny,nz,dnorm;
    double xp, yp, zp;
    double xc, yc, zc;
    double lsv;
    double psi;
    
        if(p->j_dir==0)        
        psi = 4.1*(1.0/2.0)*(p->DRM+p->DTM);
        
        if(p->j_dir==1)
        psi = 4.1*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);
    
    
    UAIRLOOP
    f(i,j,k)=0.0;
    

    UAIRLOOP
    {
    lsv = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
    
        
    
        if(lsv>-4.1*(1.0/3.0)*(p->DXN[IP] + p->DYP[JP] + p->DZP[KP]))
        {
        nx = (a->phi(i+1,j,k)-a->phi(i,j,k))/p->DXP[IP];
        ny = (0.5*(a->phi(i,j+1,k)+a->phi(i+1,j+1,k)) - 0.5*(a->phi(i,j-1,k)+a->phi(i+1,j-1,k)))/(p->DYP[JP]+p->DYP[JM1]);
        nz = (0.5*(a->phi(i,j,k+1)+a->phi(i+1,j,k+1)) - 0.5*(a->phi(i,j,k-1)+a->phi(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KM1]);

        dnorm = sqrt(nx*nx + ny*ny + nz*nz);
        
        nx/=dnorm;
        ny/=dnorm;
        nz/=dnorm;
        
        xc = p->pos1_x();
        yc = p->pos1_y();
        zc = p->pos1_z();
        
        xp = p->pos1_x() + nx*(1.0*fabs(lsv)+0.0*p->DXN[IP]);
        yp = p->pos1_y() + ny*(1.0*fabs(lsv)+0.0*p->DYP[JP]);
        zp = p->pos1_z() + nz*(1.0*fabs(lsv)+0.0*p->DZP[KP]);
        
        
        f(i,j,k) = p->ccipol1_a(f, xp, yp, zp);
        
        //if(p->mpirank==3)
        cout<<" xc: "<<xc<<" yc: "<<yc<<" zc: "<<zc<<" | "<<" xp: "<<xp<<" yp: "<<yp<<" zp: "<<zp<<" |  nx: "<<nx<<" ny: "<<ny<<" nz: "<<nz<<" lsm: "<<lsv<<" |  "<<p->ccipol1_a(f, xp, yp, zp)<<" "<<p->ccipol4_a(a->phi, xp, yp, zp)<<endl;
        }
    
    }
    
    
}
