/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"net_interface.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<sys/stat.h>

void net_interface::dlm_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int nNet)
{ 
    // Distribute net forces on surrounding cells
    int ii, jj, kk;
    double dist, D, dx, dy, dz;
    double test = 0.0;

    const EigenMat& lagrangePoints = ppnet->getLagrangePoints();
    const EigenMat& lagrangeForces = ppnet->getLagrangeForces();

    LOOP
    kernel_x(i,j,k) = 0.0;

    LOOP
    kernel_y(i,j,k) = 0.0;

    LOOP
    kernel_z(i,j,k) = 0.0;
    

    for (int pI = 0; pI < lagrangePoints.size(); pI++)
    {
        const Eigen::Vector3d& forcesI = lagrangeForces[pI];
        
        if (forcesI.norm()!=0.0)
        {
            const Eigen::Vector3d& coordI = lagrangePoints[pI];

            ii = p->posc_i(coordI(0));
            jj = p->posc_j(coordI(1));
            kk = p->posc_k(coordI(2));

            dx = p->DXN[ii + marge];
            dy = p->DYN[jj + marge];
            dz = p->DZN[kk + marge];

            for (int i_it = ii - 2; i_it <= ii + 2; i_it++)
            {
                for (int j_it = jj - 2; j_it <= jj + 2; j_it++)
                {
                    for (int k_it = kk - 2; k_it <= kk + 2; k_it++)
                    {
                        i = i_it;
                        j = j_it;
                        k = k_it;
                        
                        dist = sqrt(pow(p->XP[IP] - coordI(0), 2.0))/dx;
                        D = kernel_peskin(dist);
                        dist = sqrt(pow(p->YP[JP] - coordI(1), 2.0))/dy;
                        D *= kernel_peskin(dist);
                        dist = sqrt(pow(p->ZSP[IJK] - coordI(2), 2.0))/dz;
                        D *= kernel_peskin(dist);
                        
                        test = forcesI(0)*D/(dx*dy*dz);
                        d->Fext[IJK] -= test;

                        dist = sqrt(pow(p->XP[IP] - coordI(0), 2.0))/dx;
                        D = kernel_peskin(dist);
                        dist = sqrt(pow(p->YP[JP] - coordI(1), 2.0))/dy;
                        D *= kernel_peskin(dist);
                        dist = sqrt(pow(p->ZSP[IJK] - coordI(2), 2.0))/dz;
                        D *= kernel_peskin(dist);
                            
                        test = forcesI(1)*D/(dx*dy*dz);
                        d->Gext[IJK] -= test;
                        
                        dist = sqrt(pow(p->XP[IP] - coordI(0), 2.0))/dx;
                        D = kernel_peskin(dist);
                        dist = sqrt(pow(p->YP[JP] - coordI(1), 2.0))/dy;
                        D *= kernel_peskin(dist);
                        dist = sqrt(pow(p->ZSP[IJK] - coordI(2), 2.0))/dz;
                        D *= kernel_peskin(dist);
                         
                        test = forcesI(2)*D/(dx*dy*dz);
                        d->Hext[IJK] -= test;
                        
                    }
                }
            }
        }
    } 
    
    pgc->start1(p,d->Fext,10); 
    pgc->start2(p,d->Gext,11); 
    pgc->start3(p,d->Hext,12); 
}