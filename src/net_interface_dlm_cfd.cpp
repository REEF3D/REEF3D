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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"net_interface.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"net.h"
#include<sys/stat.h>
#include <Eigen/Dense>

void net_interface::dlm_cfd(lexer *p, fdm *a, ghostcell *pgc, int nNet)
{ 
    // Distribute net forces on surrounding cells
    int ii, jj, kk;
    double dist, D, dx, dy, dz;
    double test = 0.0;

    const EigenMat& lagrangePoints = pnet[nNet]->getLagrangePoints();
    const EigenMat& lagrangeForces = pnet[nNet]->getLagrangeForces();

    ULOOP
    kernel_x(i,j,k) = 0.0;

    VLOOP
    kernel_y(i,j,k) = 0.0;

    WLOOP
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
                        dist = sqrt(pow(p->XN[i_it + 1 + marge] - coordI(0), 2.0))/dx;
                        D = kernel_peskin(dist);
                        dist = sqrt(pow(p->YP[j_it + marge] - coordI(1), 2.0))/dy;
                        D *= kernel_peskin(dist);
                        dist = sqrt(pow(p->ZP[k_it + marge] - coordI(2), 2.0))/dz;
                        D *= kernel_peskin(dist);
                        
                        test = forcesI(0)*D/(dx*dy*dz);
                        a->Fext(i_it,j_it,k_it) -= test;

                        dist = sqrt(pow(p->XP[i_it + marge] - coordI(0), 2.0))/dx;
                        D = kernel_peskin(dist);
                        dist = sqrt(pow(p->YN[j_it + 1 + marge] - coordI(1), 2.0))/dy;
                        D *= kernel_peskin(dist);
                        dist = sqrt(pow(p->ZP[k_it + marge] - coordI(2), 2.0))/dz;
                        D *= kernel_peskin(dist);
                            
                        test = forcesI(1)*D/(dx*dy*dz);
                        a->Gext(i_it,j_it,k_it) -= test;
                        
                        dist = sqrt(pow(p->XP[i_it + marge] - coordI(0), 2.0))/dx;
                        D = kernel_peskin(dist);
                        dist = sqrt(pow(p->YP[j_it + marge] - coordI(1), 2.0))/dy;
                        D *= kernel_peskin(dist);
                        dist = sqrt(pow(p->ZN[k_it + 1 + marge] - coordI(2), 2.0))/dz;
                        D *= kernel_peskin(dist);
                         
                        test = forcesI(2)*D/(dx*dy*dz);
                        a->Hext(i_it,j_it,k_it) -= test;
                        
                    }
                }
            }
        }
    } 
    
    pgc->start1(p,a->Fext,10); 
    pgc->start2(p,a->Gext,11); 
    pgc->start3(p,a->Hext,12); 
}