/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2021 Tobias Martin

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
--------------------------------------------------------------------*/

#include"FSI_strip.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void fsi_strip::interpolate_vel(lexer* p, fdm* a, ghostcell* pgc, field& uvel, field& vvel, field& wvel)
{
    int ii, jj, kk;
    double dx, dy, dz, dist, D;

    for (int eI = 0; eI < Ne; eI++)
    {
        lagrangeVel[eI] = Eigen::MatrixXd::Zero(3,lagrangePoints[eI].cols());   
    
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            const Eigen::Vector3d& coordI = lagrangePoints[eI].col(pI);

            if 
            (
                coordI(0) >= xstart[p->mpirank] && coordI(0) < xend[p->mpirank] &&
                coordI(1) >= ystart[p->mpirank] && coordI(1) < yend[p->mpirank] &&
                coordI(2) >= zstart[p->mpirank] && coordI(2) < zend[p->mpirank]
            )
            {
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
                            dist = (p->XN[i_it + 1 + marge] - coordI(0))/dx;
                            D = kernel_roma(dist);
                            dist = (p->YP[j_it + marge] - coordI(1))/dy;
                            D *= kernel_roma(dist);
                            dist = (p->ZP[k_it + marge] - coordI(2))/dz;
                            D *= kernel_roma(dist);
                            
                            lagrangeVel[eI](0,pI) += uvel(i_it,j_it,k_it)*D;

                            dist = (p->XP[i_it + marge] - coordI(0))/dx;
                            D = kernel_roma(dist);
                            dist = (p->YN[j_it + 1 + marge] - coordI(1))/dy;
                            D *= kernel_roma(dist);
                            dist = (p->ZP[k_it + marge] - coordI(2))/dz;
                            D *= kernel_roma(dist);
                                
                            lagrangeVel[eI](1,pI) += vvel(i_it,j_it,k_it)*D;
                            
                            dist = (p->XP[i_it + marge] - coordI(0))/dx;
                            D = kernel_roma(dist);
                            dist = (p->YP[j_it + marge] - coordI(1))/dy;
                            D *= kernel_roma(dist);
                            dist = (p->ZN[k_it + 1 + marge] - coordI(2))/dz;
                            D *= kernel_roma(dist);
                             
                            lagrangeVel[eI](2,pI) += wvel(i_it,j_it,k_it)*D;
                        }
                    }
                }
            }
        }
    }
    
    for (int eI = 0; eI < Ne; eI++)
    {
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            lagrangeVel[eI].col(pI) << pgc->globalsum(lagrangeVel[eI](0,pI)), pgc->globalsum(lagrangeVel[eI](1,pI)), pgc->globalsum(lagrangeVel[eI](2,pI));
        }
    }
}
    
double fsi_strip::kernel_roma(const double& dist)
{
    double D = 0.0;

    if (fabs(dist) <= 0.5)
    {
        D = 1.0/3.0*(1.0 + sqrt(-3*dist*dist + 1));
    }
    else if (fabs(dist) <= 1.5)
    {    
        D = 1.0/6.0*(5.0 - 3.0*fabs(dist) - sqrt(-3*(1 - fabs(dist))*(1 - fabs(dist)) + 1));
    }
    
    return D;
}

void fsi_strip::coupling_vel()
{
    getTransVel(xdot_el);
    getRotPos(q_el);
    getRotVel(qdot_el);

    Eigen::Vector3d omega_el;

    for (int eI = 0; eI < Ne; eI++)
    {
        omega_el = getOmega(q_el.col(eI),qdot_el.col(eI));

        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            lagrangeVelCoup[eI].col(pI) = (xdot_el.col(eI+1) + xdot_el.col(eI))/2.0 + omega_el.cross(Xil[eI].col(pI));
        }
    }
}

void fsi_strip::coupling_force(lexer *p, double alpha)
{
    for (int eI = 0; eI < Ne; eI++)
    {
        lagrangeForceCoup[eI] = (lagrangeVelCoup[eI] - lagrangeVel[eI])/(alpha*p->dt);
    }
}

void fsi_strip::distribute_forces(lexer *p, fdm *a, ghostcell *pgc, field1& fx, field2& fy, field3& fz)
{
    int ii, jj, kk;
    double dx, dy, dz, dV, dist, D;

    for (int eI = 0; eI < Ne; eI++)
    {
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            const Eigen::Vector3d& coordI = lagrangePoints[eI].col(pI);
            const Eigen::Vector3d& forceI = lagrangeForceCoup[eI].col(pI);
            const double& areaI = lagrangeArea[eI](pI);

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
                        dV = areaI*dx_body;

                        dist = (p->XN[i_it + 1 + marge] - coordI(0))/dx;
                        D = kernel_roma(dist);
                        dist = (p->YP[j_it + marge] - coordI(1))/dy;
                        D *= kernel_roma(dist);
                        dist = (p->ZP[k_it + marge] - coordI(2))/dz;
                        D *= kernel_roma(dist);
                        
                        fx(i_it,j_it,k_it) += forceI(0)*D*dV/(dx*dy*dz);
      

                        dist = (p->XP[i_it + marge] - coordI(0))/dx;
                        D = kernel_roma(dist);
                        dist = (p->YN[j_it + 1 + marge] - coordI(1))/dy;
                        D *= kernel_roma(dist);
                        dist = (p->ZP[k_it + marge] - coordI(2))/dz;
                        D *= kernel_roma(dist);
                        
                        fy(i_it,j_it,k_it) += forceI(1)*D*dV/(dx*dy*dz);
                

                        dist = (p->XP[i_it + marge] - coordI(0))/dx;
                        D = kernel_roma(dist);
                        dist = (p->YP[j_it + marge] - coordI(1))/dy;
                        D *= kernel_roma(dist);
                        dist = (p->ZN[k_it + 1 + marge] - coordI(2))/dz;
                        D *= kernel_roma(dist);
                        
                        fz(i_it,j_it,k_it) += forceI(2)*D*dV/(dx*dy*dz);
                    }
                }
            }
        }
    }
    
    pgc->start1(p,fx,10);
    pgc->start2(p,fy,11);
    pgc->start3(p,fz,12);

    ULOOP
    a->test(i,j,k) = fx(i,j,k); 
    
    pgc->start4(p,a->test,10);
}
