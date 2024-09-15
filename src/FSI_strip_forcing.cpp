/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

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

#include"FSI_strip.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"

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

void fsi_strip::distribute_forces(lexer *p, fdm *a, ghostcell *pgc, field& fx, field& fy, field& fz)
{
    int ii, jj, kk;
    double dx, dy, dz, dV, dist, D;
    double eps_star;
    double kin;
    double turb_force_fac = 50.0;
    
    pip=4;
    
    LOOP
    eps0(i,j,k) = 0.0;
    
    pgc->start4(p,eps0,30);

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
                        
                        
                        // RANS turbulence forcing
                        if(p->T10==2)
                        if(i_it>=0 && j_it>=0 && k_it>=0 && i_it<p->knox && j_it<p->knoy && k_it<p->knoz)
                        {
                        dist = (p->XP[i_it + marge] - coordI(0))/dx;
                        D = kernel_roma(dist);
                        dist = (p->YP[j_it + marge] - coordI(1))/dy;
                        D *= kernel_roma(dist);
                        dist = (p->ZN[k_it + marge] - coordI(2))/dz;
                        D *= kernel_roma(dist);
                        
                        kin = pturb->kinval(i_it,j_it,k_it);
                        eps_star = turb_force_fac*D*pow((kin>(0.0)?(kin):(0.0)),0.5) /(0.4*0.33*(dx+dy+dz)*pow(p->cmu, 0.25));
                        
                        eps0(i_it,j_it,k_it) += eps_star;
                        }
                    }
                }
            }     
        }
    }
    
    
    pip=0;
    
    if(p->T10==2)
    LOOP
    if(eps0(i,j,k)>1.0e-8)
    {
    eps_star = eps0(i,j,k);
    pturb->epsget(i,j,k,eps_star);
    }
    
    pgc->start1(p,fx,10);
    pgc->start2(p,fy,11);
    pgc->start3(p,fz,12);
    
    pgc->start4(p,a->test,10);
    
    
}