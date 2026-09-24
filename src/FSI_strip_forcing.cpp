/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2026 Tobias Martin

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

void fsi_strip::distribute_forces(lexer *p, fdm *a, ghostcell *pgc, field& fx, field& fy, field& fz, field& eps0)
{
    // Same kernel and stencil as before, with three changes that do not alter the result:
    //  - points whose kernel footprint cannot reach this subdomain are skipped
    //    (previously every rank looped over every Lagrangian point of every strip),
    //  - separable 1D kernel weights, zero-weight stencil entries skipped,
    //  - ghost-cell exchange, eps0 reset and epsget moved to fsi_strips::forcing (once per call).
    
    int ii, jj, kk, i_it, j_it, k_it;
    double dx, dy, dz, dV, D, dxdydz;
    double eps_star;
    double kin;
    double turb_force_fac = 50.0;
    double wxF[5], wxC[5], wyF[5], wyC[5], wzF[5], wzC[5], wzT[5];

    for (int eI = 0; eI < Ne; eI++)
    {
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            const Eigen::Vector3d& coordI = lagrangePoints[eI].col(pI);
            
            if (!near_subdomain(coordI))
            continue;
            
            const Eigen::Vector3d& forceI = lagrangeForceCoup[eI].col(pI);
            const double& areaI = lagrangeArea[eI](pI);

            ii = p->posc_i(coordI(0));
            jj = p->posc_j(coordI(1));
            kk = p->posc_k(coordI(2));

            dx = p->DXN[ii + marge];
            dy = p->DYN[jj + marge];
            dz = p->DZN[kk + marge];
            
            dV = areaI*dx_body;
            dxdydz = dx*dy*dz;
            
            for (int n = 0; n < 5; n++)
            {
                wxF[n] = kernel_roma((p->XN[ii + n - 2 + 1 + marge] - coordI(0))/dx);
                wxC[n] = kernel_roma((p->XP[ii + n - 2 + marge] - coordI(0))/dx);
                wyF[n] = kernel_roma((p->YN[jj + n - 2 + 1 + marge] - coordI(1))/dy);
                wyC[n] = kernel_roma((p->YP[jj + n - 2 + marge] - coordI(1))/dy);
                wzF[n] = kernel_roma((p->ZN[kk + n - 2 + 1 + marge] - coordI(2))/dz);
                wzC[n] = kernel_roma((p->ZP[kk + n - 2 + marge] - coordI(2))/dz);
                wzT[n] = kernel_roma((p->ZN[kk + n - 2 + marge] - coordI(2))/dz);   // as before: ZN[k], not ZP[k]
            }

            for (int ni = 0; ni < 5; ni++)
            {
                i_it = ii - 2 + ni;
                
                for (int nj = 0; nj < 5; nj++)
                {
                    j_it = jj - 2 + nj;
                    
                    for (int nk = 0; nk < 5; nk++)
                    {
                        k_it = kk - 2 + nk;
                        
                        D = wxF[ni];
                        D *= wyC[nj];
                        D *= wzC[nk];
                        if (D != 0.0)
                        fx(i_it,j_it,k_it) += forceI(0)*D*dV/dxdydz;
  
                        D = wxC[ni];
                        D *= wyF[nj];
                        D *= wzC[nk];
                        if (D != 0.0)
                        fy(i_it,j_it,k_it) += forceI(1)*D*dV/dxdydz;
            
                        D = wxC[ni];
                        D *= wyC[nj];
                        D *= wzF[nk];
                        if (D != 0.0)
                        fz(i_it,j_it,k_it) += forceI(2)*D*dV/dxdydz;
                        
                        // RANS turbulence forcing
                        if(p->T10==2)
                        if(i_it>=0 && j_it>=0 && k_it>=0 && i_it<p->knox && j_it<p->knoy && k_it<p->knoz)
                        {
                        D = wxC[ni];
                        D *= wyC[nj];
                        D *= wzT[nk];
                        
                        if (D != 0.0)
                        {
                        kin = pturb->kinval(i_it,j_it,k_it);
                        eps_star = turb_force_fac*D*pow((kin>(0.0)?(kin):(0.0)),0.5) /(0.4*0.33*(dx+dy+dz)*pow(p->cmu, 0.25));
                        
                        eps0(i_it,j_it,k_it) += eps_star;
                        }
                        }
                    }
                }
            }     
        }
    }
}

bool fsi_strip::near_subdomain(const Eigen::Vector3d& coordI) const
{
    return  coordI(0) >= xlo_ext && coordI(0) <= xhi_ext &&
            coordI(1) >= ylo_ext && coordI(1) <= yhi_ext &&
            coordI(2) >= zlo_ext && coordI(2) <= zhi_ext;
}
