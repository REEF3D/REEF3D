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

void fsi_strip::interpolate_vel(lexer* p, fdm* a, ghostcell* pgc, field& uvel, field& vvel, field& wvel)
{
    // Local contribution only; the global sum over all ranks is done in fsi_strips::forcing.
    // The 3D kernel is a tensor product, so the 1D weights are evaluated once per point and
    // direction (30 kernel calls instead of 1125) and zero-weight stencil entries are skipped.
    // Products and summation order are the same as before.
    
    int ii, jj, kk, i_it, j_it, k_it;
    double dx, dy, dz, D;
    double wxF[5], wxC[5], wyF[5], wyC[5], wzF[5], wzC[5];
    double su, sv, sw;

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
                
                for (int n = 0; n < 5; n++)
                {
                    wxF[n] = kernel_roma((p->XN[ii + n - 2 + 1 + marge] - coordI(0))/dx);
                    wxC[n] = kernel_roma((p->XP[ii + n - 2 + marge] - coordI(0))/dx);
                    wyF[n] = kernel_roma((p->YN[jj + n - 2 + 1 + marge] - coordI(1))/dy);
                    wyC[n] = kernel_roma((p->YP[jj + n - 2 + marge] - coordI(1))/dy);
                    wzF[n] = kernel_roma((p->ZN[kk + n - 2 + 1 + marge] - coordI(2))/dz);
                    wzC[n] = kernel_roma((p->ZP[kk + n - 2 + marge] - coordI(2))/dz);
                }
                
                su = sv = sw = 0.0;

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
                            su += uvel(i_it,j_it,k_it)*D;

                            D = wxC[ni];
                            D *= wyF[nj];
                            D *= wzC[nk];
                            if (D != 0.0)
                            sv += vvel(i_it,j_it,k_it)*D;
                            
                            D = wxC[ni];
                            D *= wyC[nj];
                            D *= wzF[nk];
                            if (D != 0.0)
                            sw += wvel(i_it,j_it,k_it)*D;
                        }
                    }
                }
                
                lagrangeVel[eI](0,pI) = su;
                lagrangeVel[eI](1,pI) = sv;
                lagrangeVel[eI](2,pI) = sw;
            }
        }
    }
}
