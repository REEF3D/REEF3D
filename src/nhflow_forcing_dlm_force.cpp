/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_forcing::dlm_forcecalc(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double alpha, double *U, double *V, double *W, slice &WL)
{
    for(int n=0; n<p->A584; ++n)
    for(int q=0; q<Np; ++q)
    if(EL_f[n][q]==1)
    {
    EL_FX[n][q] = 0.0;
    EL_FY[n][q] = 0.0;
    EL_FZ[n][q] = 0.0;        
    }

    for(int n=0; n<p->A584; ++n)
    for(int q=0; q<Np; ++q)
    if(EL_f[n][q]==1)
    {

        ii = p->posc_i(EL_X[n][q]);
        jj = p->posc_j(EL_Y[n][q]);
        kk = p->posc_sig(ii, jj, EL_Z[n][q]);

        dx = p->DXN[ii + marge];
        dy = p->DYN[jj + marge];
        dz = p->DZN[kk + marge]*WL(ii,jj);

        for (i = ii - 2; i <= ii + 2; ++i)
        for (j = jj - 2; j <= jj + 2; ++j)
        for (k = kk - 2; k <= kk + 2; ++k)
        {
            dist = (p->XN[IP] - EL_X[n][q])/dx;
            D = kernel(dist);
            //cout<<"EL_FX[n][q]: "<<EL_FX[n][q]<<" D: "<<D<<" dist: "<<dist<<endl;
            
            dist = (p->YP[JP] - EL_Y[n][q])/dy;
            D *= kernel(dist);
            //dist = (p->ZP[KP] - EL_Z[n][q])/dz;
            //D *= kernel(dist);
                        
            EL_FX[n][q] -= U[IJK]*D*EL_V[n][q]/(dx*dy*dz);
            
            

            dist = (p->XP[IP] - EL_X[n][q])/dx;
            D = kernel(dist);
            dist = (p->YN[JP] - EL_Y[n][q])/dy;
            D *= kernel(dist);
            //dist = (p->ZP[KP] - EL_Z[n][q])/dz;
            //D *= kernel(dist);
                        
            EL_FY[n][q] -= V[IJK]*D*EL_V[n][q]/(dx*dy*dz);
                

            dist = (p->XP[IP] - EL_X[n][q])/dx;
            D = kernel(dist);
            dist = (p->YP[JP] - EL_Y[n][q])/dy;
            D *= kernel(dist);
            //dist = (p->ZN[KP] - EL_Z[n][q])/dz;
            //D *= kernel(dist);
                        
            EL_FZ[n][q] -= W[IJK]*D*EL_V[n][q]/(dx*dy*dz);
        }     
    }
    
    
    for(int n=0; n<p->A584; ++n)
    for(int q=0; q<Np; ++q)
    if(EL_f[n][q]==1)
    {
    EL_FX[n][q] /= (alpha*p->dt);
    EL_FY[n][q] /= (alpha*p->dt);
    EL_FZ[n][q] /= (alpha*p->dt);        
    }
 
}