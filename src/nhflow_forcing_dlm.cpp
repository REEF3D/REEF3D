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

void nhflow_forcing::dlm_forcing(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double alpha, double *U, double *V, double *W, slice &WL)
{
    /*
    LOOP
    eps0(i,j,k) = 0.0;
    
    pgc->start4(p,eps0,30);*/

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
                       
            dist = (p->XP[IP] - EL_X[n][q])/dx;
            D = kernel(dist);
            dist = (p->YP[JP] - EL_Y[n][q])/dy;
            D *= kernel(dist);
            dist = (p->ZSP[IJK] - EL_Z[n][q])/dz;
            D *= kernel(dist);
                        
            FX[IJK] += EL_FX[n][q]*D*EL_V[n][q]/(dx*dy*dz);
            
            d->test[IJK] += D;
      

            dist = (p->XP[IP] - EL_X[n][q])/dx;
            D = kernel(dist);
            dist = (p->YP[JP] - EL_Y[n][q])/dy;
            D *= kernel(dist);
            dist = (p->ZSP[IJK] - EL_Z[n][q])/dz;
            D *= kernel(dist);
                        
            FY[IJK] += EL_FY[n][q]*D*EL_V[n][q]/(dx*dy*dz);
                

            dist = (p->XP[IP] - EL_X[n][q])/dx;
            D = kernel(dist);
            dist = (p->YP[JP] - EL_Y[n][q])/dy;
            D *= kernel(dist);
            dist = (p->ZSP[IJK] - EL_Z[n][q])/dz;
            D *= kernel(dist);
                        
            FZ[IJK] += EL_FZ[n][q]*D*EL_V[n][q]/(dx*dy*dz);
                        
                        
            // RANS turbulence forcing
            /*if(p->A560==2)
            if(i_it>=0 && j_it>=0 && k_it>=0 && i_it<p->knox && j_it<p->knoy && k_it<p->knoz)
            {
            dist = (p->XP[IP] - EL_X[n][q])/dx;
            D = kernel(dist);
            dist = (p->YP[JP] - EL_Y[n][q])/dy;
            D *= kernel(dist);
            dist = (p->ZN[KP] - EL_Z[n][q])/dz;
            D *= kernel(dist);
                        
            kin = pturb->kinval(i_it,j_it,k_it);
            eps_star = turb_force_fac*D*pow((kin>(0.0)?(kin):(0.0)),0.5) /(0.4*0.33*(dx+dy+dz)*pow(p->cmu, 0.25));
                        
            eps0(i_it,j_it,k_it) += eps_star;
            }*/

        }     
    }
    

    /*
    if(p->T10==2)
    LOOP
    if(eps0(i,j,k)>1.0e-8)
    {
    eps_star = eps0(i,j,k);
    pturb->epsget(i,j,k,eps_star);
    }*/
    
}