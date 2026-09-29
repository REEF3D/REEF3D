/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"net_interface.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"net.h"
#include<sys/stat.h>
#include <Eigen/Dense>

void net_interface::dlm_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int nNet)
{ 
    // Distribute net forces (per unit fluid density) on surrounding cells.
    //
    // NHFLOW specifics:
    //  - collocated sigma grid: kernel evaluated at cell centres XP, YP, ZSP
    //  - Fext enters the conservative momentum equation d(UH)/dt (see irhs),
    //    so the spread acceleration is multiplied with the local water depth WL
    //  - the vertical kernel is renormalised over the wet water column, so no
    //    force is lost where the stencil is cut by the bed or the free surface
    //  - all ranks hold all Lagrangian forces (globalsum in coupling_dlm_nhflow),
    //    so each rank spreads onto its own interior cells only; ghost cells
    //    are filled by start4V
    int ii, jj, kk, kc, isave, jsave;
    int is, ie, js, je, ks, ke;
    double dx, dy, dz, Dx, Dy, sum, fac, w;
    double Dzk[7];

    const EigenMat& lagrangePoints = pnet[nNet]->getLagrangePoints();
    const EigenMat& lagrangeForces = pnet[nNet]->getLagrangeForces();

    for (int pI = 0; pI < lagrangePoints.size(); pI++)
    {
        const Eigen::Vector3d& forcesI = lagrangeForces[pI];
        
        if (forcesI.norm()==0.0)
        continue;
        
        const Eigen::Vector3d& coordI = lagrangePoints[pI];

        ii = p->posc_i(coordI(0));
        jj = p->posc_j(coordI(1));

        dx = p->DXN[MAX(MIN(ii,p->knox-1),0) + marge];
        dy = p->DYN[MAX(MIN(jj,p->knoy-1),0) + marge];
        
        // interior cells of this subdomain inside the kernel support
        is = MAX(ii-2,0);
        ie = MIN(ii+2,p->knox-1);
        js = MAX(jj-2,0);
        je = MIN(jj+2,p->knoy-1);
        
        if(p->j_dir==0)
        js = je = 0;

        for(i=is; i<=ie; ++i)
        {
            Dx = kernel_peskin(fabs(p->XP[IP] - coordI(0))/dx);
            
            if(Dx==0.0)
            continue;
            
            for(j=js; j<=je; ++j)
            {
                // 2D: the force acts over the full width of the single cell row
                Dy = (p->j_dir==0) ? 1.0 : kernel_peskin(fabs(p->YP[JP] - coordI(1))/dy);
                
                if(Dy==0.0 || p->wet[IJ]==0)
                continue;
                
                // vertical position in this column (posc_sig overwrites i,j,k)
                isave = i;
                jsave = j;
                kk = p->posc_sig(i,j,coordI(2));
                i = isave;
                j = jsave;
                
                kc = MAX(MIN(kk,p->knoz-1),0);
                dz = p->DZN[kc + marge]*d->WL(i,j);
                
                ks = MAX(kk-3,0);
                ke = MIN(kk+3,p->knoz-1);
                
                if(ke<ks || dz<1.0e-20)
                continue;
                
                sum = 0.0;
                for(k=ks; k<=ke; ++k)
                {
                    Dzk[k-ks] = kernel_peskin(fabs(p->ZSP[IJK] - coordI(2))/dz);
                    sum += Dzk[k-ks]*p->DZN[KP]*d->WL(i,j);
                }
                
                if(sum<1.0e-20)
                continue;
                
                // Dx/dx*Dy/dy*Dz/sum integrates to one over the cells; times WL for d(UH)/dt
                fac = (Dx/dx)*(p->j_dir==0 ? 1.0/dy : Dy/dy)*d->WL(i,j)/sum;
                
                for(k=ks; k<=ke; ++k)
                {
                    w = fac*Dzk[k-ks];
                    
                    d->Fext[IJK] -= forcesI(0)*w;
                    d->Gext[IJK] -= forcesI(1)*w;
                    d->Hext[IJK] -= forcesI(2)*w;
                }
            }
        }
    } 
    
    pgc->start4V(p,d->Fext,10); 
    pgc->start4V(p,d->Gext,11); 
    pgc->start4V(p,d->Hext,12); 
}
