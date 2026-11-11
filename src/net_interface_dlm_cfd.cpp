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
    dlm_cfd_spread(p,a,nNet);
    
    pgc->start1(p,a->Fext,10); 
    pgc->start2(p,a->Gext,11); 
    pgc->start3(p,a->Hext,12); 
}

// CFD mesh refinement (cfd_amr): the forces of all nets as of their last stage (the same on all
// ranks, net::getLagrangeForces) onto the grid of a patch.  Its momentum takes them from Fext,
// Gext, Hext at the faces of the patch (momentum_rk::rhs), as that of level 0 takes those of
// dlm_cfd.  The kernel of a knot is that of the cells of one level on all the grids that hold its
// force (width: the cells of that level in cells of the patch, 1, 2, 4, ... per direction; a patch
// that is covered around the knot keeps its own kernel, as level 0): the kernel of the patch
// where the patch and its coarser levels reach well around the knot, a coarser one near a patch
// edge, which the patch samples at its faces.  So the force of every knot is on the leaf faces of
// the grids once, for a net inside a patch, across a patch edge or next to one.  fmax: the widest
// kernel in cells of the patch.
void net_interface::dlm_cfd_amr(lexer *p, fdm *a, int fmax, const std::function<void(const double*,int*)> &width)
{
    if(p->X320<=0)
    return;
    
    for(int n=0; n<(int)pnet.size(); ++n)
    dlm_cfd_spread_amr(p,a,n,fmax,width);
}

void net_interface::dlm_cfd_spread_amr(lexer *p, fdm *a, int nNet, int fmax, const std::function<void(const double*,int*)> &width)
{
    const EigenMat& lagrangePoints = pnet[nNet]->getLagrangePoints();
    const EigenMat& lagrangeForces = pnet[nNet]->getLagrangeForces();

    // the box of the faces of the patch, with its margin
    const int ia = p->imin, ib = p->imin + p->imax - 1;
    const int ja = p->jmin, jb = p->jmin + p->jmax - 1;
    const int ka = p->kmin, kb = p->kmin + p->kmax - 1;
    const double xa = p->XN[ia + marge], xb = p->XN[ib + 1 + marge];
    const double ya = p->YN[ja + marge], yb = p->YN[jb + 1 + marge];
    const double za = p->ZN[ka + marge], zb = p->ZN[kb + 1 + marge];

    for (int pI = 0; pI < lagrangePoints.size(); pI++)
    {
        const Eigen::Vector3d& forcesI = lagrangeForces[pI];
        
        if (forcesI.norm()==0.0)
        continue;
        
        const Eigen::Vector3d& coordI = lagrangePoints[pI];

        int ii = p->posc_i(coordI(0));
        int jj = p->posc_j(coordI(1));
        int kk = p->posc_k(coordI(2));
        
        ii = MAX(MIN(ii,ib),ia);
        jj = MAX(MIN(jj,jb),ja);
        kk = MAX(MIN(kk,kb),ka);

        const double dxp = p->DXN[ii + marge];
        const double dyp = p->DYN[jj + marge];
        const double dzp = p->DZN[kk + marge];
        
        // out of reach of the widest kernel
        const double rx = 2.0*double(fmax)*dxp, ry = 2.0*double(fmax)*dyp, rz = 2.0*double(fmax)*dzp;
        if(coordI(0)<xa-rx || coordI(0)>xb+rx || coordI(1)<ya-ry || coordI(1)>yb+ry || coordI(2)<za-rz || coordI(2)>zb+rz)
        continue;
        
        const double xk[3] = {coordI(0),coordI(1),coordI(2)};
        int f[3] = {1,1,1};
        width(xk,f);
        if(f[0]<=0 || f[1]<=0 || f[2]<=0)
        continue;

        const double dx = dxp*double(f[0]);
        const double dy = dyp*double(f[1]);
        const double dz = dzp*double(f[2]);
        const double vol = dx*dy*dz;
        
        // the kernel per direction at the faces and at the centres
        const int i0 = MAX(ii - 2*f[0] - 1,ia), i1 = MIN(ii + 2*f[0] + 1,ib);
        const int j0 = MAX(jj - 2*f[1] - 1,ja), j1 = MIN(jj + 2*f[1] + 1,jb);
        const int k0 = MAX(kk - 2*f[2] - 1,ka), k1 = MIN(kk + 2*f[2] + 1,kb);
        if(i0>i1 || j0>j1 || k0>k1)
        continue;
        
        std::vector<double> wxn(i1-i0+1), wxp(i1-i0+1), wyn(j1-j0+1), wyp(j1-j0+1), wzn(k1-k0+1), wzp(k1-k0+1);
        for(int q=i0; q<=i1; ++q)
        {
            wxn[q-i0] = kernel_peskin(fabs(p->XN[q + 1 + marge] - coordI(0))/dx);
            wxp[q-i0] = kernel_peskin(fabs(p->XP[q + marge] - coordI(0))/dx);
        }
        for(int q=j0; q<=j1; ++q)
        {
            wyn[q-j0] = kernel_peskin(fabs(p->YN[q + 1 + marge] - coordI(1))/dy);
            wyp[q-j0] = kernel_peskin(fabs(p->YP[q + marge] - coordI(1))/dy);
        }
        for(int q=k0; q<=k1; ++q)
        {
            wzn[q-k0] = kernel_peskin(fabs(p->ZN[q + 1 + marge] - coordI(2))/dz);
            wzp[q-k0] = kernel_peskin(fabs(p->ZP[q + marge] - coordI(2))/dz);
        }

        for (int i_it = i0; i_it <= i1; i_it++)
        for (int j_it = j0; j_it <= j1; j_it++)
        for (int k_it = k0; k_it <= k1; k_it++)
        {
            a->Fext(i_it,j_it,k_it) -= forcesI(0)*(wxn[i_it-i0]*wyp[j_it-j0]*wzp[k_it-k0])/vol;
            a->Gext(i_it,j_it,k_it) -= forcesI(1)*(wxp[i_it-i0]*wyn[j_it-j0]*wzp[k_it-k0])/vol;
            a->Hext(i_it,j_it,k_it) -= forcesI(2)*(wxp[i_it-i0]*wyp[j_it-j0]*wzn[k_it-k0])/vol;
        }
    } 
}

void net_interface::dlm_cfd_spread(lexer *p, fdm *a, int nNet)
{ 
    // Distribute net forces on surrounding cells
    int ii, jj, kk;
    double dist, D, dx, dy, dz;
    double test = 0.0;

    const EigenMat& lagrangePoints = pnet[nNet]->getLagrangePoints();
    const EigenMat& lagrangeForces = pnet[nNet]->getLagrangeForces();

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
                        
                        a->Fext(i_it,j_it,k_it) -= forcesI(0)*D/(dx*dy*dz);
                        
                        
                        dist = sqrt(pow(p->XP[i_it + marge] - coordI(0), 2.0))/dx;
                        D = kernel_peskin(dist);
                        dist = sqrt(pow(p->YN[j_it + 1 + marge] - coordI(1), 2.0))/dy;
                        D *= kernel_peskin(dist);
                        dist = sqrt(pow(p->ZP[k_it + marge] - coordI(2), 2.0))/dz;
                        D *= kernel_peskin(dist);
                            
                        a->Gext(i_it,j_it,k_it) -= forcesI(1)*D/(dx*dy*dz);
                        
                        
                        dist = sqrt(pow(p->XP[i_it + marge] - coordI(0), 2.0))/dx;
                        D = kernel_peskin(dist);
                        dist = sqrt(pow(p->YP[j_it + marge] - coordI(1), 2.0))/dy;
                        D *= kernel_peskin(dist);
                        dist = sqrt(pow(p->ZN[k_it + 1 + marge] - coordI(2), 2.0))/dz;
                        D *= kernel_peskin(dist);
                         
                        a->Hext(i_it,j_it,k_it) -= forcesI(2)*D/(dx*dy*dz);
                        
                    }
                }
            }
        }
    } 
}