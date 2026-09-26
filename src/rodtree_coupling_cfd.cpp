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

#include"rodtree_coupling.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include<cmath>
#include<algorithm>

void rodtree_coupling::start_cfd(lexer *p, fdm *a, ghostcell *pgc, double alpha,
                                 field &u, field &v, field &w, field &fx, field &fy, field &fz, bool finalize)
{
    // ------------------------------------------------------------------
    // 1. sample fluid velocity at the Lagrangian points (owner rank only),
    //    water phase only (level set phi >= 0)
    // ------------------------------------------------------------------
    // halo of the stage velocity must be current for interpolation at
    // points next to a subdomain boundary
    pgc->start1(p,u,10);
    pgc->start2(p,v,11);
    pgc->start3(p,w,12);

    std::fill(buf.begin(),buf.end(),0.0);

    for(int e=0; e<rt.nelem(); ++e)
    for(int q=0; q<nq[e]; ++q)
    {
        Eigen::Vector3d x = rt.point(e,q,nq[e]);

        if(x(0)<p->originx || x(0)>=p->endx)
        continue;
        if(p->j_dir==1 && (x(1)<p->originy || x(1)>=p->endy))
        continue;
        if(x(2)<p->originz || x(2)>=p->endz)
        continue;

        if(p->ccipol4(a->phi,x(0),x(1),x(2))<0.0)
        continue;

        double *b = &buf[4*(first[e]+q)];
        b[0] = p->ccipol1(u,x(0),x(1),x(2));
        b[1] = (p->j_dir==1) ? p->ccipol2(v,x(0),x(1),x(2)) : 0.0;
        b[2] = p->ccipol3(w,x(0),x(1),x(2));
        b[3] = 1.0;
    }

    reduce_samples(pgc);
    apply_samples();
    rt.compute_hydro();

    // ------------------------------------------------------------------
    // 2. spread the reaction onto the staggered forcing fields
    //    (acceleration, applied by momentum_forcing as u += alpha*dt*fx)
    // ------------------------------------------------------------------
    for(int e=0; e<rt.nelem(); ++e)
    {
        Eigen::Vector3d F = rt.fluid_reaction(e);
        if(F.squaredNorm()==0.0)
        continue;

        int nin = 0;
        for(int q=0; q<nq[e]; ++q)
        if(buf[4*(first[e]+q)+3]>0.5) ++nin;
        if(nin==0)
        continue;

        const double cpt = rt.drag_slope(e)/double(nin);
        const Eigen::Vector3d Fp0 = F/(rho*double(nin));

        for(int q=0; q<nq[e]; ++q)
        {
            if(buf[4*(first[e]+q)+3]<0.5)
            continue;

            Eigen::Vector3d x = rt.point(e,q,nq[e]);

            const int ii = p->posc_i(x(0));
            const int jj = (p->j_dir==1) ? p->posc_j(x(1)) : 0;
            const int kk = p->posc_k(x(2));
            const int ic = std::max(0,std::min(ii,p->knox-1));
            const int jc = std::max(0,std::min(jj,p->knoy-1));
            const int kc = std::max(0,std::min(kk,p->knoz-1));

            const double dx = p->DXN[ic+marge];
            const double dy = p->DYN[jc+marge];
            const double dz = p->DZN[kc+marge];

            double S = self_weight(std::fmod(std::fabs(x(0)-p->XP[ic+marge])/dx,1.0))
                      *self_weight(std::fmod(std::fabs(x(2)-p->ZP[kc+marge])/dz,1.0));
            if(p->j_dir==1)
            S *= self_weight(std::fmod(std::fabs(x(1)-p->YP[jc+marge])/dy,1.0));
            const double relax = 1.0/(1.0 + alpha*p->dt*cpt*S/(dx*dy*dz*rho));
            const Eigen::Vector3d Fp = relax*Fp0/(dx*dy*dz);

            const int is = std::max(ii-2,0), ie = std::min(ii+2,p->knox-1);
            const int js = (p->j_dir==1) ? std::max(jj-2,0) : 0;
            const int je = (p->j_dir==1) ? std::min(jj+2,p->knoy-1) : 0;
            const int ks = std::max(kk-2,0), ke = std::min(kk+2,p->knoz-1);

            for(int ia=is; ia<=ie; ++ia)
            for(int ja=js; ja<=je; ++ja)
            for(int ka=ks; ka<=ke; ++ka)
            {
                const double DyP = (p->j_dir==0) ? 1.0 : kernel((p->YP[ja+marge]-x(1))/dy);
                const double DxP = kernel((p->XP[ia+marge]-x(0))/dx);
                const double DzP = kernel((p->ZP[ka+marge]-x(2))/dz);

                // u faces
                const double DxN = kernel((p->XN[ia+1+marge]-x(0))/dx);
                fx(ia,ja,ka) += Fp(0)*DxN*DyP*DzP;

                // v faces
                if(p->j_dir==1)
                {
                const double DyN = kernel((p->YN[ja+1+marge]-x(1))/dy);
                fy(ia,ja,ka) += Fp(1)*DxP*DyN*DzP;
                }

                // w faces
                const double DzN = kernel((p->ZN[ka+1+marge]-x(2))/dz);
                fz(ia,ja,ka) += Fp(2)*DxP*DyP*DzN;
            }
        }
    }

    pgc->start1(p,fx,10);
    pgc->start2(p,fy,11);
    pgc->start3(p,fz,12);

    // ------------------------------------------------------------------
    // 3. advance the structure once per time step
    // ------------------------------------------------------------------
    finish_step(p,pgc,finalize);
}
