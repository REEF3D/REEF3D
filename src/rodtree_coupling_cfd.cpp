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
#include<mpi.h>

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

    // mesh refinement: the samples of the patches replace those of level 0 (amr_sample)
    if(amr_sampled)
    {
        if(npts>0)
        MPI_Allreduce(MPI_IN_PLACE,abuf.data(),5*npts,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

        for(int m=0; m<npts; ++m)
        if(abuf[5*m+4]>0.5)
        for(int c=0; c<4; ++c)
        buf[4*m+c] = abuf[5*m+c];

        amr_sampled = false;
    }

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

        const double cpt = wgt[e]*rt.drag_slope(e)/double(nin);
        const Eigen::Vector3d Fp0 = wgt[e]*F/(rho*double(nin));
        const double sh = shk[e];
        const int Rh = stencil_h(sh);

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

            double S = self_weight_s(std::fmod(std::fabs(x(0)-p->XP[ic+marge])/dx,1.0),sh)
                      *self_weight(std::fmod(std::fabs(x(2)-p->ZP[kc+marge])/dz,1.0));
            if(p->j_dir==1)
            S *= self_weight_s(std::fmod(std::fabs(x(1)-p->YP[jc+marge])/dy,1.0),sh);
            const double relax = 1.0/(1.0 + alpha*p->dt*cpt*S/(dx*dy*dz*rho));
            const Eigen::Vector3d Fp = relax*Fp0/(dx*dy*dz);

            const int is = std::max(ii-Rh,0), ie = std::min(ii+Rh,p->knox-1);
            const int js = (p->j_dir==1) ? std::max(jj-Rh,0) : 0;
            const int je = (p->j_dir==1) ? std::min(jj+Rh,p->knoy-1) : 0;
            const int ks = std::max(kk-2,0), ke = std::min(kk+2,p->knoz-1);

            // the weight of a face of component comp
            const std::function<double(int,int,int,int)> D = [&](int comp, int ia, int ja, int ka)
            {
                const double Dx = kernel_s(((comp==0 ? p->XN[ia+1+marge] : p->XP[ia+marge])-x(0))/dx,sh);
                const double Dy = (p->j_dir==0) ? 1.0 : kernel_s(((comp==1 ? p->YN[ja+1+marge] : p->YP[ja+marge])-x(1))/dy,sh);
                const double Dz = kernel(((comp==2 ? p->ZN[ka+1+marge] : p->ZP[ka+marge])-x(2))/dz);
                return Dx*Dy*Dz;
            };

            // the part of the kernel in the solids and the bed goes to the faces in the fluid
            const double sx = fluid_scale(p,a,0,ii,jj,kk,Rh,Rh,2,D);
            const double sy = (p->j_dir==1) ? fluid_scale(p,a,1,ii,jj,kk,Rh,Rh,2,D) : 0.0;
            const double sz = fluid_scale(p,a,2,ii,jj,kk,Rh,Rh,2,D);

            for(int ia=is; ia<=ie; ++ia)
            for(int ja=js; ja<=je; ++ja)
            for(int ka=ks; ka<=ke; ++ka)
            {
                // u faces
                fx(ia,ja,ka) += sx*Fp(0)*D(0,ia,ja,ka);

                // v faces
                if(p->j_dir==1)
                fy(ia,ja,ka) += sy*Fp(1)*D(1,ia,ja,ka);

                // w faces
                fz(ia,ja,ka) += sz*Fp(2)*D(2,ia,ja,ka);
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

double rodtree_coupling::fluid_scale(lexer *p, fdm *a, int comp, int ii, int jj, int kk, int rx, int ry, int rz,
                                     const std::function<double(int,int,int,int)> &D) const
{
    if(p->solidread==0 && p->toporead==0)
    return 1.0;

    // the stencil as far as the cells (with their halo) reach
    const int i0 = std::max(ii-rx,-marge+1), i1 = std::min(ii+rx,p->knox+marge-2);
    const int j0 = (p->j_dir==1) ? std::max(jj-ry,-marge+1) : 0;
    const int j1 = (p->j_dir==1) ? std::min(jj+ry,p->knoy+marge-2) : 0;
    const int k0 = std::max(kk-rz,-marge+1), k1 = std::min(kk+rz,p->knoz+marge-2);

    const int di = comp==0 ? 1 : 0, dj = comp==1 ? 1 : 0, dk = comp==2 ? 1 : 0;

    double wall = 0.0, wfl = 0.0;
    for(int ia=i0; ia<=i1; ++ia)
    for(int ja=j0; ja<=j1; ++ja)
    for(int ka=k0; ka<=k1; ++ka)
    {
        const double w = D(comp,ia,ja,ka);
        if(w==0.0)
        continue;
        wall += w;

        const double sf = 0.5*(a->solid(ia,ja,ka) + a->solid(ia+di,ja+dj,ka+dk));
        const double tf = 0.5*(a->topo(ia,ja,ka) + a->topo(ia+di,ja+dj,ka+dk));
        if(sf>0.0 && tf>0.0)
        wfl += w;
    }

    if(wall<=0.0)
    return 1.0;
    if(wfl<0.05*wall)
    return 0.0;
    return wall/wfl;
}


// --------------------------------------------------------------------- mesh refinement
void rodtree_coupling::amr_sample_begin()
{
    abuf.assign(5*npts,0.0);
    amr_sampled = true;
}

// the points of this patch (mine): the fluid velocity of the patch, water phase only
void rodtree_coupling::amr_sample(lexer *p, fdm *a, const std::function<bool(const double*)> &mine)
{
    for(int e=0; e<rt.nelem(); ++e)
    for(int q=0; q<nq[e]; ++q)
    {
        const Eigen::Vector3d x = rt.point(e,q,nq[e]);
        const double xk[3] = {x(0),x(1),x(2)};
        if(!mine(xk))
        continue;

        double *b = &abuf[5*(first[e]+q)];
        b[4] = 1.0;

        if(p->ccipol4(a->phi,x(0),x(1),x(2))<0.0)
        continue;

        b[0] = p->ccipol1(a->u,x(0),x(1),x(2));
        b[1] = (p->j_dir==1) ? p->ccipol2(a->v,x(0),x(1),x(2)) : 0.0;
        b[2] = p->ccipol3(a->w,x(0),x(1),x(2));
        b[3] = 1.0;
    }
}

// the reaction of the last evaluation (start_cfd of level 0 in this stage) on the faces of the
// patch, as an acceleration into Fext, Gext, Hext (added to the right-hand side of the momentum)
void rodtree_coupling::amr_spread(lexer *p, fdm *a, const std::function<void(const double*,int*)> &width)
{
    // the box of the faces of the patch, with its margin
    const int ia0 = p->imin, ib0 = p->imin + p->imax - 1;
    const int ja0 = p->jmin, jb0 = p->jmin + p->jmax - 1;
    const int ka0 = p->kmin, kb0 = p->kmin + p->kmax - 1;
    const double xa = p->XN[ia0+marge], xb = p->XN[ib0+1+marge];
    const double ya = p->YN[ja0+marge], yb = p->YN[jb0+1+marge];
    const double za = p->ZN[ka0+marge], zb = p->ZN[kb0+1+marge];

    for(int e=0; e<rt.nelem(); ++e)
    {
        const Eigen::Vector3d F = rt.fluid_reaction(e);
        if(F.squaredNorm()==0.0)
        continue;

        int nin = 0;
        for(int q=0; q<nq[e]; ++q)
        if(buf[4*(first[e]+q)+3]>0.5) ++nin;
        if(nin==0)
        continue;

        const double cpt = wgt[e]*rt.drag_slope(e)/double(nin);
        const Eigen::Vector3d Fp0 = wgt[e]*F/(rho*double(nin));

        for(int q=0; q<nq[e]; ++q)
        {
            if(buf[4*(first[e]+q)+3]<0.5)
            continue;

            const Eigen::Vector3d x = rt.point(e,q,nq[e]);
            const double xk[3] = {x(0),x(1),x(2)};

            int ii = std::max(std::min(p->posc_i(x(0)),ib0),ia0);
            int jj = (p->j_dir==1) ? std::max(std::min(p->posc_j(x(1)),jb0),ja0) : 0;
            int kk = std::max(std::min(p->posc_k(x(2)),kb0),ka0);

            const double dxp = p->DXN[ii+marge], dyp = p->DYN[jj+marge], dzp = p->DZN[kk+marge];

            int f[3] = {1,1,1};
            width(xk,f);
            if(f[0]<=0 || f[1]<=0 || f[2]<=0)
            continue;

            // out of reach of the kernel
            const double rx = 2.0*f[0]*dxp, ry = 2.0*f[1]*dyp, rz = 2.0*f[2]*dzp;
            if(x(0)<xa-rx || x(0)>xb+rx || (p->j_dir==1 && (x(1)<ya-ry || x(1)>yb+ry)) || x(2)<za-rz || x(2)>zb+rz)
            continue;

            const double dx = dxp*f[0], dy = dyp*f[1], dz = dzp*f[2];

            // point-implicit drag on the cells of this patch (as start_cfd)
            double S = self_weight(std::fmod(std::fabs(x(0)-p->XP[ii+marge])/dx,1.0))
                      *self_weight(std::fmod(std::fabs(x(2)-p->ZP[kk+marge])/dz,1.0));
            if(p->j_dir==1)
            S *= self_weight(std::fmod(std::fabs(x(1)-p->YP[jj+marge])/dy,1.0));
            const double relax = 1.0/(1.0 + p->dt*cpt*S/(dx*dy*dz*rho));
            const Eigen::Vector3d Fp = relax*Fp0/(dx*dy*dz);

            const std::function<double(int,int,int,int)> D = [&](int comp, int ia, int ja, int ka)
            {
                const double Dx = kernel(((comp==0 ? p->XN[ia+1+marge] : p->XP[ia+marge])-x(0))/dx);
                const double Dy = (p->j_dir==0) ? 1.0 : kernel(((comp==1 ? p->YN[ja+1+marge] : p->YP[ja+marge])-x(1))/dy);
                const double Dz = kernel(((comp==2 ? p->ZN[ka+1+marge] : p->ZP[ka+marge])-x(2))/dz);
                return Dx*Dy*Dz;
            };

            const int Ri = 2*f[0]+1, Rj = 2*f[1]+1, Rk = 2*f[2]+1;
            const double sx = fluid_scale(p,a,0,ii,jj,kk,Ri,Rj,Rk,D);
            const double sy = (p->j_dir==1) ? fluid_scale(p,a,1,ii,jj,kk,Ri,Rj,Rk,D) : 0.0;
            const double sz = fluid_scale(p,a,2,ii,jj,kk,Ri,Rj,Rk,D);

            const int i0 = std::max(ii-Ri,ia0), i1 = std::min(ii+Ri,ib0);
            const int j0 = (p->j_dir==1) ? std::max(jj-Rj,ja0) : 0;
            const int j1 = (p->j_dir==1) ? std::min(jj+Rj,jb0) : 0;
            const int k0 = std::max(kk-Rk,ka0), k1 = std::min(kk+Rk,kb0);

            for(int ia=i0; ia<=i1; ++ia)
            for(int ja=j0; ja<=j1; ++ja)
            for(int ka=k0; ka<=k1; ++ka)
            {
                a->Fext(ia,ja,ka) += sx*Fp(0)*D(0,ia,ja,ka);
                if(p->j_dir==1)
                a->Gext(ia,ja,ka) += sy*Fp(1)*D(1,ia,ja,ka);
                a->Hext(ia,ja,ka) += sz*Fp(2)*D(2,ia,ja,ka);
            }
        }
    }
}
