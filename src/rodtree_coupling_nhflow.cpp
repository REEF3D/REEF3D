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
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include<cmath>
#include<algorithm>

namespace
{
    // Vertical linear interpolation in one sigma column, with the column
    // z-coordinates rebuilt from ZP*WL+bed (the ZSP halo is not reliable at
    // subdomain faces, which biases lexer::ccipol4V there).
    double column_value(lexer *p, fdm_nhf *d, double *f, slice &WL, int ic, int jc, double z)
    {
        const int marge = increment::marge;
        int &i = increment::i, &j = increment::j, &k = increment::k;
        i = ic; j = jc;
        const double wl = WL(i,j), zb = d->bed(i,j);
        if(p->knoz<2 || wl<1.0e-12)
        {
            k = 0;
            return f[IJK];
        }
        double zc0 = zb + p->ZP[0+marge]*wl;
        if(z<=zc0) {k = 0; return f[IJK];}
        for(int kk=0; kk<p->knoz-1; ++kk)
        {
            double za = zb + p->ZP[kk+marge]*wl, zbb = zb + p->ZP[kk+1+marge]*wl;
            if(z<=zbb)
            {
                double w = (z-za)/(zbb-za);
                i = ic; j = jc; k = kk;   double fa = f[IJK];
                k = kk+1;                 double fb = f[IJK];
                return (1.0-w)*fa + w*fb;
            }
        }
        k = p->knoz-1;
        return f[IJK];
    }

    double ipol_nhflow(lexer *p, fdm_nhf *d, double *f, slice &WL, double x, double y, double z)
    {
        const int marge = increment::marge;
        int i0 = p->posf_i(x);
        double wa = (p->XP[i0+1+marge]-x)/p->DXP[i0+marge];
        wa = std::max(0.0,std::min(1.0,wa));
        if(p->j_dir==0)
        return wa*column_value(p,d,f,WL,i0,0,z) + (1.0-wa)*column_value(p,d,f,WL,i0+1,0,z);
        int j0 = p->posf_j(y);
        double wb = (p->YP[j0+1+marge]-y)/p->DYP[j0+marge];
        wb = std::max(0.0,std::min(1.0,wb));
        return wa*wb*column_value(p,d,f,WL,i0,j0,z) + (1.0-wa)*wb*column_value(p,d,f,WL,i0+1,j0,z)
             + wa*(1.0-wb)*column_value(p,d,f,WL,i0,j0+1,z) + (1.0-wa)*(1.0-wb)*column_value(p,d,f,WL,i0+1,j0+1,z);
    }
}

void rodtree_coupling::start_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha,
                                    double *UH, double *VH, double *WH, slice &WL, bool finalize)
{
    // ------------------------------------------------------------------
    // 1. sample fluid velocity at the Lagrangian points (owner rank only)
    // ------------------------------------------------------------------
    // halo of the stage velocity must be current for interpolation at
    // points next to a subdomain boundary
    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);

    std::fill(buf.begin(),buf.end(),0.0);

    for(int e=0; e<rt.nelem(); ++e)
    for(int q=0; q<nq[e]; ++q)
    {
        Eigen::Vector3d x = rt.point(e,q,nq[e]);

        if(x(0)<p->originx || x(0)>=p->endx)
        continue;
        if(p->j_dir==1 && (x(1)<p->originy || x(1)>=p->endy))
        continue;

        i = std::max(0,std::min(p->posc_i(x(0)),p->knox-1));
        j = (p->j_dir==1) ? std::max(0,std::min(p->posc_j(x(1)),p->knoy-1)) : 0;

        if(p->wet[IJ]==0)
        continue;

        const double zb = p->ccslipol4(d->bed,x(0),x(1));
        const double zs = zb + p->ccslipol4(WL,x(0),x(1));

        if(x(2)<zb || x(2)>zs)
        continue;

        double *b = &buf[4*(first[e]+q)];
        b[0] = ipol_nhflow(p,d,d->U,WL,x(0),x(1),x(2));
        b[1] = (p->j_dir==1) ? ipol_nhflow(p,d,d->V,WL,x(0),x(1),x(2)) : 0.0;
        b[2] = ipol_nhflow(p,d,d->W,WL,x(0),x(1),x(2));
        b[3] = 1.0;
    }

    reduce_samples(pgc);
    apply_samples();
    rt.compute_hydro();

    // ------------------------------------------------------------------
    // 2. spread the reaction onto the sigma grid (acceleration)
    //    horizontal: 4-point Peskin kernel; vertical: kernel renormalised
    //    over the wet column (no loss at bed / free surface)
    // ------------------------------------------------------------------
    double Dzk[7];

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
        const Eigen::Vector3d Fp0 = F/(rho*double(nin));        // per point, per unit density

        for(int q=0; q<nq[e]; ++q)
        {
            if(buf[4*(first[e]+q)+3]<0.5)
            continue;

            Eigen::Vector3d x = rt.point(e,q,nq[e]);

            int ii = p->posc_i(x(0));
            int jj = (p->j_dir==1) ? p->posc_j(x(1)) : 0;
            int ic = std::max(0,std::min(ii,p->knox-1));
            int jc = std::max(0,std::min(jj,p->knoy-1));

            const double dx = p->DXN[ic+marge];
            const double dy = p->DYN[jc+marge];

            // point-implicit fluid-side drag (uniform-grid kernel self weight)
            i = ic; j = jc;
            double dzc = dx;
            {
                int kc = std::max(0,std::min(p->posc_sig(ic,jc,x(2)),p->knoz-1));
                i = ic; j = jc;
                if(WL(i,j)>1.0e-10) dzc = p->DZN[kc+marge]*WL(i,j);
            }
            double S = self_weight(std::fmod(std::fabs(x(0)-p->XP[ic+marge])/dx,1.0))*self_weight(0.0);
            if(p->j_dir==1)
            S *= self_weight(std::fmod(std::fabs(x(1)-p->YP[jc+marge])/dy,1.0));
            const double relax = 1.0/(1.0 + alpha*p->dt*cpt*S/(dx*dy*dzc*rho));
            const Eigen::Vector3d Fp = relax*Fp0;

            const int is = std::max(ii-2,0), ie = std::min(ii+2,p->knox-1);
            const int js = (p->j_dir==1) ? std::max(jj-2,0) : 0;
            const int je = (p->j_dir==1) ? std::min(jj+2,p->knoy-1) : 0;

            for(int ia=is; ia<=ie; ++ia)
            {
                i = ia;
                const double Dx = kernel((p->XP[IP]-x(0))/dx);
                if(Dx==0.0)
                continue;

                for(int ja=js; ja<=je; ++ja)
                {
                    i = ia; j = ja;
                    const double Dy = (p->j_dir==0) ? 1.0 : kernel((p->YP[JP]-x(1))/dy);

                    if(Dy==0.0 || p->wet[IJ]==0 || WL(i,j)<1.0e-10)
                    continue;

                    const int kk = p->posc_sig(ia,ja,x(2));      // overwrites i,j
                    i = ia; j = ja;

                    const int kc = std::max(0,std::min(kk,p->knoz-1));
                    const double dz = p->DZN[kc+marge]*WL(i,j);
                    const int ks = std::max(kk-3,0), ke = std::min(kk+3,p->knoz-1);

                    if(ke<ks || dz<1.0e-20)
                    continue;

                    double sum = 0.0;
                    for(k=ks; k<=ke; ++k)
                    {
                        Dzk[k-ks] = kernel((p->ZSP[IJK]-x(2))/dz);
                        sum += Dzk[k-ks]*p->DZN[KP]*WL(i,j);
                    }

                    if(sum<1.0e-20)
                    continue;

                    const double fac = (Dx/dx)*((p->j_dir==0) ? 1.0/dy : Dy/dy)/sum;

                    for(k=ks; k<=ke; ++k)
                    {
                        const double w = alpha*p->dt*fac*Dzk[k-ks];

                        d->U[IJK] += w*Fp(0);
                        UH[IJK]   += w*Fp(0)*WL(i,j);

                        if(p->j_dir==1)
                        {
                        d->V[IJK] += w*Fp(1);
                        VH[IJK]   += w*Fp(1)*WL(i,j);
                        }

                        d->W[IJK] += w*Fp(2);
                        WH[IJK]   += w*Fp(2)*WL(i,j);
                    }
                }
            }
        }
    }

    // ghost cells of U,V,W,UH,VH,WH are updated by nhflow_forcing::forcing

    // ------------------------------------------------------------------
    // 3. advance the structure once per time step
    // ------------------------------------------------------------------
    finish_step(p,pgc,finalize);
}
